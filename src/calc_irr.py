import os, sys
import pyodbc
import pandas as pd
import numpy as np
from datetime import datetime, date
from dateutil.relativedelta import relativedelta
import argparse


print(os.getcwd())
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../db_python')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../calc_irr')))

from db_python import DataUploader
from src.utils import irr_by_groups


def parse_date(s: str | None) -> date:
    if not s:
        return date.today()
    try:
        return datetime.strptime(s, "%Y-%m-%d").date()
    except ValueError:
        return pd.to_datetime(s, dayfirst=True).date()



if __name__ == '__main__':

    # take input end_date 
    # --- CLI: take input end_date ---
    ap = argparse.ArgumentParser()
    ap.add_argument("--end-date", required=False, help="End date (YYYY-MM-DD); defaults to today.")
    ap.add_argument("--ref-date", required=False, help="The date between validFrom and validTo (YYYY-MM-DD); defaults to today.")
    ap.add_argument("--add_if_exists", required=False, help="The code checks if the end_date is already in the irr data table. If 0 then the code does not add to the table if end_date already exists ")
    args = ap.parse_args()

    end_date = parse_date(args.end_date)
    ref_date = parse_date(args.ref_date)
    add_if_exists = args.add_if_exists
    print(add_if_exists)
    print(end_date)
    y6y_ago  = end_date - relativedelta(years=1)  # 6 years ago from end_date


    # --- DB setup ---
    db_server = os.getenv("DB_SERVER")
    db_name   = os.getenv("DB_NAME", "risk")
    db_user   = os.getenv("DB_USER")
    db_pass   = os.getenv("DB_USER_PW")
    db_port   = 1433

    print("starting")

    conn_str = (
        f"DRIVER={{ODBC Driver 18 for SQL Server}};"
        f"SERVER={db_server},{db_port};"
        f"DATABASE={db_name};"
        f"UID={db_user};"
        f"PWD={{{db_pass}}};"
        f"Encrypt=yes;"
        f"TrustServerCertificate=yes;"
    )
 

    # --- Execute procedure and fetch result ---
    with pyodbc.connect(conn_str) as conn:
        query = f"exec am.fetch_irr_data '{y6y_ago}', '{end_date}', '{end_date}'"
        print(query)
        df = pd.read_sql(query, conn)

    # df.to_pickle("irr_data.pkl")
    # print(f"{len(df):,} nr. of rows to")




    current_max_date = df["vidmdags"].max().date()
    if end_date > current_max_date:
        end_date = current_max_date

    print(end_date)
    print(add_if_exists)
    dw = DataUploader(conn_str)
    current_end_dates = dw.fetch("select distinct end_date from risk.am.irr")



    if (np.isin(end_date, current_end_dates)) and ((add_if_exists == 0) or (add_if_exists is None)):
        print("end_date already in irr table. Nothing to add, exiting")
        exit()





    #df = pd.read_pickle("irr_data.pkl")

    # helpers
    def start_of_year(d): return d.replace(month=1, day=1)
    def months_ago(d, n): return (pd.Timestamp(d) - pd.DateOffset(months=n)).date()
    def years_ago(d, n):  return (pd.Timestamp(d) - pd.DateOffset(years=n)).date()


    group_by_list  = [("sjodnr","gerd"), 
                    ("sjodnr","eink"), 
                    ("sjodnr","gerdeink"), 
                    ("sjodnr", "slflokkun"),
                    ("sjodnr",),
                    ("sjodnr", "innl_erl", "mynt"),
                    ("sjodnr", "innl_erl"),
                    ("sjodnr", "slflokkun", "mynt"),
                    ("sjodnr","verdtryggt",),
                    ("sjodnr", "yflokkheiti",),
                    ("sjodnr", "fmeflheiti",)
                    ]
    measure_list   = [("verdkr", "verd"),   
                    ("verdkr", "verd"),
                    ("verdkr", "verd"),    
                    ("verdkr",),
                    ("verdkr",),
                    ("verdkr", "verd"),
                    ("verdkr",),
                    ("verdkr", "verd"),
                    ("verdkr", ),
                    ("verdkr", ),
                    ("verdkr", )]
    
    keep_cols_list = [("slflokkun","mynt","gerdeink", 'eink', "yflokkheiti", "uflokkheiti", "fmeflheiti", "innl_erl", "markskrad", "verdtryggt"),
                    ("slflokkun","mynt", "yflokkheiti", "uflokkheiti", "fmeflheiti", "innl_erl", "markskrad", "verdtryggt"),
                    ("slflokkun","mynt", "yflokkheiti", "uflokkheiti", "fmeflheiti", "innl_erl", "markskrad", "verdtryggt"), 
                    (),
                    (),
                    (),
                    (),
                    (),
                    (),
                    (),
                    ()]   
    

    assert len(group_by_list) == len(measure_list), "config lists mismatch in length"
    assert len(group_by_list) == len(keep_cols_list), "config lists mismatch in length"
    
    windows = [
        ("YTD", start_of_year(end_date), end_date),
        ("1M",  months_ago(end_date, 1), end_date),
        ("3M",  months_ago(end_date, 3), end_date),
        ("6M",  months_ago(end_date, 6), end_date),
       # ("9M",  months_ago(end_date, 9), end_date),
       # ("1Y",  years_ago(end_date, 1),  end_date),
       # ("3Y",  years_ago(end_date, 3),  end_date),
       # ("5Y",  years_ago(end_date, 5),  end_date),
    ]

    all_blocks = []
    for i, group_by in enumerate(group_by_list):
        print(group_by)
        measures  = measure_list[i]
        keep_cols = keep_cols_list[i]
        group_id_name = group_by[-1]  # e.g., "gerd" or "slflokkun"

        acc = None
        key_cols = ["group_id", *group_by, *keep_cols]  # everything identical except the IRR metrics

        for irr_id, d_from, d_to in windows:
            print(irr_id)
            res = irr_by_groups(
                df,
                date_from=d_from,
                date_to=d_to,
                irr_id=irr_id,
                group_id=group_id_name,
                group_by=group_by,
                measures=measures,
                keep_cols=keep_cols,
                include_edge_movements=False
            )

            # drop start_date for every window; keep a single end_date
            res = res.drop(columns=[c for c in ("start_date",) if c in res.columns])
            if acc is None:
                # first window initializes accumulator; keep its end_date
                acc = res.copy()
            else:
                # subsequent windows: drop their end_date to avoid dup; UNION columns on keys
                res = res.drop(columns=[c for c in ("end_date",) if c in res.columns])
                acc = acc.merge(res, on=key_cols, how="outer")

        # ensure only one end_date column and no duplicate columns overall
        acc = acc.loc[:, ~acc.columns.duplicated()]
        all_blocks.append(acc)

    # final: concatenate different group_by “blocks” vertically
    final_df = pd.concat(all_blocks, ignore_index=True)
    print(final_df.head())
    # optional: put group_id first, then the union of unique columns deterministically
    ordered = ["group_id"] + [c for c in final_df.columns if c != "group_id"]
    ordered = list(dict.fromkeys(ordered))  # de-dup while preserving order
    final_df = final_df[ordered]
    final_df["insertTime"] = datetime.now()

    
    # --- desired schema & basic coercion maps ---
    SCHEMA_ORDER = [
        "group_id", "end_date", "sjodnr", "slflokkun", "innl_erl", "mynt",
        "gerd", "gerdeink", "eink", "yflokkheiti", "uflokkheiti",
        "fmeflheiti", "markskrad", "verdtryggt",
        "YTD_verdkr", "YTD_verdkr_real",
        "1M_verdkr", "1M_verdkr_real",
        "3M_verdkr", "3M_verdkr_real",
        "6M_verdkr", "6M_verdkr_real",
        "9M_verdkr", "9M_verdkr_real",
        "1Y_verdkr", "1Y_verdkr_real",
        "3Y_verdkr", "3Y_verdkr_real",
        "5Y_verdkr", "5Y_verdkr_real",
        "YTD_verd", "YTD_verd_real",
        "1M_verd", "1M_verd_real",
        "3M_verd", "3M_verd_real",
        "6M_verd", "6M_verd_real",
        "9M_verd", "9M_verd_real",
        "1Y_verd", "1Y_verd_real",
        "3Y_verd", "3Y_verd_real",
        "5Y_verd", "5Y_verd_real",
        "insertTime",
        "active"
    ]

    FLOAT_COLS = [
        "YTD_verdkr", "YTD_verdkr_real",
        "1M_verdkr", "1M_verdkr_real",
        "3M_verdkr", "3M_verdkr_real",
        "6M_verdkr", "6M_verdkr_real",
        "9M_verdkr", "9M_verdkr_real",
        "1Y_verdkr", "1Y_verdkr_real",
        "3Y_verdkr", "3Y_verdkr_real",
        "5Y_verdkr", "5Y_verdkr_real",
        "YTD_verd", "YTD_verd_real",
        "1M_verd", "1M_verd_real",
        "3M_verd", "3M_verd_real",
        "6M_verd", "6M_verd_real",
        "9M_verd", "9M_verd_real",
        "1Y_verd", "1Y_verd_real",
        "3Y_verd", "3Y_verd_real",
        "5Y_verd", "5Y_verd_real",
    ]

    def align_to_schema(final_df: pd.DataFrame) -> pd.DataFrame:
        df = final_df.copy()

        # add missing columns as NULLs
        for c in SCHEMA_ORDER:
            if c not in df.columns:
                df[c] = pd.NA

        # optional: coerce date/float types
        if "end_date" in df.columns:
            df["end_date"] = pd.to_datetime(df["end_date"], errors="coerce").dt.date
        present_float_cols = [c for c in FLOAT_COLS if c in df.columns]
        if present_float_cols:
            df[present_float_cols] = df[present_float_cols].apply(
                pd.to_numeric, errors="coerce"
            )

        # reorder to exact schema (drop extras)
        df = df.reindex(columns=SCHEMA_ORDER)

        # If you prefer to KEEP any extra columns, append them at the end:
        # extras = [c for c in final_df.columns if c not in SCHEMA_ORDER]
        # df = df.reindex(columns=SCHEMA_ORDER + extras)

        return df
    
    def prepare_for_upload(df, end_date_value):
        # 1) ensure schema/columns exist (reuse your align_to_schema if you have it)
        # df = align_to_schema(df)

        # 2) enforce end_date (fill missing + coerce to DATE)
        if "end_date" not in df.columns:
            df["end_date"] = end_date_value
        else:
            df["end_date"] = df["end_date"].fillna(end_date_value)

        df["end_date"] = pd.to_datetime(df["end_date"], errors="coerce").dt.date

        # 3) enforce other NOT NULLs (group_id in your DDL)
        if "group_id" not in df.columns:
            raise ValueError("group_id column is missing.")
        df["group_id"] = df["group_id"].astype(object)  # keep as str-like
        if df["group_id"].isna().any():
            bad = df[df["group_id"].isna()].head(10)
            raise ValueError(f"Found NULL group_id in {len(bad)} rows (showing head):\n{bad}")

        # 4) final sanity: no NULL end_date
        if df["end_date"].isna().any():
            bad = df[df["end_date"].isna()].head(10)
            raise ValueError(f"Found NULL end_date in {len(bad)} rows (showing head):\n{bad}")

        return df

    # ---- usage before to_sql ----
    # end_date is your script argument (date or 'YYYY-MM-DD')


    final_df = prepare_for_upload(final_df, end_date)
    final_df = align_to_schema(final_df)

    final_df = final_df.loc[final_df["end_date"] == end_date]

    # upload data
    dw = DataUploader(conn_str)
    dw.insert(final_df, "am.irr")




    print(final_df.head())













