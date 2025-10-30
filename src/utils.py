
import pandas as pd
import numpy as np
from math import isfinite
from datetime import datetime, date

# --- minimal IRR bits (from your code, no Joakim logic) ---
def _to_date(d):
    if isinstance(d, (date, datetime)):
        return d.date() if isinstance(d, datetime) else d
    return pd.to_datetime(d, dayfirst=True, errors="coerce").date()

def _xnpv(r, cf, T):  return sum(c / (1.0 + r) ** t for c, t in zip(cf, T))
def _dxnpv(r, cf, T): return sum(-c * t * (1.0 + r) ** (-(t + 1.0)) for c, t in zip(cf, T))

def _solve_xirr(cf, T, guess=0.1, tol=1e-12):
    r = guess
    for _ in range(50):
        f, df = _xnpv(r, cf, T), _dxnpv(r, cf, T)
        if abs(df) < 1e-18: break
        r2 = r - f/df
        if not (-0.9999 < r2 < 1e12): break
        if abs(r2 - r) < tol: return r2
        r = r2
    lo, hi = -0.9999, 10.0
    flo, fhi = _xnpv(lo, cf, T), _xnpv(hi, cf, T)
    k = 0
    while flo * fhi > 0 and hi < 1e6 and k < 60:
        hi *= 2.0; fhi = _xnpv(hi, cf, T); k += 1
    if flo * fhi > 0: return float("nan")
    for _ in range(200):
        mid = 0.5*(lo+hi); fm = _xnpv(mid, cf, T)
        if abs(fm) < tol: return mid
        if flo * fm <= 0: hi, fhi = mid, fm
        else:             lo, flo = mid, fm
    return 0.5*(lo+hi)

def xirr_from_dates(cf, dates, anchor=None, basis=365.0, guess=0.1):
    dts = [_to_date(d) for d in dates]
    if anchor is None: anchor = min(dts)
    T = [(d - anchor).days / basis for d in dts]
    return _solve_xirr(cf, T, guess=guess)

def annual_and_period_return(cf, dates, basis=365.0):
    if len(cf) < 2: return float("nan"), float("nan")
    dts = [_to_date(d) for d in dates]
    r = xirr_from_dates(cf, dts, anchor=min(dts), basis=basis)
    if not isfinite(r): return float("nan"), float("nan")
    D = (max(dts) - min(dts)).days
    pr = (1.0 + r) ** (D / basis) - 1.0
    return r, pr

# --- builder + runner ---


def irr_by_groups(
    df: pd.DataFrame,
    date_from, date_to,
    irr_id: str,
    group_id: str,
    group_by=("gerdeink",),
    measures=("verdkr", "verd"),
    keep_cols=("slflokkun","mynt"),
    cpi_col="visitala",
    date_col="vidmdags",
    value_suffix="_value",
    move_suffix="_movement",
    basis=365.0,
    include_edge_movements=False
):
    d0 = pd.to_datetime(date_from).date()
    d1 = pd.to_datetime(date_to).date()

    work = df.copy()
    work[date_col] = pd.to_datetime(work[date_col], errors="coerce").dt.date
    work = work[(work[date_col] >= d0) & (work[date_col] <= d1)]
    if work.empty:
        # Build empty frame with expected columns
        base_cols = [*group_by, "start_date", "end_date"]
        meas_cols = []
        for m in measures:
            meas_cols += [f"{irr_id}_{m}", f"{irr_id}_{m}_real"]
        out_cols = base_cols + meas_cols + (list(keep_cols) if keep_cols else []) + ["group_id"]
        return pd.DataFrame(columns=out_cols)

    # Aggregate per day per group
    agg_cols = {cpi_col: "first"}
    for m in measures:
        agg_cols[m + value_suffix] = "sum"
        agg_cols[m + move_suffix] = "sum"

    perday = (
        work.groupby([*group_by, date_col], as_index=False, dropna=False)
            .agg(agg_cols)
            .sort_values(date_col)
    )

    # Keep descriptive columns (first non-null per group)
    if keep_cols:
        kept = (
            work.groupby(list(group_by), as_index=False, dropna=False)
                .agg({c: "first" for c in keep_cols})
        )
    else:
        kept = None

    rows = []
    for keys, g in perday.groupby(list(group_by), sort=False):
        g = g.sort_values(date_col)
        start_date = g[date_col].min()
        end_date   = g[date_col].max()

        mid_mask = (
            (g[date_col] > start_date) & (g[date_col] < end_date)
            if not include_edge_movements
            else (g[date_col] >= start_date) & (g[date_col] <= end_date)
        )

        res = dict(zip(group_by, keys if isinstance(keys, tuple) else (keys,)))
        res.update({"start_date": start_date, "end_date": end_date})

        for m in measures:
            val_col, mov_col = m + value_suffix, m + move_suffix

            start_v = float(g.loc[g[date_col] == start_date, val_col].sum()) + float(g.loc[g[date_col] == start_date, mov_col].sum())
            end_v   = float(g.loc[g[date_col] == end_date,   val_col].sum())
            mids    = g.loc[mid_mask, [date_col, mov_col]]

            # Sign convention (your corrected version):
            #   start value = +start_v
            #   interior movements = as-is
            #   end value = -end_v
            cf  = [start_v]
            dts = [start_date]
            for _, r in mids.iterrows():
                mv = float(r[mov_col])
                if mv != 0.0:
                    cf.append(mv)
                    dts.append(r[date_col])
            cf.append(-end_v)
            dts.append(end_date)

            if np.all(np.abs(cf) < 5):
                r_nom = np.nan
                r_real = np.nan
            else:
                # --- nominal ---
                r_nom, _ = annual_and_period_return(cf, dts, basis)

                # --- real (deflate each CF to end-date CPI) ---
                base_cpi = float(g.loc[g[date_col] == end_date, cpi_col].iloc[0])
                cpis_per_date = [float(g.loc[g[date_col] == d, cpi_col].iloc[0]) for d in dts]
                cf_real = [float(c) * (base_cpi / float(pi)) for c, pi in zip(cf, cpis_per_date)]
                r_real, _ = annual_and_period_return(cf_real, dts, basis)

            # prefix with irr_id
            res[f"{irr_id}_{m}"] = r_nom
            res[f"{irr_id}_{m}_real"] = r_real

        # attach kept cols
        if kept is not None and len(keep_cols) > 0:
            key_ser = pd.Series(keys, index=group_by) if isinstance(keys, tuple) else pd.Series({group_by[0]: keys})
            row_match = (kept[list(group_by)] == key_ser.values).all(axis=1)
            for c in keep_cols:
                res[c] = kept.loc[row_match, c].iloc[0] if row_match.any() else np.nan

        rows.append(res)

    # Column order
    base_cols = [*group_by, "start_date", "end_date"]
    meas_cols = []
    for m in measures:
        meas_cols += [f"{irr_id}_{m}", f"{irr_id}_{m}_real"]

    cols = base_cols + meas_cols + (list(keep_cols) if keep_cols else [])
    out_df = pd.DataFrame(rows)[cols].sort_values(list(group_by)).reset_index(drop=True)
    out_df["group_id"] = group_id
    # return with group_id first for clarity
    return out_df[["group_id", *cols]]