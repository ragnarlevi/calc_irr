
from datetime import datetime, date
from math import isfinite
import pandas as pd
from pathlib import Path

def _to_date(d):
    if isinstance(d, (date, datetime)):
        return d.date() if isinstance(d, datetime) else d
    for fmt in ("%d/%m/%Y", "%Y-%m-%d"):
        try:
            return datetime.strptime(str(d), fmt).date()
        except Exception:
            pass
    return pd.to_datetime(d, dayfirst=True).date()

def _xnpv(r, cf, T):
    return sum(c / (1.0 + r) ** t for c, t in zip(cf, T))

def _dxnpv(r, cf, T):
    return sum(-c * t * (1.0 + r) ** (-(t + 1.0)) for c, t in zip(cf, T))

def _solve_xirr(cf, T, guess=0.1, tol=1e-12):
    r = guess
    for _ in range(50):
        f, df = _xnpv(r, cf, T), _dxnpv(r, cf, T)
        if abs(df) < 1e-18:
            break
        r2 = r - f/df
        if not (-0.9999 < r2 < 1e12):
            break
        if abs(r2 - r) < tol:
            return r2
        r = r2
    lo, hi = -0.9999, 10.0
    flo, fhi = _xnpv(lo, cf, T), _xnpv(hi, cf, T)
    k = 0
    while flo * fhi > 0 and hi < 1e6 and k < 60:
        hi *= 2.0; fhi = _xnpv(hi, cf, T); k += 1
    if flo * fhi > 0:
        return float("nan")
    for _ in range(200):
        mid = 0.5*(lo+hi); fm = _xnpv(mid, cf, T)
        if abs(fm) < tol:
            return mid
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5*(lo+hi)

def xirr_from_dates(cf, dates, anchor=None, basis=365.0, guess=0.1):
    if anchor is None:
        anchor = min(dates)
    T = [(d - anchor).days / basis for d in dates]
    return _solve_xirr(cf, T, guess=guess)

def annual_and_period_return(cf, dates, basis=365.0):
    if len(cf) < 2:
        return float("nan"), float("nan")
    dts = [ _to_date(d) for d in dates ]
    r = xirr_from_dates(cf, dts, anchor=min(dts), basis=basis)
    if not isfinite(r):
        return float("nan"), float("nan")
    D = (max(dts) - min(dts)).days
    pr = (1.0 + r) ** (D / basis) - 1.0
    return r, pr

def xirr_joakim(cf, dates, window_end=None, move_only_if_after=True, basis=365.0, guess=0.1):
    """Move the LAST (by date) non-zero CF to window_end (if later), then compute IRR + period return."""
    dts = [_to_date(d) for d in dates]
    if not cf or not dts:
        return float("nan"), float("nan"), False, None, None

    # default window_end = group's max date
    if window_end is None:
        window_end = max(dts)
    window_end = _to_date(window_end)

    # find index of the non-zero CF with the LATEST DATE (ties: pick last occurrence on that date)
    latest_date = None
    latest_idx = None
    for i, (c, d) in enumerate(zip(cf, dts)):
        if float(c) != 0.0:
            if latest_date is None or d > latest_date or (d == latest_date and i > latest_idx):
                latest_date, latest_idx = d, i
    if latest_idx is None:
        return float("nan"), float("nan"), False, None, None

    dts_mod = list(dts)
    moved = False
    moved_from = latest_date
    moved_to = latest_date

    if (not move_only_if_after) or (window_end > latest_date):
        dts_mod[latest_idx] = window_end
        moved = (window_end != latest_date)
        moved_to = window_end

    r = xirr_from_dates(cf, dts_mod, anchor=min(dts_mod), basis=basis, guess=guess)
    if not isfinite(r):
        return float("nan"), float("nan"), moved, moved_from, moved_to
    

    D = (max(dts_mod) - min(dts_mod)).days
    pr = (1.0 + r) ** (D / basis) - 1.0
    return r, pr, moved, moved_from, moved_to



def joakim_returns(cf, dates, window_end, basis=365.0, guess=0.1):
    """
    Move the LAST non-zero CF date to window_end, then return:
      (annual_irr, period_return)
    period_return = (1 + r)^(days/basis) - 1 using the modified dates.
    """
    r, pr = xirr_joakim(cf, dates, window_end=window_end, basis=basis, guess=guess)
    return r, pr




# --- helpers for REAL (CPI-adjusted) returns ---
def _pick_base_cpi(dates, cpis, base_date=None):
    """Pick CPI for base_date; fallback to CPI at max(dates)."""
    dts = [_to_date(d) for d in dates]
    if base_date is None:
        base_date = max(dts)
    base_date = _to_date(base_date)
    # try exact match
    for d, c in zip(dts, cpis):
        if d == base_date:
            return float(c), base_date
    # fallback: CPI at latest date we have
    idx = max(range(len(dts)), key=lambda i: dts[i])
    return float(cpis[idx]), dts[idx]

def _deflate_cashflows(cf, cpis, base_cpi):
    """Return amounts in base-date prices: cf_real_i = cf_i * (base_cpi / cpi_i)."""
    return [float(c) * (base_cpi / float(pi)) for c, pi in zip(cf, cpis)]

def annual_and_period_return_real(cf, dates, cpis, basis=365.0):
    """Standard real IRR: deflate by CPI to group end-date base, then XIRR."""
    if len(cf) < 2:
        return float("nan"), float("nan"), None, None
    dts = [_to_date(d) for d in dates]
    base_cpi, base_date = _pick_base_cpi(dts, cpis, base_date=max(dts))
    cf_real = _deflate_cashflows(cf, cpis, base_cpi)
    r = xirr_from_dates(cf_real, dts, anchor=min(dts), basis=basis)
    if not isfinite(r):
        return float("nan"), float("nan"), base_cpi, base_date
    D = (max(dts) - min(dts)).days
    pr = (1.0 + r) ** (D / basis) - 1.0
    return r, pr, base_cpi, base_date


def _cpi_at_or_before(dates, cpis, target_date):
    """Return CPI for the latest observation <= target_date; fallback to earliest."""
    pairs = sorted(((_to_date(d), float(c)) for d, c in zip(dates, cpis)), key=lambda x: x[0])
    chosen = None
    for d, c in pairs:
        if d <= _to_date(target_date):
            chosen = c
        else:
            break
    return chosen if chosen is not None else pairs[0][1]

def xirr_joakim_real(cf, dates, cpis, window_end=None, move_only_if_after=True, basis=365.0, guess=0.1):
    """
    Joakim real IRR:
      • Move the LAST (by DATE) non-zero CF to window_end (if later).
      • Replace that CF's CPI with CPI(window_end).
      • Deflate all CFs by base_cpi = CPI(window_end).
      • Compute annual IRR and period return on deflated CFs.
    Returns: (annual_real, period_real, moved, moved_from, moved_to, base_cpi, base_date)
    """
    dts = [_to_date(d) for d in dates]
    if not cf or not dts:
        return float("nan"), float("nan"), False, None, None, None, None

    # choose window_end (default group end)
    if window_end is None:
        window_end = max(dts)
    window_end = _to_date(window_end)

    # find index of latest non-zero by DATE (ties -> last occurrence)
    latest_date, latest_idx = None, None
    for i, (c, d) in enumerate(zip(cf, dts)):
        if float(c) != 0.0 and (latest_date is None or d > latest_date or (d == latest_date and i > latest_idx)):
            latest_date, latest_idx = d, i
    if latest_idx is None:
        return float("nan"), float("nan"), False, None, None, None, None

    # base CPI is CPI at window_end
    base_cpi = _cpi_at_or_before(dts, cpis, window_end)
    base_date = window_end

    # prepare modified dates + CPIs
    dts_mod = list(dts)
    cpis_mod = [float(x) for x in cpis]

    moved = False
    moved_from, moved_to = latest_date, latest_date
    if (not move_only_if_after) or (window_end > latest_date):
        dts_mod[latest_idx] = window_end
        cpis_mod[latest_idx] = _cpi_at_or_before(dts, cpis, window_end)  # << replace CPI for moved CF
        moved = (window_end != latest_date)
        moved_to = window_end

    # deflate all CFs by base_cpi
    cf_real = [float(c) * (base_cpi / float(pi)) for c, pi in zip(cf, cpis_mod)]

    # compute IRR + period on real CFs
    r = xirr_from_dates(cf_real, dts_mod, anchor=min(dts_mod), basis=basis, guess=guess)
    if not isfinite(r):
        return float("nan"), float("nan"), moved, moved_from, moved_to, base_cpi, base_date

    D = (max(dts_mod) - min(dts_mod)).days
    pr = (1.0 + r) ** (D / basis) - 1.0
    return r, pr, moved, moved_from, moved_to, base_cpi, base_date







input_path = Path("./sl_uppruni.xlsx")
df = pd.read_excel(input_path)

df.columns = [str(c).strip() for c in df.columns]
for col in ("vidmdags", "dags_neys"):
    if col in df.columns:
        df[col] = pd.to_datetime(df[col], dayfirst=True, errors='coerce').dt.date

window_end_global = df["vidmdags"].max()

group_cols = ["uppruni", "sjodnr"]
results = []

for keys, g in df.groupby(group_cols, sort=False):
    # unpack keys into the same columns you grouped by
    key_dict = dict(zip(group_cols, keys))

    g = g.sort_values("vidmdags")
    dts  = g["vidmdags"].tolist()
    cfs  = g["verdkr"].astype(float).tolist()
    cpis = g["visitala"].astype(float).tolist()

    print(keys)
    print(dts)

    # Standard (nominal)
    r_std, pr_std = annual_and_period_return(cfs, dts, basis=365.0)

    # Standard (real)
    r_std_real, pr_std_real, base_cpi_std, base_date_std = annual_and_period_return_real(
        cfs, dts, cpis, basis=365.0
    )

    # Joakim (nominal): end window per-group
    window_end_group = max(dts)
    r_joi, pr_joi, moved, moved_from, moved_to = xirr_joakim(
        cfs, dts, window_end=window_end_group, move_only_if_after=True, basis=365.0
    )

    # Joakim (real)
    r_joi_real, pr_joi_real, moved2, moved_from2, moved_to2, base_cpi_joi, base_date_joi = xirr_joakim_real(
        cfs, dts, cpis, window_end=window_end_group, move_only_if_after=True, basis=365.0
    )

    results.append({
        **key_dict,                      # <- keep original group columns
        "n_flows": len(dts),
        "start_date": min(dts),
        "end_date": max(dts),

        # nominal
        "annual_irr": r_std,
        "period_return": pr_std,

        # nominal (Joakim)
        "annual_irr_joakim": r_joi,
        "period_return_joakim": pr_joi,

        # real
        "annual_irr_real": r_std_real,
        "period_return_real": pr_std_real,

        # real (Joakim)
        "annual_irr_joakim_real": r_joi_real,
        "period_return_joakim_real": pr_joi_real,
    })

out_df = (
    pd.DataFrame(results)
      .sort_values(group_cols)
      .reset_index(drop=True)
)


xlsx_path = Path("./irr_slflokkun_uppruni.xlsx")
with pd.ExcelWriter(xlsx_path, engine="xlsxwriter") as wr:
    out_df.to_excel(wr, index=False, sheet_name="IRR_by_gerdeink")


print(str(xlsx_path))
