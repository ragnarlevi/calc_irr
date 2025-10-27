

from datetime import datetime, date
from math import isfinite

# ---------- helpers ----------
def _to_date(d):
    if isinstance(d, (date, datetime)): 
        return d.date() if isinstance(d, datetime) else d
    return datetime.strptime(d, "%d/%m/%Y").date()

def _xnpv(r, cf, T):  # T in years
    return sum(c / (1.0 + r) ** t for c, t in zip(cf, T))

def _dxnpv(r, cf, T):
    return sum(-c * t * (1.0 + r) ** (-(t + 1.0)) for c, t in zip(cf, T))

def _solve_xirr(cf, T, guess=0.1, tol=1e-12):
    # Newton
    r = guess
    for _ in range(50):
        f, df = _xnpv(r, cf, T), _dxnpv(r, cf, T)
        if abs(df) < 1e-18: break
        r2 = r - f/df
        if not (-0.9999 < r2 < 1e12): break
        if abs(r2 - r) < tol: return r2
        r = r2
    # Bisection fallback
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

# ---------- main ----------
def period_irr(
    cashflows,            # iterable of (date|"%d/%m/%Y", amount)
    start, end,           # date|str bounds, inclusive
    basis=365.0,          # day-count divisor
    opening=None,         # optional start valuation (adds -opening at start)
    closing=None,         # optional end valuation (adds +closing at end)
    guess=0.1
):
    """
    Returns (annual_IRR, period_return) for cashflows in [start, end].
    - annual_IRR is XIRR (per year)
    - period_return is non-annualized over the window length
    """
    start = _to_date(start); end = _to_date(end)
    assert start <= end, "start must be <= end"

    # filter to window
    dts, cf = [], []
    for d, v in cashflows:
        d0 = _to_date(d)
        if start <= d0 <= end:
            dts.append(d0); cf.append(float(v))

    # synthetic valuation flows (optional)
    if opening is not None:
        dts.append(start); cf.append(-float(opening))
    if closing is not None:
        dts.append(end);   cf.append(+float(closing))

    # need at least one pos and one neg to have an IRR root
    if not dts or min(cf) >= 0 or max(cf) <= 0:
        return float("nan"), float("nan")

    T = [(d - start).days / basis for d in dts]
    r = _solve_xirr(cf, T, guess=guess)

    if not isfinite(r): 
        return float("nan"), float("nan")

    D = (end - start).days
    pr = (1.0 + r) ** (D / basis) - 1.0  # de-annualize to YTD/period return
    return r, pr

# ---------- example ----------
if __name__ == "__main__":


    rows = [
        ("01/01/2025", -1_375_102_501.00),
        ("05/05/2025",    2_697_324.00),
        ("15/05/2025", -103_737_022.00),
        ("16/06/2025",   14_924_566.00),
        ("22/09/2025",    2_206_929.00),
        ("01/10/2025", 1_538_239_368.00),
    ]


    rows = [
        ("01/01/2025", -186346780.00),
        ("05/05/2025",    2697324.00),
        ("01/10/2025", 195382531.00),
    ]
    r, pr = period_irr(rows, start="01/01/2025", end="01/10/2025")
    print(f"Annual IRR: {r*100:.4f}%")
    print(f"Period return: {pr*100:.4f}%")


    import datetime 
    from math import pow
    from scipy.optimize import newton

    # def xirr(cashflows, dates, guess=0.1):
    #     # Convert to days from first date
    #     days = [(d - dates[0]).days for d in dates]
    #     print(days)

    #     def npv(r):
    #         return sum(cf / pow(1 + r, d / 365.0) for cf, d in zip(cashflows, days))

    #     return newton(npv, guess)


    import datetime as _dt
    import math
    from collections import defaultdict
    from scipy.optimize import brentq

    def _prep(cashflows, dates):
        # pair, aggregate same day, drop zeros, sort
        agg = defaultdict(float)
        for cf, d in zip(cashflows, dates):
            if isinstance(d, _dt.datetime):  # normalize to date
                d = d.date()
            agg[d] += float(cf)
        items = sorted((d, cf) for d, cf in agg.items() if cf != 0.0)
        if not items:
            raise ValueError("All cash flows are zero.")
        d0 = items[0][0]
        days = [(d - d0).days for d, _ in items]
        cfs  = [cf for _, cf in items]
        if not (any(cf > 0 for cf in cfs) and any(cf < 0 for cf in cfs)):
            raise ValueError("Need at least one positive and one negative cash flow.")
        return cfs, days

    def _xnpv(rate, cfs, days):
        # Only defined for rate > -1. Use exp/log1p for stability.
        if rate <= -1.0:
            return math.copysign(float("inf"), -1.0)  # push out of domain
        lr = math.log1p(rate)
        return sum(cf / math.exp(lr * (d / 365.0)) for cf, d in zip(cfs, days))

    def xirr(cashflows, dates, guess=0.1):
        cfs, days = _prep(cashflows, dates)

        f = lambda r: _xnpv(r, cfs, days)

        # Find a bracket [lo, hi] with opposite signs
        lo = -0.999999999  # just inside domain
        hi = max(guess, 0.1)

        flo = f(lo + 1e-12)  # nudge inside domain
        fhi = f(hi)
        # Expand until sign change or limits reached
        for _ in range(80):
            if math.isfinite(flo) and math.isfinite(fhi) and flo * fhi < 0:
                break
            # expand outward: move the side with smaller |f| less aggressively
            hi = hi * 2.0 + 0.05
            fhi = f(hi)
            # also try moving lo closer to -1 (can’t cross -1)
            lo = -1.0 + (lo + 1.0) * 0.5
            flo = f(lo + 1e-12)
        else:
            # As a last resort, try a wide fixed bracket in plausible range
            lo, hi = -0.999999999, 10.0
            flo, fhi = f(lo + 1e-12), f(hi)
            if not (math.isfinite(flo) and math.isfinite(fhi) and flo * fhi < 0):
                raise RuntimeError("Failed to bracket the IRR root. Check cash flows/dates.")

        # Solve with Brent's method on a valid domain
        root = brentq(lambda r: f(r), lo + 1e-12, hi, maxiter=200, xtol=1e-12, rtol=1e-12)
        return root

    # ---- Example ----
    cfs = [-153961188128.7, 0.0, -522900.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 29514282.0, 0.0, 0.0, 0.0, -19664365.0, 0.0, 0.0, 0.0, 0.0, -11590721.0, -3394688.0, 0.0, 0.0, 0.0, 0.0, -75930777.0, 0.0, 0.0, -269567961.0, 0.0, 0.0, 92578418.0, -15898840.0, 0.0, -23457971.0, -38142900.0, 0.0, 0.0, -21986705.0, 40017478.0, 0.0, 0.0, -61545783.0, 0.0, 87058.0, 0.0, -48098484.0, -27664289.0, 0.0, 0.0, 0.0, 4000099.0, 0.0, 0.0, 56667492.0, -13405671.0, 0.0, 0.0, 3017617.0, 583958.0, -33290242.0, 13360251.59, 0.0, 0.0, 0.0, 0.0, 52807608.52, -43649721.61, 0.0, 162914152.0, 0.0, 0.0, 0.0, 115458914.0, -3217880.0, 29731033.0, 0.0, 0.0, 0.0, 190341652.0, -57727435.0, 384632949.0, -51196501.0, 0.0, 0.0, 0.0, -30301987.0, -499500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 60151841.0, 0.0, 0.0, 0.0, 0.0, 0.0, -71803028.0, 0.0, -48284011.0, 0.0, 0.0, -94480854.0, -6640162.0, -43579276.0, -28582114.0, -35626783.0, 0.0, 0.0, -8544587.0, 0.0, 0.0, -38503601.0, -24951340.0, 0.0, 0.0, 0.0, 0.0, 0.0, -11981342.0, -3447670.0, 0.0, 41626152.0, 0.0, 0.0, -9566250.0, 38806283.0, 0.0, 0.0, 12966732.0, 4500046.0, 176082312.0, 0.0, -8543219.0, 0.0, 0.0, 3133998.0, 619965000.0, 0.0, 35769584.0, 0.0, 0.0, 0.0, -5183950.0, 43904000.0, 0.0, 455701785.0, 0.0, 0.0, 0.0, 0.0, 13478530.0, 82576692.0, -33606998.0, 0.0, 0.0, 0.0, 0.0, -53573146.0, -181107574.0, 0.0, 0.0, 0.0, 501214039.0, -24870905.0, 0.0, -1605988.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -14972144.0, 8789981.0, -30547600.0, -4128525.0, 0.0, 41650363.0, 0.0, -34484095.0, 0.0, 0.0, 0.0, 0.0, -7134919.0, 614150000.0, 30462185.0, 0.0, 0.0, 0.0, 0.0, 0.0, 259189525.52, -172290.0, 2290891.0, 0.0, 0.0, 0.0, -1326748.84, 0.0, -22637195.0, -508620.0, -100376988.0, 0.0, 0.0, 0.0, 205589310.25, 47472625.0, -33807877.0, 0.0, 0.0, 0.0, -18607369.0, -81699635.0, 0.0, -13568500.0, 0.0, 0.0, 0.0, 0.0, 11436482.0, -24510515.0, 0.0, 0.0, 0.0, 0.0, -93374303.0, 0.0, -1568040.0, -14881046.0, 107785495.0, 0.0, 0.0, 0.0, 0.0, 0.0, -39895769.0, -70700051.0, 0.0, -2875559.0, 0.0, -27342430.0, -77464321.0, -109128517.0, 0.0, 0.0, -38864540.0, -12873334.0, -456113.0, 0.0, -21352970.0, 0.0, 0.0, 87134400.0, 65578199.0, 0.0, 10566423.0, 264009999.0, 0.0, 0.0, 148263635.0, 35220000.0, 0.0, 0.0, 0.0, 306972048.0, 155537587901.6]

    dates = [datetime.date(2025, 1, 1), datetime.date(2025, 1, 1), datetime.date(2025, 1, 2), datetime.date(2025, 1, 3), datetime.date(2025, 1, 4), datetime.date(2025, 1, 5), datetime.date(2025, 1, 6), datetime.date(2025, 1, 7), datetime.date(2025, 1, 8), datetime.date(2025, 1, 9), datetime.date(2025, 1, 10), datetime.date(2025, 1, 11), datetime.date(2025, 1, 12), datetime.date(2025, 1, 13), datetime.date(2025, 1, 14), datetime.date(2025, 1, 15), datetime.date(2025, 1, 16), datetime.date(2025, 1, 17), datetime.date(2025, 1, 18), datetime.date(2025, 1, 19), datetime.date(2025, 1, 20), datetime.date(2025, 1, 21), datetime.date(2025, 1, 22), datetime.date(2025, 1, 23), datetime.date(2025, 1, 24), datetime.date(2025, 1, 25), datetime.date(2025, 1, 26), datetime.date(2025, 1, 27), datetime.date(2025, 1, 28), datetime.date(2025, 1, 29), datetime.date(2025, 1, 30), datetime.date(2025, 1, 31), datetime.date(2025, 2, 1), datetime.date(2025, 2, 2), datetime.date(2025, 2, 3), datetime.date(2025, 2, 4), datetime.date(2025, 2, 5), datetime.date(2025, 2, 6), datetime.date(2025, 2, 7), datetime.date(2025, 2, 8), datetime.date(2025, 2, 9), datetime.date(2025, 2, 10), datetime.date(2025, 2, 11), datetime.date(2025, 2, 12), datetime.date(2025, 2, 13), datetime.date(2025, 2, 14), datetime.date(2025, 2, 15), datetime.date(2025, 2, 17), datetime.date(2025, 2, 18), datetime.date(2025, 2, 19), datetime.date(2025, 2, 20), datetime.date(2025, 2, 21), datetime.date(2025, 2, 22), datetime.date(2025, 2, 23), datetime.date(2025, 2, 24), datetime.date(2025, 2, 25), datetime.date(2025, 2, 26), datetime.date(2025, 2, 27), datetime.date(2025, 2, 28), datetime.date(2025, 3, 1), datetime.date(2025, 3, 2), datetime.date(2025, 3, 3), datetime.date(2025, 3, 4), datetime.date(2025, 3, 5), datetime.date(2025, 3, 6), datetime.date(2025, 3, 7), datetime.date(2025, 3, 8), datetime.date(2025, 3, 9), datetime.date(2025, 3, 10), datetime.date(2025, 3, 11), datetime.date(2025, 3, 12), datetime.date(2025, 3, 13), datetime.date(2025, 3, 14), datetime.date(2025, 3, 15), datetime.date(2025, 3, 16), datetime.date(2025, 3, 17), datetime.date(2025, 3, 18), datetime.date(2025, 3, 19), datetime.date(2025, 3, 20), datetime.date(2025, 3, 21), datetime.date(2025, 3, 22), datetime.date(2025, 3, 23), datetime.date(2025, 3, 24), datetime.date(2025, 3, 25), datetime.date(2025, 3, 26), datetime.date(2025, 3, 27), datetime.date(2025, 3, 28), datetime.date(2025, 3, 29), datetime.date(2025, 3, 30), datetime.date(2025, 3, 31), datetime.date(2025, 4, 1), datetime.date(2025, 4, 2), datetime.date(2025, 4, 3), datetime.date(2025, 4, 4), datetime.date(2025, 4, 5), datetime.date(2025, 4, 6), datetime.date(2025, 4, 7), datetime.date(2025, 4, 8), datetime.date(2025, 4, 9), datetime.date(2025, 4, 10), datetime.date(2025, 4, 11), datetime.date(2025, 4, 12), datetime.date(2025, 4, 13), datetime.date(2025, 4, 14), datetime.date(2025, 4, 15), datetime.date(2025, 4, 16), datetime.date(2025, 4, 17), datetime.date(2025, 4, 21), datetime.date(2025, 4, 22), datetime.date(2025, 4, 23), datetime.date(2025, 4, 24), datetime.date(2025, 4, 25), datetime.date(2025, 4, 26), datetime.date(2025, 4, 27), datetime.date(2025, 4, 28), datetime.date(2025, 4, 29), datetime.date(2025, 4, 30), datetime.date(2025, 5, 1), datetime.date(2025, 5, 2), datetime.date(2025, 5, 3), datetime.date(2025, 5, 4), datetime.date(2025, 5, 5), datetime.date(2025, 5, 6), datetime.date(2025, 5, 7), datetime.date(2025, 5, 8), datetime.date(2025, 5, 9), datetime.date(2025, 5, 10), datetime.date(2025, 5, 11), datetime.date(2025, 5, 12), datetime.date(2025, 5, 13), datetime.date(2025, 5, 14), datetime.date(2025, 5, 15), datetime.date(2025, 5, 16), datetime.date(2025, 5, 17), datetime.date(2025, 5, 19), datetime.date(2025, 5, 20), datetime.date(2025, 5, 21), datetime.date(2025, 5, 22), datetime.date(2025, 5, 23), datetime.date(2025, 5, 24), datetime.date(2025, 5, 25), datetime.date(2025, 5, 26), datetime.date(2025, 5, 27), datetime.date(2025, 5, 28), datetime.date(2025, 5, 29), datetime.date(2025, 5, 30), datetime.date(2025, 5, 31), datetime.date(2025, 6, 1), datetime.date(2025, 6, 2), datetime.date(2025, 6, 3), datetime.date(2025, 6, 4), datetime.date(2025, 6, 5), datetime.date(2025, 6, 6), datetime.date(2025, 6, 7), datetime.date(2025, 6, 8), datetime.date(2025, 6, 9), datetime.date(2025, 6, 10), datetime.date(2025, 6, 11), datetime.date(2025, 6, 12), datetime.date(2025, 6, 13), datetime.date(2025, 6, 14), datetime.date(2025, 6, 15), datetime.date(2025, 6, 16), datetime.date(2025, 6, 18), datetime.date(2025, 6, 19), datetime.date(2025, 6, 20), datetime.date(2025, 6, 21), datetime.date(2025, 6, 22), datetime.date(2025, 6, 23), datetime.date(2025, 6, 24), datetime.date(2025, 6, 25), datetime.date(2025, 6, 26), datetime.date(2025, 6, 27), datetime.date(2025, 6, 28), datetime.date(2025, 6, 29), datetime.date(2025, 6, 30), datetime.date(2025, 7, 1), datetime.date(2025, 7, 2), datetime.date(2025, 7, 3), datetime.date(2025, 7, 4), datetime.date(2025, 7, 5), datetime.date(2025, 7, 6), datetime.date(2025, 7, 7), datetime.date(2025, 7, 8), datetime.date(2025, 7, 9), datetime.date(2025, 7, 10), datetime.date(2025, 7, 11), datetime.date(2025, 7, 12), datetime.date(2025, 7, 14), datetime.date(2025, 7, 15), datetime.date(2025, 7, 16), datetime.date(2025, 7, 17), datetime.date(2025, 7, 18), datetime.date(2025, 7, 21), datetime.date(2025, 7, 22), datetime.date(2025, 7, 23), datetime.date(2025, 7, 24), datetime.date(2025, 7, 25), datetime.date(2025, 7, 26), datetime.date(2025, 7, 27), datetime.date(2025, 7, 28), datetime.date(2025, 7, 29), datetime.date(2025, 7, 30), datetime.date(2025, 7, 31), datetime.date(2025, 8, 1), datetime.date(2025, 8, 2), datetime.date(2025, 8, 3), datetime.date(2025, 8, 4), datetime.date(2025, 8, 5), datetime.date(2025, 8, 6), datetime.date(2025, 8, 7), datetime.date(2025, 8, 8), datetime.date(2025, 8, 9), datetime.date(2025, 8, 10), datetime.date(2025, 8, 11), datetime.date(2025, 8, 12), datetime.date(2025, 8, 13), datetime.date(2025, 8, 14), datetime.date(2025, 8, 15), datetime.date(2025, 8, 16), datetime.date(2025, 8, 17), datetime.date(2025, 8, 18), datetime.date(2025, 8, 19), datetime.date(2025, 8, 20), datetime.date(2025, 8, 21), datetime.date(2025, 8, 22), datetime.date(2025, 8, 23), datetime.date(2025, 8, 24), datetime.date(2025, 8, 25), datetime.date(2025, 8, 26), datetime.date(2025, 8, 27), datetime.date(2025, 8, 28), datetime.date(2025, 8, 29), datetime.date(2025, 8, 30), datetime.date(2025, 8, 31), datetime.date(2025, 9, 1), datetime.date(2025, 9, 2), datetime.date(2025, 9, 3), datetime.date(2025, 9, 4), datetime.date(2025, 9, 5), datetime.date(2025, 9, 6), datetime.date(2025, 9, 7), datetime.date(2025, 9, 8), datetime.date(2025, 9, 9), datetime.date(2025, 9, 10), datetime.date(2025, 9, 11), datetime.date(2025, 9, 12), datetime.date(2025, 9, 13), datetime.date(2025, 9, 14), datetime.date(2025, 9, 15), datetime.date(2025, 9, 16), datetime.date(2025, 9, 17), datetime.date(2025, 9, 18), datetime.date(2025, 9, 19), datetime.date(2025, 9, 21), datetime.date(2025, 9, 22), datetime.date(2025, 9, 23), datetime.date(2025, 9, 24), datetime.date(2025, 9, 25), datetime.date(2025, 9, 26), datetime.date(2025, 9, 27), datetime.date(2025, 9, 28), datetime.date(2025, 9, 29), datetime.date(2025, 9, 30), datetime.date(2025, 10, 1), datetime.date(2025, 10, 2), datetime.date(2025, 10, 3), datetime.date(2025, 10, 4), datetime.date(2025, 10, 5), datetime.date(2025, 10, 6), datetime.date(2025, 10, 7), datetime.date(2025, 10, 8), datetime.date(2025, 10, 9), datetime.date(2025, 10, 10), datetime.date(2025, 10, 11), datetime.date(2025, 10, 12), datetime.date(2025, 10, 13), datetime.date(2025, 10, 14), datetime.date(2025, 10, 15), datetime.date(2025, 10, 16), datetime.date(2025, 10, 17), datetime.date(2025, 10, 20), datetime.date(2025, 10, 21)]
    
    print(len(cfs))
    print(len(dates))
    
    r_excel = xirr(cfs, dates)
    r_machine = xirr(cfs[:-1], dates[:-1])  # drop zero
    print(f"Excel-style (with zero): {r_excel:.4%}")
    print(f"Machine-style (truncated): {r_machine:.4%}")

