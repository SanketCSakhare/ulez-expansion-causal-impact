"""Did the August 2023 London-wide ULEZ expansion reduce NO2 in outer London?

Difference-in-differences on monthly NO2 at London Air Quality Network monitors:
outer-London sites (newly charged from 29 Aug 2023) against inner-London sites
(inside the zone since Oct 2021). Run: python pipeline.py
"""
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import statsmodels.formula.api as smf

ROOT = Path(__file__).parent
DATA, REPORTS = ROOT / "data", ROOT / "reports"
API = "https://api.erg.ic.ac.uk/AirQuality"
START, END = "2022-09-01", "2024-09-01"  # 12 months before and 12 after, same calendar months each side
POLICY_MONTH = "2023-09"                 # expansion began 29 Aug 2023; September is the first full month
MIN_CAPTURE = 0.75                       # share of hourly readings a site-month needs to count
MIN_MONTHS = 20                          # valid months (of 24) a site needs to stay in the panel

# Boroughs wholly outside / wholly inside the North and South Circular roads (the 2021 boundary).
# Boroughs the boundary cuts through are left out because site-level treatment status is ambiguous.
OUTER = {"Barking and Dagenham", "Bexley", "Bromley", "Croydon", "Harrow", "Havering",
         "Hillingdon", "Kingston", "Merton", "Sutton"}
INNER = {"Camden", "City of London", "Hackney", "Hammersmith and Fulham", "Islington",
         "Kensington and Chelsea", "Lambeth", "Southwark", "Tower Hamlets", "Westminster"}


def assign_group(borough: str):
    return "outer" if borough in OUTER else "inner" if borough in INNER else None


def candidate_sites() -> pd.DataFrame:
    js = requests.get(f"{API}/Information/MonitoringSiteSpecies/GroupName=London/Json", timeout=120).json()
    rows = []
    for s in js["Sites"]["Site"]:
        species = s["Species"] if isinstance(s["Species"], list) else [s["Species"]]
        no2 = [sp for sp in species if sp["@SpeciesCode"] == "NO2"]
        group = assign_group(s["@LocalAuthorityName"])
        if not no2 or group is None:
            continue
        started, finished = no2[0]["@DateMeasurementStarted"], no2[0]["@DateMeasurementFinished"]
        if started <= START and (not finished or finished >= END):
            rows.append({"site": s["@SiteCode"], "borough": s["@LocalAuthorityName"],
                         "site_type": s["@SiteType"], "group": group})
    return pd.DataFrame(rows)


def monthly_no2(site: str) -> pd.DataFrame:
    """Hourly NO2 for one site, aggregated to months that meet the data-capture threshold."""
    cache = DATA / "sites" / f"{site}.csv"  # per-site cache so an interrupted run resumes
    if cache.exists():
        return pd.read_csv(cache)
    edges = [str(d.date()) for d in pd.date_range(START, END, freq="MS")]
    data = []
    for a, b in zip(edges, edges[1:]):  # one month per request: the API stalls on long date ranges
        url = f"{API}/Data/SiteSpecies/SiteCode={site}/SpeciesCode=NO2/StartDate={a}/EndDate={b}/Json"
        data += requests.get(url, timeout=120).json()["RawAQData"].get("Data", [])
    h = pd.DataFrame(data, columns=["@MeasurementDateGMT", "@Value"]).drop_duplicates("@MeasurementDateGMT")
    h["value"] = pd.to_numeric(h["@Value"], errors="coerce")
    h["month"] = h["@MeasurementDateGMT"].str[:7]
    m = h.groupby("month").value.agg(no2="mean", capture=lambda v: v.notna().mean()).reset_index()
    m = m[m.capture >= MIN_CAPTURE].assign(site=site)[["site", "month", "no2"]]
    cache.parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(cache, index=False)
    print(f"fetched {site}", flush=True)
    return m


def build_panel() -> pd.DataFrame:
    cache = DATA / "panel.csv"
    if cache.exists():
        return pd.read_csv(cache)
    sites = candidate_sites()
    with ThreadPoolExecutor(4) as pool:  # the API is slow per request; four at a time keeps the load modest
        frames = list(pool.map(monthly_no2, sites.site))
    panel = pd.concat(frames, ignore_index=True).merge(sites, on="site")
    panel = panel[panel.groupby("site").month.transform("count") >= MIN_MONTHS]
    DATA.mkdir(exist_ok=True)
    panel.to_csv(cache, index=False)
    return panel


def add_design(panel: pd.DataFrame) -> pd.DataFrame:
    p = panel.astype({"no2": float})  # empty site files make the column load as object
    p["treated"] = (p.group == "outer").astype(int)
    p["post"] = (p.month >= POLICY_MONTH).astype(int)
    p["did"] = p.treated * p.post
    months = sorted(p.month.unique())
    p["rel"] = p.month.map({m: i for i, m in enumerate(months)}) - months.index(POLICY_MONTH)
    return p


def did(p: pd.DataFrame, outcome: str = "no2") -> dict:
    """Two-way fixed effects DiD with standard errors clustered by monitoring site."""
    fit = smf.ols(f"{outcome} ~ did + C(site) + C(month)", p).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(p.site)[0]})
    lo, hi = fit.conf_int().loc["did"]
    return {"estimate": round(float(fit.params["did"]), 3), "ci95": [round(float(lo), 3), round(float(hi), 3)],
            "p_value": round(float(fit.pvalues["did"]), 4), "sites": int(p.site.nunique()),
            "treated_sites": int(p[p.treated == 1].site.nunique()), "obs": int(len(p))}


def event_study(p: pd.DataFrame):
    """Month-by-month treated-vs-control gaps relative to the month before the policy (rel = -1)."""
    p = p.copy()
    terms = []
    for k in sorted(p.rel.unique()):
        if k == -1:
            continue
        name = f"ev_m{-k}" if k < 0 else f"ev_p{k}"
        p[name] = ((p.rel == k) & (p.treated == 1)).astype(int)
        terms.append((k, name))
    fit = smf.ols("no2 ~ " + " + ".join(n for _, n in terms) + " + C(site) + C(month)", p).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(p.site)[0]})
    ci = fit.conf_int()
    coefs = pd.DataFrame([{"rel": k, "coef": fit.params[n], "lo": ci.loc[n, 0], "hi": ci.loc[n, 1]}
                          for k, n in terms])
    return coefs


def pre_trend(p: pd.DataFrame) -> dict:
    """Pre-period only: does the outer-inner gap drift per month before the policy? Should be near zero."""
    pre = p[p.post == 0].copy()
    pre["drift"] = pre.treated * pre.rel
    fit = smf.ols("no2 ~ drift + C(site) + C(month)", pre).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(pre.site)[0]})
    return {"gap_drift_ugm3_per_month": round(float(fit.params["drift"]), 3),
            "p_value": round(float(fit.pvalues["drift"]), 4)}


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = add_design(build_panel())
    p["log_no2"] = np.log(p.no2)
    means = p.groupby(["group", "post"]).no2.mean().unstack().round(2)
    coefs = event_study(p)
    roadside = p[p.site_type.isin(["Roadside", "Kerbside"])]
    results = {
        "window": f"{p.month.min()} to {p.month.max()}", "policy_month": POLICY_MONTH,
        "mean_no2_ugm3": {g: {"pre": float(means.loc[g, 0]), "post": float(means.loc[g, 1])} for g in means.index},
        "did_levels_ugm3": did(p),
        "did_log_points": did(p, "log_no2"),
        "did_roadside_only_ugm3": did(roadside) if roadside.treated.nunique() == 2 else None,
        "pre_trend": pre_trend(p),
        "sites_by_group_and_type": {f"{g} / {t}": int(n) for (g, t), n in
                                    p.drop_duplicates("site").groupby(["group", "site_type"]).size().items()},
    }
    REPORTS.mkdir(exist_ok=True)
    (REPORTS / "results.json").write_text(json.dumps(results, indent=2))
    coefs.round(3).to_csv(REPORTS / "event_study.csv", index=False)

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    trend = p.groupby(["month", "group"]).no2.mean().unstack()
    trend.plot(ax=ax[0], marker="o", ms=3)
    ax[0].axvline(list(trend.index).index(POLICY_MONTH) - 0.5, color="black", ls="--", lw=1)
    ax[0].set(ylabel="Mean NO2 (ug/m3)", xlabel="", title="Monthly mean NO2 by group")
    ax[0].tick_params(axis="x", rotation=45)
    ax[1].errorbar(coefs.rel, coefs.coef, yerr=[coefs.coef - coefs.lo, coefs.hi - coefs.coef], fmt="o", capsize=3)
    ax[1].axhline(0, color="black", lw=1)
    ax[1].axvline(-0.5, color="black", ls="--", lw=1)
    ax[1].set(xlabel="Months relative to expansion", ylabel="Outer minus inner gap vs month -1 (ug/m3)",
              title="Event study with 95% confidence intervals")
    fig.tight_layout()
    fig.savefig(REPORTS / "did_results.png", dpi=130)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
