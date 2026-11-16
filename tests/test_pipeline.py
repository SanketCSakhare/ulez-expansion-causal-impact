import numpy as np
import pandas as pd
from pipeline import add_design, assign_group, did


def test_ambiguous_boroughs_are_excluded():
    assert assign_group("Bromley") == "outer"
    assert assign_group("Camden") == "inner"
    assert assign_group("Ealing") is None  # North Circular cuts through it


def test_did_recovers_a_known_effect():
    rng = np.random.default_rng(0)
    months = [str(m) for m in pd.period_range("2022-09", "2024-08", freq="M")]
    rows = [{"site": f"s{i}", "group": "outer" if i < 10 else "inner", "month": m,
             "no2": 30 + i + 0.3 * j + rng.normal(0, 0.5) - (3 if i < 10 and m >= "2023-09" else 0)}
            for i in range(20) for j, m in enumerate(months)]
    res = did(add_design(pd.DataFrame(rows)))
    assert abs(res["estimate"] - (-3)) < 0.3
    assert res["ci95"][0] < -3 < res["ci95"][1]
