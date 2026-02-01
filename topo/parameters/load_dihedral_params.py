"""
Load dihedral parameters from CSV for the TOPO model.

Returns a dict keyed by "('res1', 'res2', period)" with values [period, delta, k_D]
as expected by topo.core.system.addPeriodicTorsionForce().
"""
import csv
from pathlib import Path


def load_dihedral_params():
    """Load dihedral_params from topo/parameters/data/dihedral_params.csv."""
    data_dir = Path(__file__).resolve().parent / "data"
    csv_path = data_dir / "dihedral_params.csv"
    result = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            res1 = row["res1"]
            res2 = row["res2"]
            period = int(row["period"])
            delta = float(row["delta"])
            k_D = float(row["k_D"])*0.756 #scaled by 0.756 for consistent with Yang code:        https://github.com/obrien-lab/cg_simtk_protein_folding/blob/dfa3a6f9a732cd0eb908934ba2acd4e1c9f6bbf5/CG_protein_parameterization/create_cg_protein_model.py#L706
            key = f"('{res1}', '{res2}', {period})"
            result[key] = [period, delta, k_D]
    return result
