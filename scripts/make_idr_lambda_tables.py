#!/usr/bin/env python
"""Generate the shipped IDR parameter sets: ``topo/parameters/data/idr/<model>.yaml``.

Each file holds one model's residue-pair table ``lambda_ij`` (20 x 20, 4 decimals) and
its energy scale ``epsilon_kj`` for the Ashbaugh-Hatch IDR force

    U = 4 eps [(s/r)^12 - (s/r)^6] + (1 - lambda_ij) eps      r <= 2^(1/6) s
    U = lambda_ij * 4 eps [(s/r)^12 - (s/r)^6]                r >  2^(1/6) s

The files are DERIVED. Do not edit them by hand: change the inputs below and re-run

    python scripts/make_idr_lambda_tables.py

``tests/test_idr_models.py`` checks the shipped files against :func:`build_tables`.

bt
--
``lambda_ij = BT_IDR_SCALE * eps_BT(i, j) / EPSILON_KJ`` with ``eps_BT`` the shifted
Betancourt-Thirumalai energy topo uses for folded contacts (``4.184 * |BT - 0.6|``
kJ/mol). ``BT_IDR_SCALE = 0.10`` was calibrated on the 18 genuinely disordered
proteins of a 24-protein SAXS set (the other 6 are foldable globular proteins whose
published Rg is a folded-state value). Fitting Rg = R0 N^nu:

    bt (scale 0.10)  nu 0.566   R0 0.223   RMS 12.0%   r 0.89
    experiment       nu 0.551   R0 0.244

The theta point sits near scale 0.32; above it the chain collapses.

hps_urry
--------
``lambda_ij = (lambda_i + lambda_j) / 2`` from the per-residue Urry hydropathy of
Regy, Thompson, Kim & Mittal, Protein Sci. 30, 1371 (2021),
https://doi.org/10.1002/pro.4094 (mu = 1, shift 0.08 already applied). GLU is -0.08, so
GLU-GLU is negative: that pair is purely repulsive, as published.

calvados2
---------
``lambda_ij = (lambda_i + lambda_j) / 2`` from the per-residue CALVADOS 2 stickiness of
Tesei & Lindorff-Larsen, Open Res. Europe 2, 94 (2022),
https://doi.org/10.12688/openreseurope.14967.2, copied from the ``lambdas`` column of
``examples/single_RNA/residues_C2RNA.csv`` in https://github.com/KULL-Centre/CALVADOS
(the 20 amino acids only). Every value is positive; GLU is 0.0007, nearly repulsive.
Only lambda is taken: the CALVADOS sigma, charges, bond length and cutoff are not.

All sets use ``epsilon = 0.8368 kJ/mol`` (0.2 kcal/mol). Bead sizes are not part of a
model: they come from topo's per-residue ``Rmin/2``.
"""
from pathlib import Path

import yaml

from topo.utils.nonbonded import load_bt_potential

EPSILON_KJ = 0.8368
BT_IDR_SCALE = 0.10
DECIMALS = 4

RESIDUES = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

HPS_URRY_LAMBDA = {
    "ALA": 0.522942, "ARG": 0.478824, "ASN": 0.508236, "ASP": 0.214119,
    "CYS": 0.56706, "GLN": 0.478824, "GLU": -0.08, "GLY": 0.49353,
    "HIS": 0.684707, "ILE": 0.625883, "LEU": 0.640589, "LYS": 0.302354,
    "MET": 0.596471, "PHE": 0.74353, "PRO": 0.678824, "SER": 0.508236,
    "THR": 0.508236, "TRP": 0.92, "TYR": 0.817059, "VAL": 0.584707,
}

CALVADOS2_LAMBDA = {
    "ALA": 0.2743297969, "ARG": 0.7307624768, "ASN": 0.4255859010,
    "ASP": 0.0416040481, "CYS": 0.5615435099, "GLN": 0.3934318551,
    "GLU": 0.0006935461, "GLY": 0.7058843734, "HIS": 0.4663667291,
    "ILE": 0.5423623611, "LEU": 0.6440005008, "LYS": 0.1790211739,
    "MET": 0.5308481134, "PHE": 0.8672358982, "PRO": 0.3593126576,
    "SER": 0.4625416812, "THR": 0.3713162976, "TRP": 0.9893764740,
    "TYR": 0.9774611449, "VAL": 0.2083769608,
}

OUT_DIR = Path(__file__).resolve().parents[1] / "topo" / "parameters" / "data" / "idr"


def _pair_table(value):
    """20 x 20 nested dict of ``round(value(a, b), DECIMALS)``; ``+ 0.0`` clears -0.0."""
    return {a: {b: round(float(value(a, b)), DECIMALS) + 0.0 for b in RESIDUES}
            for a in RESIDUES}


def build_tables():
    """Return ``{model: {name, citation, epsilon_kj, source, lambda}}``."""
    bt = load_bt_potential()          # kJ/mol, 4.184 * |BT - 0.6|
    return {
        "bt": {
            "name": "bt",
            "citation": "Betancourt & Thirumalai, Protein Sci. 8, 361 (1999); "
                        "scaled for topo disordered regions",
            "epsilon_kj": EPSILON_KJ,
            "source": f"lambda_ij = {BT_IDR_SCALE} * 4.184 * |BT_ij - 0.6| / {EPSILON_KJ}",
            "lambda": _pair_table(lambda a, b: BT_IDR_SCALE * bt.loc[a, b] / EPSILON_KJ),
        },
        "hps_urry": {
            "name": "hps_urry",
            "citation": "Regy, Thompson, Kim & Mittal, Protein Sci. 30, 1371 (2021), "
                        "https://doi.org/10.1002/pro.4094",
            "epsilon_kj": EPSILON_KJ,
            "source": "lambda_ij = (lambda_i + lambda_j) / 2, Urry scale (shift 0.08 applied)",
            "lambda": _pair_table(
                lambda a, b: 0.5 * (HPS_URRY_LAMBDA[a] + HPS_URRY_LAMBDA[b])),
        },
        "calvados2": {
            "name": "calvados2",
            "citation": "Tesei & Lindorff-Larsen, Open Res. Europe 2, 94 (2022), "
                        "https://doi.org/10.12688/openreseurope.14967.2",
            "epsilon_kj": EPSILON_KJ,
            "source": "lambda_ij = (lambda_i + lambda_j) / 2, CALVADOS 2 stickiness",
            "lambda": _pair_table(
                lambda a, b: 0.5 * (CALVADOS2_LAMBDA[a] + CALVADOS2_LAMBDA[b])),
        },
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for model, table in build_tables().items():
        path = OUT_DIR / f"{model}.yaml"
        with open(path, "w") as fh:
            fh.write("# GENERATED by scripts/make_idr_lambda_tables.py -- do not edit.\n")
            yaml.safe_dump(table, fh, sort_keys=False, default_flow_style=None, width=200)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
