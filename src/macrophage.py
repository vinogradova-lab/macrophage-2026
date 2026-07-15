"""Shared macrophage-project constants and small loaders.

"""

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]

# Manuscript folder tree (read-only source data)
MANUSCRIPT = Path(
    "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Macrophage project/"
    "Manuscript 1"
)
PROTEOMICS = MANUSCRIPT / "01_Data Analysis/01_Proteomics"
# The Figure 1 "Panels/" folder is the de-facto home for the cross-modal source tables
# (volcano long-format, per-condition RNA DEG, western blot / qPCR source data).
FIGURE1_PANELS = (
    MANUSCRIPT / "02_Figures/6-figure format_EV/Figure 1_RNAseq, TMT-exp/Panels"
)
RNA_DEG_DIR = FIGURE1_PANELS / "rna_seq_deg"

# Reference databases for the cysteine cross-reference (CORUM, ComplexPortal, conservation,
# AlphaMissense, pPSE) — see src/cross_ref.py.
REFERENCE_DBS = PROTEOMICS / "00_reference lists/reference_databases"
# ClinVar keys on gene symbol, not uniprot. Lives outside the project tree.
CLINVAR_TXT = Path("/Users/henrysanford/dev/cys-base-edit-design/variant_summary.txt")

# Current reactivity-changes run. ``reactivity_changes_long_format.csv`` is the canonical output
# (one row per uniprot/residue/condition); the sibling wide ``reactivity_changes.csv`` is a QC
# artifact that drops proteins — see src/reactivity.py.
RC_OUTPUT = (
    PROTEOMICS / "02_reactivity/03_rc_analysis/reactivity_output"
    / "20251023_ratio_directionality/median_for_5+_peptides/output"
)
RC_LONG_FORMAT = RC_OUTPUT / "reactivity_changes_long_format.csv"

# Stimulus conditions

CONDITIONS = ["TLR1-2", "STING", "TLR3", "TLR4", "TLR7", "TLR8", "TLR9"]
# The reactivity experiment covers six of the seven CONDITIONS (no TLR9), in the run's own order.
RC_CONDITIONS = ["TLR1-2", "TLR3", "TLR4", "TLR7", "TLR8", "STING"]

# Per-condition DESeq2 result file (sheet ``DEG_SHEET``) under ``RNA_DEG_DIR``.
DEG_FILE = {
    "TLR1-2": "DE_TLR1_rd.xlsx",
    "STING": "DE_STING_rd.xlsx",
    "TLR3": "DE_TLR3_rd.xlsx",
    "TLR4": "DE_TLR4_rd.xlsx",
    "TLR7": "DE_TLR7_rd.xlsx",
    "TLR8": "DE_TLR8_rd.xlsx",
    "TLR9": "DE_TLR9_rd.xlsx",
}
DEG_SHEET = "DEseq2, coding"


def require(path, hint=None):
    """Return ``path`` if it exists, else exit with an error (optionally appending ``hint``)."""
    path = Path(path)
    if not path.exists():
        msg = f"ERROR: missing input {path}"
        if hint:
            msg += f"\n{hint}"
        raise SystemExit(msg)
    return path


def load_deg(condition):
    """Read one condition's DESeq2 ``coding`` sheet as a DataFrame (columns as in the file)."""
    return pd.read_excel(require(RNA_DEG_DIR / DEG_FILE[condition]), sheet_name=DEG_SHEET)
