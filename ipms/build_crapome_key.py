"""Regenerate the CRAPome contaminant key for the macrophage IP-MS proteins.

This rebuilds the key against the macrophage protein list.

Three steps:
    1. export_uniprot_list  -- write the quantified uniprots.
    2. build_gp_from_api    -- query CRAPome v2.0 per gene for a "*_gp.txt"
       export (skipped if a manually downloaded one already sits in crapome/).
    3. build_key            -- turn that export into macrophage_crapome_key.csv.
"""

import re
import ssl
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl

analysis_directory = Path(
    "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Macrophage project/Manuscript 1/01_Data Analysis/01_Proteomics/04_IP-MS/03_analysis"
)
data_directory = analysis_directory / "census-out-processing"

crapome_dir = Path(__file__).parent / "crapome"
crapome_dir.mkdir(exist_ok=True)

UNIPROT_LIST_PATH = crapome_dir / "macrophage_uniprots_for_crapome.txt"
KEY_PATH = crapome_dir / "macrophage_crapome_key.csv"
API_GP_PATH = crapome_dir / "macrophage_from_api_gp.txt"

# CRAPome "found in > 200 of 716 experiments" defines a contaminant.
# Reproduces mouse_fasta_to_FP_TP_tables.py (the original key builder).
CONTAMINANT_FOUND_CUTOFF = 200

# CRAPome v2.0 web service. Human controls = 716; the response lists one
# <protein ... expt="CCnnn"> element per (isoform, control) detection, so the
# number of distinct expt names is the "found" count.
CRAPOME_TOTAL = 716
API_URL = "https://reprint-apms.org/?q=ws/proteindetail/{gene}/human/singleStep/2.0"
_EXPT_RE = re.compile(r'expt="([^"]+)"')
# reprint-apms.org serves an unverifiable cert chain; skip verification.
_SSL_CTX = ssl._create_unverified_context()


def export_uniprot_list() -> None:
    """Step 1: write every unique uniprot in the IP-MS data, one per line."""
    sources = ["long_volcano_data.csv", "unfiltered_data.csv"]
    uniprots: set[str] = set()
    for name in sources:
        path = data_directory / name
        if not path.exists():
            print(f"  (skipping missing source {name})")
            continue
        uniprots |= set(pl.read_csv(path)["uniprot"].drop_nulls().to_list())

    uniprots = sorted(u for u in uniprots if u and u != "None")
    UNIPROT_LIST_PATH.write_text("\n".join(uniprots) + "\n")
    print(f"Wrote {len(uniprots)} uniprots to {UNIPROT_LIST_PATH}")


def _uniprot_to_gene() -> dict[str, str]:
    """Map each quantified uniprot to its gene symbol from the IP-MS tables.

    Both source CSVs carry uniprot + protein (an uppercase human gene symbol),
    so no external lookup is needed; unfiltered_data.csv is the fuller list.
    """
    mapping: dict[str, str] = {}
    for name in ["unfiltered_data.csv", "long_volcano_data.csv"]:
        path = data_directory / name
        if not path.exists():
            continue
        pairs = (
            pl.read_csv(path)
            .select(["uniprot", "protein"])
            .drop_nulls()
            .unique()
        )
        for uni, gene in zip(pairs["uniprot"], pairs["protein"]):
            mapping.setdefault(uni, gene)
    return mapping


def _crapome_found(gene: str) -> int:
    """Distinct CRAPome control experiments the gene is detected in (0-716).

    A gene absent from CRAPome returns the site HTML (no expt= tokens) -> 0.
    """
    url = API_URL.format(gene=urllib.parse.quote(gene))
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30, context=_SSL_CTX) as resp:
                body = resp.read().decode("utf-8", "replace")
            return len(set(_EXPT_RE.findall(body)))
        except Exception:  # transient network / server hiccup -> back off
            if attempt == 2:
                raise
            time.sleep(1.5 * (attempt + 1))


def build_gp_from_api() -> None:
    """Step 2a: query CRAPome v2.0 per gene and emit a website-schema gp.txt."""
    uni2gene = _uniprot_to_gene()
    if not UNIPROT_LIST_PATH.exists():
        export_uniprot_list()
    uniprots = [
        u for u in UNIPROT_LIST_PATH.read_text().splitlines() if u.strip()
    ]

    # query once per unique gene (several uniprots can share a symbol),
    # concurrently -- the per-gene requests are independent and I/O-bound.
    genes = sorted({uni2gene[u] for u in uniprots if u in uni2gene})
    print(f"Querying CRAPome v2.0 for {len(genes)} genes "
          f"({len(uniprots)} uniprots)...")
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=16) as pool:
        counts = list(pool.map(_crapome_found, genes))
    found_by_gene = dict(zip(genes, counts))
    print(f"  done in {time.time() - t0:.0f}s")

    rows = []
    unmapped = 0
    for uni in uniprots:
        gene = uni2gene.get(uni)  # None -> unmapped, found_by_gene.get(None) == 0
        if gene is None:
            unmapped += 1
        rows.append(
            {
                "User Input": uni,
                "Mapped Gene Symbol": gene if gene is not None else "Invalid identifier",
                "Num of Expt. (found/total)": f"{found_by_gene.get(gene, 0)} / {CRAPOME_TOTAL}",
                "Ave SC": "",
                "Max SC": "",
            }
        )

    pl.DataFrame(rows).write_csv(API_GP_PATH, separator="\t")
    n_contam = sum(f > CONTAMINANT_FOUND_CUTOFF for f in found_by_gene.values())
    print(
        f"Wrote {len(rows)} rows to {API_GP_PATH} "
        f"({unmapped} uniprots without a gene symbol; "
        f"{n_contam} genes over the contaminant cutoff)."
    )


def build_key() -> None:
    """Step 3: turn a Crapome "*_gp.txt" export into the macrophage key."""
    exports = sorted(crapome_dir.glob("*_gp.txt"))
    if not exports:
        print(
            f"No '*_gp.txt' Crapome export found in {crapome_dir}.\n"
            "Submit macrophage_uniprots_for_crapome.txt to Crapome.org, download "
            "the gp.txt export into that folder, then re-run this script."
        )
        return

    export_path = exports[-1]  # most recent by name (timestamped)
    print(f"Building key from {export_path.name}")

    key = (
        pl.read_csv(export_path, separator="\t")
        # drop identifiers Crapome could not map
        .filter(
            pl.col("Mapped Gene Symbol")
            .str.contains("Invalid identifier")
            .not_()
        )
        .with_columns(
            pl.col("Num of Expt. (found/total)")
            .str.split(" / ")
            .list.get(0)
            .str.strip_chars()
            .replace("", "0")
            .cast(pl.Float64, strict=False)
            .fill_null(0)
            .alias("found")
        )
        .with_columns(
            (pl.col("found") > CONTAMINANT_FOUND_CUTOFF).alias("contaminant")
        )
        .select(
            "User Input",
            "Mapped Gene Symbol",
            "Num of Expt. (found/total)",
            "found",
            "contaminant",
        )
    )

    key.write_csv(KEY_PATH)
    n_contam = int(key["contaminant"].sum())
    print(f"Wrote {key.height} rows ({n_contam} contaminants) to {KEY_PATH}")


if __name__ == "__main__":
    export_uniprot_list()
    # If no gp.txt export is present yet, build one automatically from the
    # CRAPome v2.0 API; a manually downloaded "*_gp.txt" takes precedence.
    if not sorted(crapome_dir.glob("*_gp.txt")):
        build_gp_from_api()
    build_key()
