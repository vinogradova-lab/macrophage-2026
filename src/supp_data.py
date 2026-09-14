"""Shared writer for the manuscript's supplementary-data workbooks ("Data S<n>.xlsx").

Every assay's `build_data_s<n>.py` assembles its tables as pandas frames and hands them here,
so the house style lives in one place. That style is the sister repo's, matched sheet for
sheet against its delivered workbooks (`../t-cell-dysfunction-2026/supp_data/Data S3_*.xlsx`,
written by `bin/build_data_s3_2.py::build_sheet` and `bin/renumber_data_s2.py::build_sheet`):

    contents        the first tab: `#` / `Title` / `Description`, entry n on row n + 1
    row 1           `S<n>-<m> <title> (related to Fig. ...)` - the tab's number, then the
                    title, bold
    row 2           the methods sentence
    row 3           the column header
    row 4+          the data, under a striped Excel table

The three text slots are the `SuppSheet` fields, and the contents `Description` column repeats
row 1 without its `S<n>-<m> ` prefix rather than holding the row-2 sentence - which is what the
sister workbooks do, so the index reads as a list of titles. Tab names keep their `S<n>-<m> `
form because Excel caps a tab name at 31 characters, which a spelled-out
"Supplementary Data <n>-<m> ..." would exceed.

An assay script imports the two names it needs:

    from src.supp_data import SuppSheet, write_supplementary_workbook

    sheets = [SuppSheet(1, SHEET_NAME, TITLE, DESCRIPTION, build_table())]
    write_supplementary_workbook(REPO_ROOT / "Data S1.xlsx", sheets)

and is run in the `macrophage-2026` env, as `rna/build_data_s1.py` and its siblings are:

    conda run -n macrophage-2026 python rna/build_data_s1.py

Run this module directly to write a two-sheet demo workbook and verify it - a formatting smoke
test that needs no assay data, and the quickest way to eyeball the style in Excel:

    conda run -n macrophage-2026 python src/supp_data.py --out "Data S0 demo.xlsx"
"""

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo

CONTENTS_SHEET = "contents"
CONTENTS_HEADER = ["#", "Title", "Description"]

TITLE_ROW = 1
DESCRIPTION_ROW = 2
HEADER_ROW = 3

# Excel's own limit on a tab name; asserted rather than silently truncated, because openpyxl
# would otherwise write a name the contents index no longer matches.
MAX_TAB_NAME = 31

# The sister workbooks are Arial throughout, with a bold 14 pt row-1 title over a 12 pt row-2
# methods sentence. Nothing else is styled by hand: the striped table carries the header and
# body formatting, which is why no font is set below row 2.
FONT = "Arial"
TITLE_FONT = Font(name=FONT, size=14, bold=True)
DESCRIPTION_FONT = Font(name=FONT, size=12)
CONTENTS_FONT = Font(name=FONT, size=14)
CONTENTS_HEADER_FONT = Font(name=FONT, size=14, bold=True)

# Only the `Title` column is widened, as in the sister's contents tab; `#` and `Description`
# sit at the default width and overflow, which is how that index reads there.
CONTENTS_TITLE_WIDTH = 43.83203125

# `Table Style Light 8` is the style the sister repo names wherever it names one (the polars
# `write_excel` calls in `low_input.py` / `solubility.py`); its openpyxl builders ask only for
# row stripes and let Excel supply the rest. Naming it here gets the same striping without
# depending on what Excel fills in on first open.
TABLE_STYLE = "Table Style Light 8"


@dataclass
class SuppSheet:
    """One sheet in a supplementary workbook.

    ``number`` is its 1-based position in the contents index; ``sheet_name`` is the tab, whose
    first token (`S2-4`) prefixes the row-1 title; ``title`` is the row-1 body, reused verbatim
    as the contents ``Description``; ``description`` is the row-2 methods sentence; ``df`` is
    the flat table, written from row 3 (header) down.
    """

    number: int
    sheet_name: str
    title: str
    description: str
    df: pd.DataFrame

    @property
    def prefix(self):
        """The `S<n>-<m>` token the tab name leads with."""
        return self.sheet_name.split(" ", 1)[0]

    @property
    def row1(self):
        """Row 1: the tab's number, then the title."""
        return f"{self.prefix} {self.title}"


def _cell_rows(df):
    """Rows of plain Python scalars - openpyxl cannot write numpy types or NaN."""
    frame = df.astype(object).where(df.notna(), None)
    return frame.itertuples(index=False, name=None)


def _same_cell(got, want):
    """Whether a cell read back from Excel is the value that was written.

    Not `==`: Excel stores a float to 15 significant digits, so 1.2345678901234567 reads back
    as ...66, a pandas `Timestamp` comes back as a plain `datetime`, and an empty string comes
    back as None because Excel has no empty string -- openpyxl writes it as a blank cell. All
    three are faithful round trips, and none compares equal.
    """
    if want == "" and got is None:
        return True
    if isinstance(want, pd.Timestamp):
        want = want.to_pydatetime()
    if isinstance(want, bool) or isinstance(got, bool):
        return got is want
    if isinstance(want, (int, float)) and isinstance(got, (int, float)):
        return math.isclose(got, want, rel_tol=1e-12)
    return got == want


def _write_data_sheet(workbook, sheet, freeze_cols):
    """Add one data sheet: title, description, header on row 3, striped table below."""
    assert len(sheet.sheet_name) <= MAX_TAB_NAME, (
        f"tab name {sheet.sheet_name!r} is {len(sheet.sheet_name)} characters "
        f"(Excel allows {MAX_TAB_NAME})"
    )
    # An Excel table names each of its columns, so a repeated header writes a file Excel
    # reports as corrupt and offers to repair - and it writes it silently. Worth guarding for:
    # these tables are per-comparison blocks of the same metrics, so a block left untagged
    # duplicates a header rather than producing anything that looks wrong.
    columns = list(sheet.df.columns)
    duplicated = sorted({c for c in columns if columns.count(c) > 1})
    assert not duplicated, f"{sheet.sheet_name}: duplicate column headers {duplicated}"

    ws = workbook.create_sheet(title=sheet.sheet_name)

    ws.cell(row=TITLE_ROW, column=1, value=sheet.row1).font = TITLE_FONT
    ws.cell(row=DESCRIPTION_ROW, column=1, value=sheet.description).font = DESCRIPTION_FONT
    # Writing row 2 leaves openpyxl's append cursor there, so the header lands on row 3 - the
    # layout every sheet shares, with the two title rows above the header.
    ws.append(list(sheet.df.columns))
    for row in _cell_rows(sheet.df):
        ws.append(list(row))

    last_col = get_column_letter(sheet.df.shape[1])
    ws.add_table(
        Table(
            # Excel requires a workbook-unique name with no spaces; the sister workbooks call
            # theirs `Frame<nn>`.
            displayName=f"Frame{sheet.number:02d}",
            ref=f"A{HEADER_ROW}:{last_col}{HEADER_ROW + len(sheet.df)}",
            tableStyleInfo=TableStyleInfo(name=TABLE_STYLE, showRowStripes=True),
        )
    )

    # Off by default: no sister sheet freezes anything. A caller with wide value blocks passes
    # its lead-identifier count to keep those columns and the header in view.
    if freeze_cols:
        ws.freeze_panes = f"{get_column_letter(freeze_cols + 1)}{HEADER_ROW + 1}"
    return ws


def _write_contents_sheet(workbook, sheets):
    """Add the index tab: one row per sheet, in `number` order, header on row 1."""
    ws = workbook.create_sheet(title=CONTENTS_SHEET, index=0)
    ws.append(CONTENTS_HEADER)
    for sheet in sorted(sheets, key=lambda s: s.number):
        # The description column repeats the sheet's row-1 title without its `S<n>-<m> `
        # prefix - not the row-2 methods sentence.
        ws.append([sheet.number, sheet.sheet_name, sheet.title])

    for row in ws.iter_rows(min_row=1, max_row=ws.max_row):
        for cell in row:
            cell.font = CONTENTS_HEADER_FONT if cell.row == 1 else CONTENTS_FONT
    ws.column_dimensions["B"].width = CONTENTS_TITLE_WIDTH
    return ws


def write_supplementary_workbook(path, sheets, freeze_cols=0, contents=True):
    """Write a supplementary workbook: a ``contents`` index + one sheet per ``SuppSheet``.

    Sheets are written in ``number`` order, and the result is reopened and checked before this
    returns - `verify_workbook` is not optional, since a workbook that reads back wrong is the
    one failure a build script cannot see.
    """
    sheets = sorted(sheets, key=lambda s: s.number)
    numbers = [sheet.number for sheet in sheets]
    assert numbers == list(range(1, len(sheets) + 1)), f"contents numbering is {numbers}"

    workbook = Workbook()
    # A new Workbook starts with one empty sheet; every sheet here is created explicitly.
    workbook.remove(workbook.active)
    for sheet in sheets:
        _write_data_sheet(workbook, sheet, freeze_cols)
    if contents:
        _write_contents_sheet(workbook, sheets)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(path)
    workbook.close()

    verify_workbook(path, sheets, contents=contents)
    return path


def write_supplementary_sheet(df, path, sheet_name, title, description, freeze_cols=0):
    """Write a single-sheet workbook in the same style, with no contents index.

    For a table delivered on its own - the shape `bin/build_data_s3_2.py` writes in the sister
    repo - before it is folded into the assay's `Data S<n>.xlsx`.
    """
    sheet = SuppSheet(1, sheet_name, title, description, df)
    return write_supplementary_workbook(
        path, [sheet], freeze_cols=freeze_cols, contents=False
    )


def verify_workbook(path, sheets, contents=True):
    """Reopen a written workbook and check its layout against the sheets it was built from.

    Checks the structure every sheet shares - tab order, the two title rows, the header on row
    3, the table over the data - plus the first and last data row of each sheet, which is what
    catches a frame written with its columns or rows misaligned. It does not re-read every
    cell: the sheets run to tens of thousands of rows.
    """
    workbook = load_workbook(path)
    try:
        expected_tabs = [sheet.sheet_name for sheet in sheets]
        if contents:
            expected_tabs = [CONTENTS_SHEET] + expected_tabs
        assert workbook.sheetnames == expected_tabs, workbook.sheetnames

        for sheet in sheets:
            ws = workbook[sheet.sheet_name]
            assert ws[f"A{TITLE_ROW}"].value == sheet.row1, ws[f"A{TITLE_ROW}"].value
            assert ws[f"A{DESCRIPTION_ROW}"].value == sheet.description
            header = [
                ws.cell(row=HEADER_ROW, column=i + 1).value
                for i in range(sheet.df.shape[1])
            ]
            assert header == list(sheet.df.columns), f"{sheet.sheet_name}: header {header}"
            assert ws.max_row == HEADER_ROW + len(sheet.df), (
                f"{sheet.sheet_name}: {ws.max_row} rows, expected "
                f"{HEADER_ROW + len(sheet.df)}"
            )

            table = next(iter(ws.tables.values()), None)
            last_col = get_column_letter(sheet.df.shape[1])
            assert table is not None, f"{sheet.sheet_name}: no table"
            assert table.ref == f"A{HEADER_ROW}:{last_col}{HEADER_ROW + len(sheet.df)}", (
                f"{sheet.sheet_name}: table over {table.ref}"
            )

            for offset in {0, len(sheet.df) - 1} - {-1}:
                written = [
                    ws.cell(row=HEADER_ROW + 1 + offset, column=i + 1).value
                    for i in range(sheet.df.shape[1])
                ]
                expected = list(next(_cell_rows(sheet.df.iloc[[offset]])))
                for column, got, want in zip(sheet.df.columns, written, expected):
                    assert _same_cell(got, want), (
                        f"{sheet.sheet_name} data row {offset + 1}, {column}: "
                        f"wrote {want!r}, read back {got!r}"
                    )

        if contents:
            ws = workbook[CONTENTS_SHEET]
            assert [cell.value for cell in ws[1]] == CONTENTS_HEADER
            # Entry n sits on row n + 1, and its description is the sheet's row-1 title minus
            # the `S<n>-<m> ` prefix - i.e. the `title` field, verbatim.
            for sheet in sheets:
                row = [
                    ws.cell(row=sheet.number + 1, column=i + 1).value for i in range(3)
                ]
                assert row == [sheet.number, sheet.sheet_name, sheet.title], row
            assert ws.max_row == len(sheets) + 1, ws.max_row
    finally:
        workbook.close()

    print(
        f"verified: {len(sheets)} sheet(s)"
        + (f" + {CONTENTS_SHEET}, 1-{len(sheets)} consistent with the tabs" if contents else "")
    )


def demo_sheets():
    """Two throwaway sheets carrying the wording pattern the real ones use."""
    return [
        SuppSheet(
            1,
            "S0-1 Demo proteomics",
            "Demonstration table one (related to Fig. 1 and Extended Data Fig. 1)",
            "A stand-in for an assay table: identifier columns, then a per-condition value "
            "block, then statistics. Written only to check the workbook's formatting.",
            pd.DataFrame(
                {
                    "protein": ["AAAA", "BBBB", "CCCC"],
                    "uniprot": ["P00001", "P00002", "P00003"],
                    "LPS vs M0 log2FC": [1.5, -2.25, float("nan")],
                    "LPS vs M0 p_value": [0.001, 0.04, 0.6],
                }
            ),
        ),
        SuppSheet(
            2,
            "S0-2 Demo transcriptomics",
            "Demonstration table two (related to Extended Data Fig. 2)",
            "A second stand-in, so the demo covers the contents index and the per-sheet table "
            "naming as well as one sheet's layout.",
            pd.DataFrame(
                {
                    "gene.symbol": ["Actb", "Gapdh"],
                    "ensembl.gene.id": ["ENSMUSG00000029580", "ENSMUSG00000057666"],
                    "LPS vs M0 log2FoldChange": [0.2, -0.1],
                }
            ),
        ),
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "Data S0 demo.xlsx",
        help="where to write the demo workbook",
    )
    args = parser.parse_args()

    sheets = demo_sheets()
    write_supplementary_workbook(args.out, sheets)
    print(f"wrote {args.out}")
    for sheet in sheets:
        print(
            f"  {sheet.sheet_name!r}: {sheet.df.shape[0]} rows x {sheet.df.shape[1]} columns"
        )


if __name__ == "__main__":
    main()
