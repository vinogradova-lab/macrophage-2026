"""Reusable writer for manuscript supplementary-data workbooks ("Data S*")."""

from dataclasses import dataclass

from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter
import pandas as pd

TITLE_ROW = 1
DESCRIPTION_ROW = 2
HEADER_ROW = 3  # 1-indexed Excel row; == pandas startrow (2) + 1


@dataclass
class SuppSheet:
    """One sheet in a supplementary workbook. ``number`` is its 1-based contents position;
    ``sheet_name``'s first token prefixes the row-1 title; ``title`` is the row-1 body (reused as
    the contents Description); ``description`` is the row-2 blurb; ``df`` is the flat table."""

    number: int
    sheet_name: str
    title: str
    description: str
    df: pd.DataFrame


def _decorate_data_sheet(ws, n_cols, title, description, freeze_cols):
    """Add the title/description rows, bold header, and freeze panes to a written sheet."""
    last_col = get_column_letter(n_cols)

    ws.merge_cells(f"A{TITLE_ROW}:{last_col}{TITLE_ROW}")
    title_cell = ws.cell(row=TITLE_ROW, column=1, value=title)
    title_cell.font = Font(bold=True, size=13)
    title_cell.alignment = Alignment(vertical="center")

    ws.merge_cells(f"A{DESCRIPTION_ROW}:{last_col}{DESCRIPTION_ROW}")
    desc_cell = ws.cell(row=DESCRIPTION_ROW, column=1, value=description)
    desc_cell.font = Font(italic=True, size=10)
    desc_cell.alignment = Alignment(wrap_text=True, vertical="top")

    ws.row_dimensions[TITLE_ROW].height = 20
    ws.row_dimensions[DESCRIPTION_ROW].height = 60

    for cell in ws[HEADER_ROW]:
        cell.font = Font(bold=True)

    # AutoFilter dropdowns on the header row, matching the sibling repo's workbooks.
    ws.auto_filter.ref = f"A{HEADER_ROW}:{last_col}{ws.max_row}"

    ws.freeze_panes = f"{get_column_letter(freeze_cols + 1)}{HEADER_ROW + 1}"


def _write_contents_sheet(writer, sheets, contents_name):
    contents = pd.DataFrame(
        [
            {"#": s.number, "Title": s.sheet_name, "Description": s.title}
            for s in sheets
        ]
    )
    contents.to_excel(writer, sheet_name=contents_name, index=False)
    ws = writer.sheets[contents_name]
    for cell in ws[1]:
        cell.font = Font(bold=True)
    ws.column_dimensions["A"].width = 4
    ws.column_dimensions["B"].width = 30
    ws.column_dimensions["C"].width = 100
    for row in ws.iter_rows(min_row=2, min_col=3, max_col=3):
        row[0].alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "A2"


def write_supplementary_workbook(
    path, sheets, contents_name="contents", freeze_cols=2
):
    """Write a multi-sheet workbook: a ``contents`` index + one sheet per table."""
    with pd.ExcelWriter(path, engine="openpyxl", mode="w") as writer:
        _write_contents_sheet(writer, sheets, contents_name)
        for s in sheets:
            prefix = s.sheet_name.split(" ", 1)[0]
            row1_title = f"{prefix} {s.title}"
            s.df.to_excel(
                writer,
                sheet_name=s.sheet_name,
                index=False,
                startrow=HEADER_ROW - 1,
            )
            _decorate_data_sheet(
                writer.sheets[s.sheet_name],
                s.df.shape[1],
                row1_title,
                s.description,
                freeze_cols,
            )
    return path


def write_supplementary_sheet(
    df,
    path,
    sheet_name,
    title,
    description,
    index=False,
    mode="w",
    freeze_cols=2,
):
    """Write one titled sheet. ``mode="w"`` creates/overwrites; ``mode="a"`` appends,
    replacing a same-named sheet."""
    kwargs = {}
    if mode == "a":
        kwargs["if_sheet_exists"] = "replace"

    with pd.ExcelWriter(path, engine="openpyxl", mode=mode, **kwargs) as writer:
        df.to_excel(writer, sheet_name=sheet_name, index=index, startrow=HEADER_ROW - 1)
        n_cols = df.shape[1] + (df.index.nlevels if index else 0)
        _decorate_data_sheet(
            writer.sheets[sheet_name], n_cols, title, description, freeze_cols
        )
    return path
