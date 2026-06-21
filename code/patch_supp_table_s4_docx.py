#!/usr/bin/env python3
"""Patch Supplementary Table S4 in the R2 supplementary docx from R export CSV."""

from __future__ import annotations

import csv
import os
import shutil
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
FINAL_FILES = Path(
    os.environ.get(
        "FINAL_FILES_DIR",
        ROOT / ".." / ".." / "revision" / "revesion_2" / "Final_files",
    )
).resolve()
DOCX = Path(
    os.environ.get(
        "SUPP_DOCX",
        FINAL_FILES / "Supplementary_Information_revision_clean.docx",
    )
)
R2_COPY = ROOT / "text" / "manuscript" / "Supplementary_Information_revision_clean.docx"
CSV = ROOT / "output" / "revision" / "tables" / "supp_table_s4_word_paste.csv"


def find_s4_table(doc: Document):
    for idx, table in enumerate(doc.tables):
        if not table.rows:
            continue
        header = [c.text.strip() for c in table.rows[0].cells]
        if header[:3] == ["Contrast", "Session", "Subscale"]:
            return idx, table
    raise SystemExit("Could not find Supplementary Table S4 in docx")


def main() -> int:
    if not CSV.exists():
        raise SystemExit(f"Missing {CSV}; run: make export-s4")
    if not DOCX.exists():
        raise SystemExit(f"Missing {DOCX}")

    rows = list(csv.DictReader(CSV.open(encoding="utf-8")))
    if len(rows) != 12:
        raise SystemExit(f"Expected 12 S4 data rows, got {len(rows)}")

    backup = DOCX.with_suffix(".docx.bak")
    if not backup.exists():
        shutil.copy2(DOCX, backup)

    doc = Document(DOCX)
    idx, table = find_s4_table(doc)
    if len(table.rows) != 13:
        raise SystemExit(f"S4 table {idx} has {len(table.rows)} rows, expected 13")

    cols = ["Contrast", "Session", "Subscale", "Estimate", "SE", "df", "t-ratio", "p-value"]
    for ri, data in enumerate(rows, start=1):
        for ci, col in enumerate(cols):
            table.rows[ri].cells[ci].text = str(data[col])

    doc.save(DOCX)
    print(f"Patched Table S4 (table index {idx}) in {DOCX}")
    if R2_COPY.parent.exists():
        shutil.copy2(DOCX, R2_COPY)
        print(f"Copied to {R2_COPY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
