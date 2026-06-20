#!/usr/bin/env python3
"""
Extract numeric claims from R2 submission docx files.

Usage:
  python code/extract_manuscript_claims.py
  python code/extract_manuscript_claims.py --check   # regen and report counts only

Writes docs/extracted_claims_r2.json and docs/CLAIM_COUNTS.md
"""

from __future__ import annotations

import argparse
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

STAT_PATTERNS = [
    ("F_test", re.compile(r"F\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*([\d.]+)\s*,?\s*p\s*([<.=]+)\s*([\d.]+)", re.I)),
    ("chi2", re.compile(r"χ\s*²\s*\(\s*(\d+)\s*\)\s*=\s*([\d.]+)\s*,?\s*p\s*=\s*([\d.]+)", re.I)),
    ("pearson", re.compile(r"Pearson r\s*=\s*(-?[\d.]+)\s*,?\s*p\s*=\s*([\d.]+)", re.I)),
    ("spearman", re.compile(r"(?:Spearman )?(?:ρ|rho)\s*=\s*([\d.]+)\s*,?\s*p\s*([<.=]+)\s*([\d.]+)", re.I)),
    ("mean_sd", re.compile(r"(\d+\.\d+)\s*±\s*(\d+\.\d+)")),
    ("t_test", re.compile(r"t\s*\(\s*(\d+)\s*\)\s*=\s*([\d.]+)\s*,?\s*p\s*([<.=]+)\s*([\d.]+)", re.I)),
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def docx_paths(root: Path) -> dict[str, Path]:
    r2 = root / "text" / "npp_revision_2026_r2"
    return {
        "main": r2 / "Main_R2_tc.docx",
        "supplementary": r2 / "Supplementary_Information_revision_clean.docx",
        "rebuttal": r2 / "Rebuttal_R2.docx",
    }


def paragraph_text(p_el) -> str:
    parts = []
    for t in p_el.iter(W + "t"):
        if t.text:
            parts.append(t.text)
    return "".join(parts).strip()


def table_to_rows(tbl_el) -> list[list[str]]:
    rows = []
    for tr in tbl_el.findall("w:tr", NS):
        cells = []
        for tc in tr.findall("w:tc", NS):
            texts = [t.text or "" for t in tc.iter(W + "t")]
            cells.append("".join(texts).strip())
        if any(cells):
            rows.append(cells)
    return rows


def merge_table_rows(existing: list[list[str]], new_rows: list[list[str]]) -> list[list[str]]:
    """Append rows when Word splits one supplementary table across multiple tbl elements."""
    if not existing:
        return new_rows
    if not new_rows:
        return existing
    if existing[0] == new_rows[0]:
        return existing + new_rows[1:]
    return existing + new_rows


def parse_docx_tables(path: Path) -> tuple[str, list[list[list[str]]], dict[str, list[list[str]]]]:
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml")
    root = ET.fromstring(xml)
    body = root.find("w:body", NS)
    full_text_parts: list[str] = []
    all_tables: list[list[list[str]]] = []
    supp_tables: dict[str, list[list[str]]] = {}
    pending_caption: str | None = None
    last_unassigned: list[list[str]] | None = None

    for child in body:
        tag = child.tag
        if tag == W + "p":
            txt = paragraph_text(child)
            if txt:
                full_text_parts.append(txt)
            m = re.search(r"Supplementary Table S(\d+)", txt, re.I)
            if m:
                sid = f"S{m.group(1)}"
                if last_unassigned is not None:
                    # Caption follows table (common in this supplementary docx).
                    supp_tables[sid] = last_unassigned
                    last_unassigned = None
                    pending_caption = None
                elif sid not in supp_tables:
                    pending_caption = sid
                # Else: duplicate caption line for an already-assigned table — ignore.
        elif tag == W + "tbl":
            rows = table_to_rows(child)
            all_tables.append(rows)
            if pending_caption:
                if pending_caption in supp_tables:
                    supp_tables[pending_caption] = merge_table_rows(supp_tables[pending_caption], rows)
                else:
                    supp_tables[pending_caption] = rows
                pending_caption = None
                last_unassigned = None
            else:
                last_unassigned = rows

    return "\n".join(full_text_parts), all_tables, supp_tables


def extract_in_text(text: str, source: str) -> list[dict]:
    claims = []
    for kind, pat in STAT_PATTERNS:
        for m in pat.finditer(text):
            claims.append({
                "source": source,
                "kind": kind,
                "raw": m.group(0),
                "groups": list(m.groups()),
            })
    return claims


def numeric_cells_in_table(rows: list[list[str]]) -> list[dict]:
    out = []
    for ri, row in enumerate(rows):
        for ci, cell in enumerate(row):
            cell = cell.strip()
            if not cell or cell in ("***", "ns", "NA"):
                continue
            # Skip label cells that merely contain digits (e.g. "S1 (R)").
            if re.search(r"[A-Za-z]", cell):
                continue
            if re.match(r"^[\d.+\-eE<>\s]+$", cell.replace("±", "").replace(",", "")) or re.search(r"\d", cell):
                out.append({"row": ri, "col": ci, "raw": cell})
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    root = repo_root()
    paths = docx_paths(root)
    payload: dict = {"documents": {}, "supplementary_tables": {}, "summary": {}}

    total_in_text = 0
    total_cells = 0

    for name, path in paths.items():
        if not path.exists():
            raise SystemExit(f"Missing docx: {path}")
        text, all_tables, supp_tables = parse_docx_tables(path)
        in_text = extract_in_text(text, name)
        total_in_text += len(in_text)
        payload["documents"][name] = {
            "path": str(path.relative_to(root)),
            "in_text": in_text,
            "n_tables": len(all_tables),
        }
        if name == "supplementary":
            for sid, rows in supp_tables.items():
                cells = numeric_cells_in_table(rows)
                total_cells += len(cells)
                payload["supplementary_tables"][sid] = {
                    "n_rows": len(rows),
                    "n_numeric_cells": len(cells),
                    "rows": rows,
                    "numeric_cells": cells,
                }
            # Export S1 connectivity ROI atlas as tidy CSV (label, x, y, z, network).
            s1_path = root / "data" / "connectivity_roi_list.csv"
            if "S1" in supp_tables:
                rows = supp_tables["S1"]
                if len(rows) > 2:
                    out_lines = ["roi_label,x,y,z,network"]
                    for row in rows[2:]:
                        if len(row) < 5 or not row[0] or row[0] == "ROI Label":
                            continue
                        label = row[0].replace(",", ";")
                        out_lines.append(",".join([label, row[1], row[2], row[3], row[4]]))
                    s1_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")

    payload["summary"] = {
        "in_text_claims": total_in_text,
        "supp_table_ids": sorted(payload["supplementary_tables"].keys()),
        "supp_numeric_cells": total_cells,
    }

    out_json = root / "docs" / "extracted_claims_r2.json"
    out_md = root / "docs" / "CLAIM_COUNTS.md"

    if not args.check:
        out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md = [
        "# R2 claim extraction summary",
        "",
        f"- In-text stat patterns: **{total_in_text}**",
        f"- Supplementary table IDs: **{', '.join(payload['summary']['supp_table_ids'])}**",
        f"- Supplementary numeric cells: **{total_cells}**",
        "",
    ]
    for sid, meta in sorted(payload["supplementary_tables"].items()):
        md.append(f"- Table {sid}: {meta['n_rows']} rows, {meta['n_numeric_cells']} numeric cells")
    out_md.write_text("\n".join(md) + "\n", encoding="utf-8")

    print(json.dumps(payload["summary"], indent=2))
    if not args.check:
        print(f"Wrote {out_json}")
    print(f"Wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
