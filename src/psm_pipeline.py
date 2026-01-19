#!/usr/bin/env python3
"""
Full PSM pipeline:
- Merge Excel/CSV files (all sheets)
- Deduplicate by accession (keep max PSM)
- Produce two sheets (by_psm, by_psm_summary)
- Map every source sheet (filename__sheet) to master accessions (psm columns only)
- Produce long mapping + summary
- Visualize PSM (heatmap top N, bar plot sums)
Outputs:
  merged_all_rows.xlsx
  merged_unique_by_psm.xlsx
  final_two_sheets.xlsx
  final_master_with_sheet_mapping.xlsx
  accession_sheet_psm_map.xlsx
  heatmap_psm_top50.png
  psm_sum_per_sample.png
"""
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import logging
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ---------------- CONFIG ----------------
OUT_MERGED_ALL = "merged_all_rows.xlsx"
OUT_DEDUP = "merged_unique_by_psm.xlsx"
OUT_TWO = "final_two_sheets.xlsx"
OUT_WIDE = "final_master_with_sheet_mapping.xlsx"
OUT_MAP = "accession_sheet_psm_map.xlsx"
HEATMAP_PNG = "heatmap_psm_top50.png"
BAR_PNG = "psm_sum_per_sample.png"
TOP_N = 50  # top N accessions in heatmap
# files to exclude from scanning (output names)
EXCLUDE_NAMES = {
    OUT_MERGED_ALL,
    OUT_DEDUP,
    OUT_TWO,
    OUT_WIDE,
    OUT_MAP,
    HEATMAP_PNG,
    BAR_PNG,
    "psm_pipeline.py"  # script name
}

# ---------------- helpers ----------------
def find_best_column(cols, candidates):
    """Find the best matching column from candidates."""
    cols_list = list(cols)
    lowers = [c.lower() for c in cols_list]
    for cand in candidates:
        c_low = cand.lower()
        # exact or substring match
        for orig, low in zip(cols_list, lowers):
            if c_low == low or c_low in low or low in c_low:
                return orig
    return None

def safe_read_excel_sheets(path: Path):
    """Read all sheets from Excel file and return list of (sheetname, df)."""
    out = []
    try:
        xls = pd.ExcelFile(path)
    except Exception as e:
        logging.error(f"Could not open Excel {path}: {e}")
        return out
    for sheet in xls.sheet_names:
        try:
            df = pd.read_excel(xls, sheet_name=sheet, dtype=str)
            if df is None or df.shape[0] == 0:
                continue
            out.append((sheet, df))
        except Exception as e:
            logging.error(f"Error reading {path.name}::{sheet}: {e}")
    return out

def read_csv_file(path: Path):
    """Read CSV file."""
    try:
        df = pd.read_csv(path, dtype=str)
    except Exception:
        try:
            df = pd.read_csv(path, dtype=str, engine="python")
        except Exception as e:
            logging.error(f"Could not read CSV {path.name}: {e}")
            return None
    return df

def merge_files(files):
    """Merge all files into a single dataframe."""
    dfs = []
    for f in files:
        logging.info(f"Scanning file: {f.name}")
        if f.suffix.lower() in [".xlsx", ".xls"]:
            sheets = safe_read_excel_sheets(f)
            for sheetname, df in sheets:
                acc_col = find_best_column(df.columns, ["accession", "accessions", "protein accession", "protein id", "id"])
                df["_source_file"] = f.name
                df["_source_sheet"] = sheetname
                if acc_col:
                    df = df.rename(columns={acc_col: "accession"})
                dfs.append(df)
        else:  # csv
            df = read_csv_file(f)
            if df is not None and df.shape[0] > 0:
                acc_col = find_best_column(df.columns, ["accession", "accessions", "protein accession", "protein id", "id"])
                df["_source_file"] = f.name
                df["_source_sheet"] = "csv"
                if acc_col:
                    df = df.rename(columns={acc_col: "accession"})
                dfs.append(df)
    if not dfs:
        logging.error("No data found in scanned files.")
        sys.exit(1)
    merged = pd.concat(dfs, ignore_index=True, sort=False)
    return merged

def deduplicate_by_psm(merged):
    """Deduplicate by accession, keeping max PSM."""
    cols = merged.columns.tolist()
    acc_col = find_best_column(cols, ["accession", "accessions", "protein accession", "protein id", "id"])
    desc_col = find_best_column(cols, ["description", "protein name", "protein description", "names"])
    psm_col = find_best_column(cols, ["psm", "psms", "peptide spectrum matches", "psm count", "peptide count", "num. of matches", "num of matches", "number of matches"])

    if acc_col is None:
        logging.error("No accession column detected in merged data.")
        sys.exit(1)

    # Normalize canonical names
    if acc_col != "accession":
        merged = merged.rename(columns={acc_col: "accession"})
    if desc_col and desc_col != "description":
        merged = merged.rename(columns={desc_col: "description"})
    if psm_col and psm_col != "psm":
        merged = merged.rename(columns={psm_col: "psm"})

    # Ensure types
    merged["accession"] = merged["accession"].astype(str).str.strip()
    if "psm" in merged.columns:
        merged["psm"] = pd.to_numeric(merged["psm"], errors="coerce")
    else:
        merged["psm"] = pd.NA

    # Deduplicate
    if merged["psm"].notna().any():
        merged["_psm_fill"] = merged["psm"].fillna(float("-inf"))
        dedup = merged.sort_values(["accession", "_psm_fill"], ascending=[True, False]).drop_duplicates(subset=["accession"], keep="first")
        dedup = dedup.drop(columns=["_psm_fill"])
    else:
        dedup = merged.drop_duplicates(subset=["accession"], keep="first")
    return dedup

def create_two_sheets(dedup):
    """Create final_two_sheets.xlsx with by_psm and by_psm_summary."""
    by_psm_cols = ["accession"]
    if "description" in dedup.columns:
        by_psm_cols.append("description")
    if "psm" in dedup.columns:
        by_psm_cols.append("psm")
    by_psm = dedup.loc[:, [c for c in by_psm_cols if c in dedup.columns]].copy()
    if "psm" in by_psm.columns:
        by_psm = by_psm.sort_values("psm", ascending=False)
    by_psm_summary = by_psm.copy()
    return by_psm, by_psm_summary

def simplify_sheet_name(sheetname: str) -> str:
    """Simplify sheet name to a clean identifier."""
    import re
    match = re.search(r'(VB\d+)', sheetname.upper())
    if match:
        return match.group(1)
    return sheetname.split(",")[0].split("(")[0].strip().replace(" ", "_")

def create_wide_mapping(master, files):
    """Create wide mapping with psm columns for each sample."""
    out_wide = master.copy()
    out_wide["accession"] = out_wide["accession"].astype(str).str.strip()
    sample_psm_cols = []

    for f in files:
        if f.suffix.lower() in [".xlsx", ".xls"]:
            sheets = safe_read_excel_sheets(f)
            for sheetname, df in sheets:
                acc_col_local = find_best_column(df.columns, ["accession", "accessions", "protein accession", "protein id", "id"])
                if acc_col_local is None:
                    continue
                df = df.rename(columns={acc_col_local: "accession"})
                df["accession"] = df["accession"].astype(str).str.strip()
                psm_local = find_best_column(df.columns, ["psm", "psms", "peptide spectrum matches", "psm count", "peptide count", "num. of matches", "num of matches", "number of matches"])
                clean_sheet = simplify_sheet_name(sheetname)
                psm_colname = f"{clean_sheet}_psm"
                if psm_colname not in sample_psm_cols:
                    sample_psm_cols.append(psm_colname)
                if psm_local and psm_local in df.columns:
                    mapping_psm = df.set_index("accession")[psm_local].apply(lambda x: pd.to_numeric(x, errors="coerce"))
                    out_wide[psm_colname] = out_wide["accession"].map(lambda a: mapping_psm.get(str(a), np.nan)).fillna(0)
                else:
                    out_wide[psm_colname] = 0
        else:  # CSV
            df = read_csv_file(f)
            if df is None:
                continue
            acc_col_local = find_best_column(df.columns, ["accession", "accessions", "protein accession", "protein id", "id"])
            if acc_col_local is None:
                continue
            df = df.rename(columns={acc_col_local: "accession"})
            df["accession"] = df["accession"].astype(str).str.strip()
            psm_local = find_best_column(df.columns, ["psm", "psms", "peptide spectrum matches", "psm count", "peptide count", "num. of matches", "num of matches", "number of matches"])
            psm_colname = "CSV_psm"
            if psm_colname not in sample_psm_cols:
                sample_psm_cols.append(psm_colname)
            if psm_local and psm_local in df.columns:
                mapping_psm = df.set_index("accession")[psm_local].apply(lambda x: pd.to_numeric(x, errors="coerce"))
                out_wide[psm_colname] = out_wide["accession"].map(lambda a: mapping_psm.get(str(a), np.nan)).fillna(0)
            else:
                out_wide[psm_colname] = 0

    # Rename key columns
    rename_dict = {"accession": "Accession", "description": "Description", "psm": "PSMs"}
    out_wide = out_wide.rename(columns={k: v for k, v in rename_dict.items() if k in out_wide.columns})

    # Clean Description
    if "Description" in out_wide.columns:
        out_wide["Description"] = out_wide["Description"].astype(str).str.replace(r"\s+OS=.*", "", regex=True).str.strip()

    # Compute counts
    existing_sample_cols = [c for c in sample_psm_cols if c in out_wide.columns]
    if existing_sample_cols:
        numeric_matrix = out_wide[existing_sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        out_wide["Count_in_samples"] = (numeric_matrix > 0).sum(axis=1).astype(int)
    else:
        out_wide["Count_in_samples"] = 0

    # Reorder columns
    cols_order = []
    for key in ["Accession", "Description", "PSMs"]:
        if key in out_wide.columns:
            cols_order.append(key)
    cols_order += ["Count_in_samples"]
    cols_order += [c for c in out_wide.columns if c not in cols_order]
    out_wide = out_wide[cols_order]

    return out_wide, sample_psm_cols

def create_long_mapping(merged, master, dedup):
    """Create long mapping and summary."""
    merged_for_map = merged.copy()
    psm_col_map = find_best_column(merged_for_map.columns, ["psm", "psms", "peptide spectrum matches", "psm count", "peptide count", "num. of matches", "num of matches", "number of matches"])
    if psm_col_map and psm_col_map != "psm_in_sheet":
        merged_for_map = merged_for_map.rename(columns={psm_col_map: "psm_in_sheet"})
    elif "psm" in merged_for_map.columns:
        merged_for_map = merged_for_map.rename(columns={"psm": "psm_in_sheet"})
    else:
        merged_for_map["psm_in_sheet"] = pd.NA

    merged_for_map["accession"] = merged_for_map["accession"].astype(str)
    master_set = set(master["accession"].astype(str))
    long_map = merged_for_map[merged_for_map["accession"].isin(master_set)].copy()
    dedup_lookup = dedup.set_index("accession")
    long_map = long_map.merge(dedup_lookup[["description"]].reset_index(), on="accession", how="left")
    cols_long = ["accession", "description", "_source_file", "_source_sheet", "psm_in_sheet"]
    long_map_out = long_map[[c for c in cols_long if c in long_map.columns]]

    def build_summary(g):
        items = []
        for _, r in g.iterrows():
            sf = r.get("_source_file", "")
            ss = r.get("_source_sheet", "")
            sc = r.get("psm_in_sheet", "")
            items.append(f"{sf}::{ss} (psm={sc})")
        seen = set()
        out = []
        for it in items:
            if it not in seen:
                out.append(it)
                seen.add(it)
        return "; ".join(out)

    summary = long_map_out.groupby("accession").apply(build_summary).reset_index().rename(columns={0: "sheet_psm_list"})
    summary = dedup[["accession", "description"]].merge(summary, on="accession", how="left")
    summary["sheet_psm_list"] = summary["sheet_psm_list"].fillna("")
    return long_map_out, summary

def visualize(out_wide, sample_psm_cols, top_n, heatmap_png, bar_png):
    """Create visualizations."""
    # Heatmap top N
    if "PSMs" in out_wide.columns:
        top_accessions = out_wide.sort_values("PSMs", ascending=False).head(top_n)["Accession"]
    else:
        top_accessions = out_wide.head(top_n)["Accession"]

    heatmap_data = out_wide[out_wide["Accession"].isin(top_accessions)][sample_psm_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    heatmap_data.index = out_wide[out_wide["Accession"].isin(top_accessions)]["Accession"]

    plt.figure(figsize=(10, 8))
    plt.imshow(heatmap_data.values, aspect='auto', cmap='viridis')
    plt.colorbar(label='PSMs')
    plt.xticks(range(len(sample_psm_cols)), sample_psm_cols, rotation=90)
    plt.yticks(range(len(heatmap_data)), heatmap_data.index)
    plt.title(f'Heatmap of Top {top_n} Accessions PSMs')
    plt.tight_layout()
    plt.savefig(heatmap_png)
    plt.close()
    logging.info(f"Saved heatmap -> {heatmap_png}")

    # Bar plot sums
    sums = out_wide[sample_psm_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum()
    plt.figure(figsize=(10, 6))
    sums.plot(kind='bar')
    plt.title('Sum of PSMs per Sample')
    plt.ylabel('Total PSMs')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(bar_png)
    plt.close()
    logging.info(f"Saved bar plot -> {bar_png}")

def main():
    cwd = Path(".")
    files = [f for f in cwd.glob("*") if f.suffix.lower() in [".xlsx", ".xls", ".csv"] and f.name not in EXCLUDE_NAMES]
    if not files:
        logging.error("No Excel/CSV files found in current directory. Place your files here and rerun.")
        sys.exit(1)

    logging.info(f"Files to scan: {[f.name for f in files]}")

    # 1) Merge
    merged = merge_files(files)
    merged.to_excel(OUT_MERGED_ALL, index=False)
    logging.info(f"Saved merged all rows -> {OUT_MERGED_ALL} (rows: {len(merged)})")

    # 2) Dedup
    dedup = deduplicate_by_psm(merged)
    dedup.to_excel(OUT_DEDUP, index=False)
    logging.info(f"Saved deduplicated file -> {OUT_DEDUP} (unique accessions: {dedup['accession'].nunique()})")

    # 3) Two sheets
    by_psm, by_psm_summary = create_two_sheets(dedup)
    with pd.ExcelWriter(OUT_TWO, engine="openpyxl") as writer:
        by_psm.to_excel(writer, sheet_name="by_psm", index=False)
        by_psm_summary.to_excel(writer, sheet_name="by_psm_summary", index=False)
    logging.info(f"Saved two-sheet master -> {OUT_TWO}")

    # 4) Wide mapping
    master = by_psm.copy()
    out_wide, sample_psm_cols = create_wide_mapping(master, files)
    out_wide.to_excel(OUT_WIDE, index=False)
    logging.info(f"Saved wide mapping -> {OUT_WIDE}")

    # 5) Long mapping
    long_map_out, summary = create_long_mapping(merged, master, dedup)
    with pd.ExcelWriter(OUT_MAP, engine="openpyxl") as writer:
        long_map_out.to_excel(writer, sheet_name="long_mapping", index=False)
        summary.to_excel(writer, sheet_name="summary_per_accession", index=False)
    logging.info(f"Saved mapping workbook -> {OUT_MAP} (long rows: {len(long_map_out)}, summary rows: {len(summary)})")

    # 6) Visualize
    visualize(out_wide, sample_psm_cols, TOP_N, HEATMAP_PNG, BAR_PNG)

    logging.info("PSM pipeline finished successfully.")

if __name__ == "__main__":
    main()