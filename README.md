# MN_data_pipeline

Mascot data pipeline for membranous nephropathy with categorization based on scores.

This project processes Mascot score data from Excel and CSV files, merges them, deduplicates by accession, and generates various outputs including mappings and visualizations.

## Features

- Merge Excel/CSV files (all sheets)
- Deduplicate by accession (keep max Score)
- Produce two sheets (by_score, by_score_summary)
- Map every source sheet to master accessions (score columns only)
- Produce long mapping + summary
- Visualize Score (heatmap top N, bar plot sums)

## Outputs

- `merged_all_rows.xlsx`: All merged rows
- `merged_unique_by_score.xlsx`: Deduplicated by accession
- `final_two_sheets.xlsx`: By score and summary sheets
- `final_master_with_sheet_mapping.xlsx`: Wide mapping with sample scores
- `accession_sheet_score_map.xlsx`: Long mapping and summary
- `heatmap_score_top50.png`: Heatmap of top 50 accessions
- `score_sum_per_sample.png`: Bar plot of score sums per sample

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/abbhinayy7/MN_data_pipeline.git
   cd MN_data_pipeline
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. Place your Excel (.xlsx, .xls) and CSV files in the root directory.
2. Run the pipeline:
   ```bash
   python src/mascot_pipeline.py
   ```
3. Outputs will be generated in the root directory.

## Requirements

- Python 3.6+
- pandas
- openpyxl
- matplotlib
- numpy

## License

[Add license if applicable]
