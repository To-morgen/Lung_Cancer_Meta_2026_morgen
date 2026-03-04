# Setup (HPC)

This repo keeps code/metadata public while storing large data privately on HPC.

## Create data symlink
```bash
# template
ln -s /path/to/Lung_Cancer_2026 data

# example (HPC)
ln -s /data1/morgen/data_center/Lung_Cancer_2026 data
```

Notes
- Do NOT commit raw/processed data to GitHub.
- All pipelines should read/write via data/ (symlink).