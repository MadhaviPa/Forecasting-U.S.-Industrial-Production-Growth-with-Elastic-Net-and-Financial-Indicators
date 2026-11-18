# Results and provenance

The three CSV files in this folder are **transcriptions of published tables** from the uploaded seminar report, not newly executed model outputs. They reflect the rounded numbers reported in the PDF.

- `msfe_results.csv`: report Table 2, 180-month rolling window.
- `dm_results.csv`: report Table 3, HLN-adjusted point-forecast DM tests.
- `quantile_results.csv`: report Table 5, best quantile model by horizon and quantile (not every candidate model).

Running `Rscript run_all.R` creates the corresponding full recomputed outputs in `results/generated/`. The analysis also writes detailed diagnostics under `results/details/`. Recomputed values may differ if data vintages or upstream series change.
