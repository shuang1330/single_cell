# 01_association_metrics


In *setting_files_for_grnboost2* there are setting files for running GRNBoost2 on scRNAseq and BIOS data

*GRNBoost2.ipynb*: Prepare the input files for BEELINE, and examine the GRNBoost2 results for scRNAseq and for BIOS data

*scorpius_and_slingshot_clean.R*: calculates the pseudotime ordering for Oelen v2 classical monocytes, using SCORPIUS and Slingshot algorithms

*scvelo_analysis_dm.py*: runs RNA velocity analysis on Oelen v3 dataset classical monocytes after creating loom files using [velocyto](http://velocyto.org/velocyto.py/tutorial/cli.html) to get both spliced and unspliced gene count matrices

*compare_cell_classification.ipynb*: compares the aximuth cell type classification with the marker gene cell type classification in Oelen v2 and v3 dataset for untreated cells
