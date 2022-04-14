# 02_correlation_evaluation

*compare_blueprint_cutoffs_CD4T.py* : Compare correlation between Blueprint and single cell (Oelen v3 dataset) for different expression thresholds (number of cells expressing the gene), implemented for UT and CD4+ T cells here

*compare_immunexut_cutoffs_CD4T.py*: Same approach as in *compare_blueprint_cutoffs_CD4T.py*, but comparing correlation between ImmuNexUT and Oelen v3 dataset instead

*figure2_barplot_cutoffs.R*: create barplots from the results of *compare_blueprint_cutoffs_CD4T.py* and *compare_immunexut_cutoffs_CD4T.py*

*normalize_ImmuNexUT.R*: preprocessingImmuNexUT data (separately for each cell type with a matching single-cell cell type) following the description in the corresponding publication (filtering lowly expressed genes, TMM normalization and batch correction) followed by correlation calculation for all genes expressed in 50% of the cells of the Oelen v3 dataset (for comparison with single cell data)
