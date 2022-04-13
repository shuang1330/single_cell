# 02_correlation_evaluation

*compare_individuals_variance.R* : explore for all genes expressed in at least 50% of the cells the variance across individuals

*correlation_between_celltypes.R* : calculates the Pearson correlation of genepair wise Spearman correlation for all pairwise combinations of cell types within each dataset (for Oelen dataset (V2) and (V3)), taking only genes expressed in 50% of the cells in both cell types; plots results in heatmap afterwards

*correlation_celltype.py* : calculates Spearman correlation for each genepair expressed in 50% of the cells for Oelen dataset (V2) and (V3), separately per cell type, but combing all individuals; provides so the input csv files for *correlation_between_celltypes.R*

*correlation_correlation_distribution_celltypes_and_individuals.R* :

*correlation_subsampling.py* :

*plot_indiv_subsampling_effect.R* :
