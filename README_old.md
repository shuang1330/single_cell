# Old readme (to be deleted at the end)


## General instructions (for Gearshift)

All required packages are installed in two different conda environments, one for R scritps called `r40` and one for python scripts called `velocyto`

To activate the conda environment on gearshift:

```
source /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/tools/miniconda3/etc/profile.d/conda.sh
conda activate r40
```

All scripts have relative pathes, dependent on the route directory `/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing`. For this reason, the working directory is set in the beginning of all scripts (both R and python) to this directory.

Warning: because some of the scripts were implemented for Calculon, it might happen that the working directory is wrong (`tmp04` instead of `tmp01`) or that a package is missing on conda.

## Evaluation of correlation metrics (Figure 2)

All scripts are saved in `02_correlation_evaluation`:
*  `compare_blueprint_cutoffs.py`: compare correlation between Blueprint and single cell for different expression cutoffs (expressed in x% of the cells); currently only implemented only for Blueprint monocytes and single cell 1 Mio cell v3 (Figure 2a)
*  `figure2_barplot_cutoffs.R`: plotting script for Figure 2a
*  `correlation_timepoint_combined_indivs.py`: calculate correlation per timepoint and cell type for 1 Mio cell data set version 2 and 3, using either a specified expression cutoff or a predefined gene set
*  `figure2_scatterplots.R`: plotting script arranging all the different scatter plots (Figure 2c,d and supplementary figures)
*  `perturbation_wilcoxon.R`: calculate enrichment of correlated genes among KO DE genes usign the Wilcoxon rank sum text (Figure 2e)

## Network construction (Figure 4)

All scripts are saved in `04_network_reconstruction`:

*  `expressed_genes_overlap.R`: evaluate number of expressed genes for each cell type and condition (expressed in 50% of the cells) plus the intersect between them (supplementary figure)
*  `correlation_distribution_cts.R` : check correlation values for each cell type and condition (Figure 4a)
*  `ct_networks_threshold.R` : find optimal correlation cutoff for single cell data by looking at MCC with bulk data (supplementary figure)
* `evaluate_large_networks_compare_edges.R`: compare overlap between edges between cell types and condtions (Figure 4b)
*  `evaluate_large_networks`: evaluate degree distribution (supplementary figure) as well as louvain clusters plus GO term enrichment (Figure 4c) for each network


## Co-Expression QTL Mapping

All scripts are saved in `05_coeqtl_mapping`:

*  `launch_sbatch_files.sh`: pipeline for each step in the coeQTL mapping process
*  `individual_networks.py` : calculate individual networks (zscores) for given dataset, stimulation condition and celltype
*  `merge_coexpression_for_betaeqtl.py` : merge the individual networks from different datasets used for meta analysis for given stimulation condition and celltype
*  `prepare_genelist_and_annotation_for_betaqtl.py`: prepare the input files for coeqtl mapping software
*  `concat_betaqtl_results.py`: merge all output files from the mapping software
*  `screen_permutation_p_values.py`: merge all permutation p values for multiple testing purpose
*  `multipletesting_correction.py`: perform multiple testing correction procedures
*  `annotate_betaqtl_results.py`: merge the coeQTLs results with other information such as non-zero ratio and mean and variances of the co-expression patterns

