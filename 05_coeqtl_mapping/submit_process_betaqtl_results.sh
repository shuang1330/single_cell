#!/usr/bin/env bash
#SBATCH --time=02:00:00
#SBATCH --mem=36gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module purge

conda init bash
source /groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/tools/Beeline/miniconda/etc/profile.d/conda.sh
conda activate scpy3.8

condition=$1
celltype=$2

workdir="/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping"
python ${workdir}/output/concat_betaqtl_results.py \
--prefix ${workdir}/output/${condition}_${celltype} \
--savepath ${workdir}/output/${condition}_${celltype}/concated_alltests_output.tsv

python ${workdir}/output/screen_permutation_p_values.py \
--condition ${condition} \
--celltype ${celltype}

python ${workdir}/output/multipletesting_correction.py \
--permutation_prefix ${workdir}/output/${condition}_${celltype}/concated_alltests_permutations \
--coeqtl_path ${workdir}/output/${condition}_${celltype}/concated_alltests_output.tsv \
--eqtl_path ${workdir}/input/snp_selection/eqtl/${condition}_${celltype}_eQTLProbesFDR0.05-ProbeLevel.tsv \
--save_prefix ${workdir}/output/${condition}_${celltype}/coeqtls_fullresults

python ${workdir}/output/annotate_betaqtl_results.py \
--coeqtl ${workdir}/output/${condition}_${celltype}/coeqtls_fullresults.tsv.gz \
--saveprefix ${workdir}/output/${condition}_${celltype}/coeqtls_fullresults