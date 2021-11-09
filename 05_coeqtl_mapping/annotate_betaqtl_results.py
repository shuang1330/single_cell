import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import itertools
from tqdm import tqdm


# def annotate_snp_genepair_with_coeqtl_information(coeqtl_path,
#                                                  snp_genepair_path,
#                                                  savepath,
#                                                   alltests_path):
#     snp_genepair = pd.read_csv(snp_genepair_path, sep='\t')
#     snp_genepair['snp_sorted_genepair'] = ['_'.join(item) for item in
#                                            snp_genepair[['snp', 'genepair_sorted']].values]
#     coeqtl = pd.read_csv(coeqtl_path, compression='gzip', sep='\t')
#     coeqtl['snp_genepair'] = ['_'.join(item) for item in coeqtl[['SNPName', 'ProbeName']].values]
#     coeqtl_dict = coeqtl.set_index('snp_genepair')[['AlleleAssessed', 'OverallZScore']].T.to_dict()
#     find_zscore = lambda x:coeqtl_dict.get(x)['OverallZScore'] if x in coeqtl_dict else np.nan
#     snp_genepair['coeqtl_zscore'] = [find_zscore(item) for item in
#                                      snp_genepair['snp_sorted_genepair']]
#     find_allele = lambda x: coeqtl_dict.get(x)['AlleleAssessed'] if x in coeqtl_dict else np.nan
#     snp_genepair['coeqtl_allele'] = [find_allele(item) for item in snp_genepair['snp_sorted_genepair']]
#     alltests_df = pd.read_csv(alltests_path, sep=" ")
#     alltests_set = set(['_'.join(item) for item in alltests_df[['SNPName', 'ProbeName']].values])
#     tested = lambda x: True if x in alltests_set else False
#     snp_genepair['tested'] = [tested(pair) for pair in snp_genepair['snp_sorted_genepair']]
#     snp_genepair[snp_genepair['tested']].to_csv(savepath, sep='\t')
#     return snp_genepair

def argumentsparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--coeqtl', type=str, dest='coeqtl')
    parser.add_argument('--saveprefix', type=str, dest='saveprefix')
    return parser

if __name__ == '__main__':
    args = argumentsparser().parse_args()
    workdir = Path("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping")
    annotated_snp_genepair_path = workdir / 'input/snp_genepair_selection/annotations/UT_monocyte.baseline.annotatedeQTL.individualVarMean.nonZeroRatio.tsv'
    # coeqtl_path = workdir / 'output/summary/coeqtls_betaqtl.tsv'
    # load all coeQTL path
    coeqtl_path = args.coeqtl
    coeqtl = pd.read_csv(coeqtl_path, sep='\t', index_col=0, compression='gzip')
    annotated_snp_genepair_df = pd.read_csv(annotated_snp_genepair_path, sep='\t')
    annotated_snp_genepair_df['snp_genepair'] = ['_'.join(item) for item in
                                                 annotated_snp_genepair_df[['snp', 'genepair_sorted']].values]
    annotated_coeqtl = coeqtl.join(annotated_snp_genepair_df.set_index('snp_genepair'), how='inner')
    # save all coeQTLs results
    annotated_coeqtl.to_csv(f'{args.saveprefix}.annotatedeQTL.individualVarMean.nonZeroRatio.withallcoeqtl.tsv', sep='\t')
    # significant_coeqtl = annotated_coeqtl[annotated_coeqtl['multipletestP'] <= 0.05]
    # save the significant results in both tsv and excel
    significant_coeqtl = annotated_coeqtl[(annotated_coeqtl['snp_qval']<=0.05) & (annotated_coeqtl['gene2_isSig'])]
    significant_coeqtl.to_csv(f'{args.saveprefix}.annotatedeQTL.individualVarMean.nonZeroRatio.tsv', sep='\t')
    significant_coeqtl.to_excel(f'{args.saveprefix}.annotatedeQTL.individualVarMean.nonZeroRatio.xlsx')



