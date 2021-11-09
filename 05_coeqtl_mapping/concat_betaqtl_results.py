import pandas as pd
from pathlib import Path
import os
import argparse
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests


def concat_results(prefix, savepath):
    concated_df = pd.DataFrame()
    for filename in tqdm(os.listdir(prefix/'noduplicated/output')):
        if filename.endswith("-TopEffects.txt"):
            df = pd.read_csv(prefix/'noduplicated/output'/filename, sep='\t')
            concated_df = pd.concat([concated_df, df], axis=0)
    concated_df['snp_genepair'] = ['_'.join(item) for item in concated_df[['SNP', 'Gene']].values]
    version1 = pd.DataFrame()
    for filename in tqdm(os.listdir(prefix/'duplicatedversion1/output')):
        if filename.endswith("-TopEffects.txt"):
            df = pd.read_csv(prefix/'duplicatedversion1/output'/filename, sep='\t')
            version1 = pd.concat([version1, df], axis=0)
    version1['snp_genepair'] = ['_'.join(item) for item in version1[['SNP', 'Gene']].values]
    version2 = pd.DataFrame()
    for filename in tqdm(os.listdir(prefix/'duplicatedversion2/output')):
        if filename.endswith("-TopEffects.txt"):
            df = pd.read_csv(prefix/'duplicatedversion2/output'/filename, sep='\t')
            version2 = pd.concat([version2, df], axis=0)
    version2['snp_genepair'] = ['_'.join(item) for item in version2[['SNP', 'Gene']].values]
    concated_versions = pd.concat([concated_df, version1, version2], axis=0)
    concated_versions = concated_versions.sort_values(by=['GeneChr', 'GenePos'])
    concated_versions = concated_versions.set_index('snp_genepair')
    # add multiple test significance
    concated_versions['multipletestP'] = multipletests(concated_versions['BetaAdjustedMetaP'],
                                                       alpha=0.05, method='fdr_bh',
                                                       is_sorted=False, returnsorted=False)[1]
    concated_versions.to_csv(savepath, sep='\t')
    return concated_versions

def argumentsparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=str, dest='prefix')
    parser.add_argument('--savepath', type=str, dest='savepath')
    return parser

if __name__ == '__main__':
    args = argumentsparser().parse_args()
    concat_results(Path(args.prefix), args.savepath)