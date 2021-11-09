import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm
import gzip


def update_perm(old_p_list, new_p_list):
    return [np.min([x, y]) for (x, y) in zip(old_p_list, new_p_list)]


def loop_through_one_batch_perm(batch_perm_path, snpgene1_minpvalues_dict, eqtl_dict):
    with gzip.open(batch_perm_path, 'rb') as f:
        f.readline()
        while True:
            line = f.readline().decode('utf-8')
            if not line:
                break
            else:
                linecontent = line.strip().split('\t')
                perm_ps = [float(ele) for ele in linecontent[2:]]
                snp_gene1 = '_'.join([linecontent[1], eqtl_dict.get(linecontent[1])])
                snpgene1_minpvalues_dict[snp_gene1] = update_perm(snpgene1_minpvalues_dict[snp_gene1],
                                                                  perm_ps)
    return snpgene1_minpvalues_dict


def save_numpy(data_df, prefix):
    np.save(f'{prefix}.npy', data_df.values)
    with open(f'{prefix}.cols.txt', 'w') as f:
        f.write('\n'.join([str(ele) for ele in data_df.columns]))
    with open(f'{prefix}.rows.txt', 'w') as f:
        f.write('\n'.join([str(ele) for ele in data_df.index]))
    return None


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--condition', dest='condition')
    parser.add_argument('--celltype', dest='celltype')
    return parser


def loop_through_and_grab_certain_snps(snp_gene1_list, batch_perm_path, snpgene1_minpvalues_dict, eqtl_dict):
    with gzip.open(batch_perm_path, 'rb') as f:
        _ = f.readline()
        while True:
            line = f.readline().decode('utf-8')
            if not line:
                break
            else:
                linecontent = line.strip().split('\t')
                perm_ps = [float(ele) for ele in linecontent[2:]]
                snp_gene1 = '_'.join([linecontent[1], eqtl_dict.get(linecontent[1])])
                if snp_gene1 in snp_gene1_list:
                    snpgene1_minpvalues_dict[snp_gene1].append(perm_ps)
    return snpgene1_minpvalues_dict


def main():
    args = arguments().parse_args()
    condition, celltype = args.condition, args.celltype
    workdir = Path('/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/')
    eqtl_path = workdir/f"input/snp_selection/eqtl/{condition}_{celltype}_eQTLProbesFDR0.05-ProbeLevel.tsv"
    # load eqtl path
    eqtl_df = pd.read_csv(eqtl_path, sep='\t')
    eqtl_df['snp_gene1'] = ['_'.join(item) for item in eqtl_df[['SNPName', 'genename']].values]
    unique_snpgene1 = eqtl_df['snp_gene1'].valuese
    eqtl_dict = eqtl_df.set_index('SNPName')['genename'].T.to_dict()
    # initialize the dict to contain the
    snpgene1_minpvalues_dict = {item: np.ones(100) for item in unique_snpgene1}
    # loop through all batch permutation files
    results_path = workdir/f'output/{condition}_{celltype}'
    for filename in tqdm(os.listdir(results_path / 'noduplicated/output')):
        if '-Permutations.txt.gz' in filename:
            snpgene1_minpvalues_dict = loop_through_one_batch_perm(results_path / 'noduplicated/output'/filename,
                                                                   snpgene1_minpvalues_dict, eqtl_dict)
    for filename in tqdm(os.listdir(results_path / 'duplicatedversion1/output')):
        if '-Permutations.txt.gz' in filename:
            snpgene1_minpvalues_dict = loop_through_one_batch_perm(results_path / 'duplicatedversion1/output'/filename,
                                                                   snpgene1_minpvalues_dict, eqtl_dict)
    for filename in tqdm(os.listdir(results_path / 'duplicatedversion2/output')):
        if '-Permutations.txt.gz' in filename:
            snpgene1_minpvalues_dict = loop_through_one_batch_perm(results_path / 'duplicatedversion2/output'/filename,
                                                                   snpgene1_minpvalues_dict, eqtl_dict)
    snpgene1_minpvalues_df = pd.DataFrame.from_dict(snpgene1_minpvalues_dict)
    save_numpy(snpgene1_minpvalues_df.T, workdir/f'output/{condition}_{celltype}/concated_alltests_permutations')
    return snpgene1_minpvalues_df


if __name__ == '__main__':
    _ = main()
