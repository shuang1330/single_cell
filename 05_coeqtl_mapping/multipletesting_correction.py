import pandas as pd
from statsmodels.stats.multitest import multipletests
import numpy as np
import argparse
from scipy.optimize import minimize
from scipy.stats import beta
from scipy import special
from pathlib import Path


def read_numpy(prefix):
    data = np.load(f'{prefix}.npy')
    columns = [f'perm{item.strip()}' for item in open(f'{prefix}.cols.txt', 'r').readlines()]
    rows = [item.strip() for item in open(f'{prefix}.rows.txt', 'r').readlines()]
    return pd.DataFrame(data=data, columns=columns, index=rows)


def beta_distribution_mle_function(x, p):
    k, n = x
    ll = (k - 1) * np.sum(np.log(p)) + (n - 1) * np.sum(np.log(1 - p)) - np.size(p) * special.betaln(k, n)
    return -1 * ll


def beta_distribution_initial_guess(x):
    """
    https://stats.stackexchange.com/questions/13245/which-is-a-good-tool-to-compute-parameters-for-a-beta-distribution
    """
    mean = np.mean(x)
    var = np.var(x)
    a = mean * ((mean * (1 - mean) / var) - 1)
    b = (1 - mean) * ((mean * (1 - mean) / var) - 1)
    return a, b


def fit_beta_distribution(p, a_bnd=(0.1, 10), b_bnd=(1, 1000000)):
    a, b = beta_distribution_initial_guess(p)
    x0 = np.array([min(max(a, a_bnd[0]), a_bnd[1]), min(max(b, b_bnd[0]), b_bnd[1])])
    res = minimize(beta_distribution_mle_function,
                   x0=x0,
                   args=(p, ),
                   method='nelder-mead',
                   bounds=(a_bnd, b_bnd),
                   options={"maxiter": 10000, "disp": True})
    return res.x, res.nfev, res.nit


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--permutation_prefix', dest='permutation_pvalue_prefix')
    parser.add_argument('--coeqtl_path', dest='coeqtl_path')
    parser.add_argument('--eqtl_path', dest='eqtl_path')
    parser.add_argument('--save_prefix', dest='saveprefix')
    return parser


def main():
    # workdir = Path("/groups/umcg-lld/tmp01/projects/1MCellRNAseq/GRN_reconstruction/ongoing/coeqtl_mapping/")
    # # fit beta distribution per each eQTL
    # permutation_pvalue_prefix = workdir/'output/UT_monocytes/concated_alltests_permutations'
    args = arguments().parse_args()
    permutation_pvalue_prefix = args.permutation_pvalue_prefix
    permutation_pvalues_df = read_numpy(permutation_pvalue_prefix)
    permutation_cols = [f'perm{ind}' for ind in range(0, 100)]
    permutation_pvalues_df['beta_shape1'], permutation_pvalues_df['beta_shape2'] = \
    zip(*[fit_beta_distribution(x)[0] for x in permutation_pvalues_df[permutation_cols].values])
    # for snp_gene in permutation_pvalues_df.index:
    #     perms_p = permutation_pvalues_df[permutation_cols].loc[snp_gene]
    #     print(snp_gene)
    #     _ = beta.fit(perms_p)[:2]
    # load the merged results and find the lowest nominalP per each eQTL
    # coeqtl_path = workdir/"output/UT_monocytes/concated_alltests_output.tsv"
    coeqtl_path = args.coeqtl_path
    # first find the eqtl gene
    # eqtls_path = workdir/'input/snp_selection/eqtl/UT_monocyte_eQTLProbesFDR0.05-ProbeLevel.tsv'
    eqtls_path = args.eqtl_path
    eqtl_dic = pd.read_csv(eqtls_path, sep='\t').set_index('SNPName')['ProbeName'].T.to_dict()
    coeqtls = pd.read_csv(coeqtl_path, sep='\t', index_col=0)
    coeqtls['eqtlgene'] = [eqtl_dic.get(snp) for snp in coeqtls['SNP']]
    coeqtls_lowest_nominalP = coeqtls.sort_values(by='MetaP', ascending=True).drop_duplicates(subset=['SNP'])
    coeqtls_lowest_nominalP_dict = coeqtls_lowest_nominalP.set_index('SNP')['MetaP'].T.to_dict()
    permutation_pvalues_df['SNP'] = [item.split('_')[0] for item in permutation_pvalues_df.index]
    permutation_pvalues_df['nominalP'] = [coeqtls_lowest_nominalP_dict.get(snp) for snp in
                                          permutation_pvalues_df['SNP']] # todo: check column name
    permutation_pvalues_df = permutation_pvalues_df.dropna(subset=['nominalP'])
    permutation_pvalues_df['pval_beta'] = [1-beta.sf(x[0], x[1], x[2]) for x in
                                               permutation_pvalues_df[['nominalP', 'beta_shape1', 'beta_shape2']].values]

    count_freq = lambda x:x[0][x[0]<x[1]].shape[0] / len(permutation_cols)
    permutation_pvalues_df['freq'] = [count_freq(x) for x in
                                               zip(permutation_pvalues_df[permutation_cols].values,
                                                   permutation_pvalues_df['nominalP'])]

    permutation_pvalues_df['scipy_beta_shape1'], permutation_pvalues_df['scipy_beta_shape2'] = \
    zip(*[beta.fit(x)[:2] for x in permutation_pvalues_df[permutation_cols].values])
    permutation_pvalues_df['scipy_pval_beta'] = [1-beta.sf(x[0], x[1], x[2]) for x in
                                               permutation_pvalues_df[['nominalP', 'scipy_beta_shape1', 'scipy_beta_shape2']].values]
    # over all eqtls, perform BH-FDR
    permutation_pvalues_df['qval'] = multipletests(permutation_pvalues_df['pval_beta'].values, method='fdr_bh')[1]
    # # determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
    # ub <- sort(fastqtl.df[fastqtl.df$qval > args$fdr, 'pval_beta'])[1]  # smallest p-value above FDR
    ub = permutation_pvalues_df[permutation_pvalues_df['qval']>=0.05].sort_values(by=['pval_beta'], ascending=True)['pval_beta'].values[0]
    # lb <- -sort(-fastqtl.df[fastqtl.df$qval <= args$fdr, 'pval_beta'])[1]  # largest p-value below FDR
    lb = permutation_pvalues_df[permutation_pvalues_df['qval']<=0.05].sort_values(by=['pval_beta'], ascending=False)['pval_beta'].values[0]
    # pthreshold <- (lb+ub)/2
    pthreshold = (ub + lb) / 2
    # cat("  * min p-value threshold @ FDR ", args$fdr, ": ", pthreshold, "\n", sep="")
    print('Minimum p-value threshold', pthreshold)
    # fastqtl.df[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold,
    #     fastqtl.df[, 'beta_shape1'], fastqtl.df[, 'beta_shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)
    permutation_pvalues_df['threshold_per_betadistribution'] = [beta.ppf(pthreshold, x[0], x[1]) for x in
                                                                permutation_pvalues_df[['beta_shape1',
                                                                                        'beta_shape2']].values]
    permutation_pvalue_threshold_dict = permutation_pvalues_df.set_index('SNP').T.to_dict()
    coeqtls['snp_beta_shape1'] = [permutation_pvalue_threshold_dict.get(snp)['beta_shape1'] for snp in coeqtls['SNP']]
    coeqtls['snp_beta_shape2'] = [permutation_pvalue_threshold_dict.get(snp)['beta_shape2'] for snp in coeqtls['SNP']]
    coeqtls['snp_pvalbeta'] = [permutation_pvalue_threshold_dict.get(snp)['pval_beta'] for snp in coeqtls['SNP']]
    coeqtls['snp_qval'] = [permutation_pvalue_threshold_dict.get(snp)['qval'] for snp in coeqtls['SNP']]
    coeqtls['gene2_pthreshold'] = [permutation_pvalue_threshold_dict.get(snp)['threshold_per_betadistribution']
                                   for snp in coeqtls['SNP']]
    issig = lambda x:True if x[0] <= x[1] else False
    coeqtls['gene2_isSig'] = [issig(item) for item in coeqtls[['MetaP', 'gene2_pthreshold']].values]
    significant_coeqtls = coeqtls[(coeqtls['snp_qval']<=0.05) & (coeqtls['gene2_isSig'])]
    # savepath = workdir/'output/UT_monocytes/concated_alltests_output.withP.tsv'
    saveprefix = args.saveprefix
    coeqtls.to_csv(f'{saveprefix}.all.tsv.gz', sep='\t', compression='gzip')
    significant_coeqtls.to_csv(f'{saveprefix}.sig.tsv.gz', sep='\t', compression='gzip')
    return coeqtls


if __name__ == '__main__':
    _ = main()
