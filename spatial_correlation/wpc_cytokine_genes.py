from python_scripts.spatial_correlation import corr_statistics as corr_stats

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from scipy.stats.stats import pearsonr
import scipy.stats.stats as stats


def plot_wpc_cyto_genes(df_cytos_wpc, genes, cyto, save_folder):
    df_wpc = df_cytos_wpc[cyto][genes].transpose()
    mask_poscorr = (df_wpc['p-value'] < 0.05) & (df_wpc['WPC'] > 0.7)
    mask_negcorr = (df_wpc['p-value'] < 0.05) & (df_wpc['WPC'] < -0.7)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.grid(False)
    ax.scatter(df_wpc['WPC'], -np.log10(df_wpc['p-value']), s=2)
    ax.scatter(df_wpc[mask_poscorr]['WPC'], -np.log10(df_wpc[mask_poscorr]['p-value']),
               s=2, c='blue')
    ax.scatter(df_wpc[mask_negcorr]['WPC'], -np.log10(df_wpc[mask_negcorr]['p-value']),
               s=2, c='orange')
    # label
    ax.set_xlabel('Weighted Pearson Correlation')
    ax.set_ylabel('p-value')
    ax.hlines(y=-np.log10(0.05), xmin=-1, xmax=1)
    # remove borders from plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    fig.savefig(os.path.join(save_folder, 'WPC__{}_correlated_genes.pdf'.format(cyto)))
    plt.close()

    df_wpc[mask_poscorr & mask_negcorr].to_excel(os.path.join(save_folder, '{}_WPC_genes.xlsx'.format(cyto)))


def get_correlation_cyto_genes(df_counts, genes, save_folder):
    # remove zero count genes
    df_counts = df_counts.loc[:, (df_counts.sum(axis=0) != 0).values]
    genes = np.intersect1d(genes, df_counts.columns)
    # Taking into account the number of cytokines in a cluster
    weighted_correlations = {key: [] for key in ['IL17A', 'IL13', 'IFNG']}
    for cyto in ['IL17A', 'IL13', 'IFNG']:
        weighted_cytoname = "_".join(['weighted', cyto])
        # Apply weights to cytokine counts and add it to df: Multiply with cluster size
        weighted_cytocounts = df_counts[cyto] * df_counts['Cluster_size']
        df_counts[weighted_cytoname] = weighted_cytocounts.values.astype(np.float64)
        df_corr_temp = pd.DataFrame(columns=genes, index=['WPC', 'p-value'])
        for col in genes:
            # weighted spatial correlation
            df_corr_temp.loc['WPC', col] = corr_stats.weighted_corr(
                x=df_counts[cyto], y=df_counts[col], w=df_counts['Cluster_size'])
            # Calculate p-value
            df_corr_temp.loc['p-value', col] = corr_stats.apply_statstest(
                df_highcounts=df_counts[col], df_lowcounts=weighted_cytocounts,
                correlation_value=df_corr_temp.loc['WPC', col])

        weighted_correlations[cyto] = df_corr_temp

    weighted_correlations_new = pd.concat(weighted_correlations, axis=1)
    df_cyto_wc = weighted_correlations_new.loc[['WPC', 'p-value'], :]

    del weighted_correlations_new

    df_cyto_wc.to_csv(os.path.join(save_folder, 'Cytokine_weighted_correlation.csv'))

    for cyto in ['IL17A', 'IL13', 'IFNG']:
        plot_wpc_cyto_genes(df_cytos_wpc=df_cyto_wc, genes=genes, cyto=cyto, save_folder=save_folder)


def get_semibulk_correlation_cyto_genes(adata):
    correlations = {}
    spearman_correlations = {}
    columns = adata.to_df().columns.tolist()
    columns.remove('IL17A')
    df_cyto = adata[(adata.obs['IL17A_in_sdcc'] == 1) | (adata.obs['IL17A_in_sdcc'] == 2)].to_df()['IL17A']
    # Set 0 to NaN
    df_cyto[df_cyto == 0] = np.nan
    df_others = adata[(adata.obs['IL17A_in_sdcc'] == 1) | (adata.obs['IL17A_in_sdcc'] == 2)].to_df()
    # Gives result for semi-bulk approach
    for col in columns:
        # Pearson correlation
        correlations[col] = pearsonr(df_cyto.values, df_others[col].values)
        # Spearman correlation
        spearman_correlations[col] = stats.spearmanr(df_cyto.values, df_others[col].values)

    result = pd.DataFrame.from_dict(correlations, orient='index')
    result.columns = ['PCC', 'p-value']

    spearman_result = pd.DataFrame.from_dict(spearman_correlations, orient='index')
    spearman_result.columns = ['PCC', 'p-value']

    fig, ax = plt.subplots()
    ax.scatter(result['PCC'], -np.log10(result['p-value']), s=1)
    ax.hlines(y=-np.log10(0.05), xmin=-0.2, xmax=0.5)

    print(result.sort_index())

