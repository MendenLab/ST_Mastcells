from python_scripts.utils import rename_observables

import scanpy as sc
import numpy as np
import pandas as pd
from copy import copy

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

import os
from datetime import date


def plot_UMAP_counts(adata_tmp, disease, obs_names_tmp, obs_counts_tmp, save_folder, gene_symbol='TPSB2'):
    reds = copy(mpl.cm.Reds)
    reds.set_under("lightgray")
    if not isinstance(obs_names_tmp, list):
        gene_name = obs_names_tmp.split('_')[0]
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        sc.pl.umap(adata_tmp[adata_tmp.obs[obs_names_tmp] == 'Others'],
                   color=obs_names_tmp, frameon=True, ncols=1, palette=['lightgrey'],
                   ax=ax, show=False, legend_loc='right upper', size=10, title='', alpha=0.3)
        scpl = ax.scatter(
            adata_tmp[adata_tmp.obs[obs_names_tmp] == gene_name].obsm['X_umap'][:, 0],
            adata_tmp[adata_tmp.obs[obs_names_tmp] == gene_name].obsm['X_umap'][:, 1],
            cmap=reds, s=5, c=adata_tmp[adata_tmp.obs[obs_names_tmp] == gene_name].obs[obs_counts_tmp],
            edgecolor='k', lw=0.05)
        ax.set_ylabel('UMAP2', fontsize=18)
        ax.set_xlabel('UMAP1', fontsize=18)

        cb = fig.colorbar(scpl, orientation='vertical')
        cb.ax.set_ylabel('{} transcript level / spot'.format(gene_symbol), rotation=90, labelpad=15, fontsize=15)
        sns.despine(fig=fig, ax=ax, trim=False)
        plt.savefig(os.path.join(
            save_folder, 'UMAP_{}_counts_{}.pdf'.format(gene_symbol, disease)), bbox_inches='tight')
        plt.close()

        # cap colorbar to visually increase location of TPSB2
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        sc.pl.umap(adata_tmp[adata_tmp.obs[obs_names_tmp] == 'Others'],
                   color=obs_names_tmp, frameon=True, ncols=1, palette=['lightgrey'],
                   ax=ax, show=False, legend_loc='right upper', size=10, title='', alpha=0.3)
        scpl = ax.scatter(
            adata_tmp[adata_tmp.obs[obs_names_tmp] == gene_name].obsm['X_umap'][:, 0],
            adata_tmp[adata_tmp.obs[obs_names_tmp] == gene_name].obsm['X_umap'][:, 1],
            cmap=reds, s=5, c=adata_tmp[adata_tmp.obs[obs_names_tmp] == gene_name].obs[obs_counts_tmp],
            edgecolor='k', lw=0.05, vmax=5)
        ax.set_ylabel('UMAP2', fontsize=18)
        ax.set_xlabel('UMAP1', fontsize=18)

        cb = fig.colorbar(scpl, orientation='vertical')
        cb.ax.set_ylabel('{} transcript level / spot'.format(gene_symbol), rotation=90, labelpad=15, fontsize=15)
        sns.despine(fig=fig, ax=ax, trim=False)
        plt.savefig(os.path.join(
            save_folder, 'UMAP_{}_counts_{}_capped_5.pdf'.format(gene_symbol, disease)), bbox_inches='tight')
        plt.close()


def main(adata, save_folder, paired):

    if paired == 'paired':
        # Read out paired samples
        adata = adata[adata.obs['project'].isin(['P15509', 'P16357'])].copy()

    print('Number of samples: ', len(adata.obs['specimen'].unique()))
    print('Number of Lesional samples: ', len(adata.obs['sample'][adata.obs['biopsy_type'] == 'LESIONAL'].unique()))
    print('Number of Non lesional samples: ', len(
        adata.obs['sample'][adata.obs['biopsy_type'].str.contains('NON LESIONAL')].unique()))
    print('Number of patients: ', len(adata.obs['patient'].unique()))
    # Number of samples:  12
    # Number of Lesional samples:  6
    # Number of Non lesional samples:  6
    # Number of patients:  3

    # Get raw counts
    df_all = adata.to_df(layer='counts')

    # select spots which contain both TPSB2 (TPSAB1) and CMA1 (for all analyses)
    mast_cell_marker = ['TPSB2']

    # Visualisation on UMAP + H&E image (with all spots and only Mast cells)
    min_counts_threshold = 1
    mask_mastcells = np.all(df_all[mast_cell_marker] >= min_counts_threshold, axis=1)
    print('Number of Mast cells: ', np.sum(mask_mastcells))  # 1192 Mast cells for counts threshold >= 1

    obs_names = 'Mast cell_others'
    obs_counts = 'Mast cell_counts'

    mask_spots = mask_mastcells
    # add mast cell annotation to adata
    adata.obs[obs_names] = 'Others'
    adata.obs.loc[mask_spots, obs_names] = 'Mast cell'
    adata.obs[obs_counts] = 0
    adata.obs.loc[mask_spots, obs_counts] = np.sum(df_all.loc[mask_spots, mast_cell_marker], axis=1)

    # remove junction and merge dermis depth into dermis
    adata = rename_observables.remove_junction(adata=adata)

    # Plot counts of Mast cell genes
    plot_UMAP_counts(adata_tmp=adata, disease='Pso', obs_names_tmp=obs_names,
                     obs_counts_tmp=obs_counts, save_folder=save_folder,  gene_symbol='TPSB2')

    # read out adata plot infos
    df = pd.DataFrame.from_dict({
        'x': adata.obsm['X_umap'].T[0], 'y': adata.obsm['X_umap'].T[1], obs_counts: adata.obs[obs_counts],
        obs_names: adata.obs[obs_names], 'specimen': adata.obs['specimen']})
    df.to_excel(os.path.join(save_folder, 'Fig2iii_UMAP_biopsy_type.xlsx'), sheet_name='UMAP_biopsy_type')


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    do_paired = 'paired'
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/MendenLab_projects/psoriasis__Graham_Mackay/output",
        "Figure_2iii_TPSB2".format(do_paired), str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')
    # Keep Pso samples
    spatial_adata = spatial_adata[spatial_adata.obs['DISEASE'] == 'Pso'].copy()

    main(adata=spatial_adata, save_folder=savepath, paired=do_paired)
