from python_scripts.utils import rename_observables

import scanpy as sc
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

import os
from datetime import date


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

    # select spots which contain TPSB2 (for all analyses)
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

    fig, ax = plt.subplots(figsize=(5, 5))
    sc.pl.umap(adata=adata, color='biopsy_type', ax=ax, show=False, title='', size=15, alpha=0.3)
    sc.pl.umap(adata[adata.obs[obs_names] == 'Mast cell'], color='biopsy_type', frameon=True, ncols=1,
               ax=ax, show=False, legend_loc='right upper', size=25, title='')
    sns.despine(fig=fig, ax=ax, top=True, right=True, left=False, bottom=False)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), title='',
              fancybox=True, shadow=True, ncol=2, prop={'size': 12}, frameon=False, title_fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(
        save_folder, 'UMAP_biopsy_type_Pso.pdf'), bbox_inches='tight')
    plt.close()

    # read out data
    df = pd.DataFrame.from_dict({
        'x': adata.obsm['X_umap'].T[0], 'y': adata.obsm['X_umap'].T[1], 'biospy type': adata.obs['biopsy_type'],
        obs_counts: adata.obs[obs_counts], obs_names: adata.obs[obs_names], 'specimen': adata.obs['specimen']})
    # map color
    df['biopsy_type_colors'] = adata.uns['biopsy_type_colors']

    df.to_excel(os.path.join(save_folder, 'Fig2i_UMAP_biopsy_type.xlsx'), sheet_name='UMAP_biopsy_type')


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    do_paired = 'paired'
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/MendenLab_projects/psoriasis__Graham_Mackay/output",
        "Figure_2Bi_TPSB2".format(do_paired), str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')
    # Keep Pso samples
    spatial_adata = spatial_adata[spatial_adata.obs['DISEASE'] == 'Pso'].copy()

    main(adata=spatial_adata, save_folder=savepath, paired=do_paired)
