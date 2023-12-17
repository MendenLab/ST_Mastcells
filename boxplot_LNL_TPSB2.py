import scanpy as sc
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

import os
from datetime import date


def plot_boxplot(df_melted, y, title, save_folder, log=False):
    fig, ax = plt.subplots(figsize=(4, 4))
    sns.boxplot(data=df_melted, x='biopsy_type', y=y, ax=ax, palette=['lightcoral', 'teal'], showfliers=False)
    sns.stripplot(data=df_melted, x='biopsy_type', y=y, ax=ax, color='k')
    if len(df_melted['biopsy_type'].unique()) > 1:
        add_stat_annotation(ax, data=df_melted, x='biopsy_type', y=y,
                            box_pairs=[("NON LESIONAL", "LESIONAL")],
                            test='Mann-Whitney-gt', text_format='star', loc='outside', verbose=2)
    sns.despine(fig=fig, ax=ax, trim=False)
    if log:
        ax.set_yscale('log')
    ax.set_ylim([1, df_melted[y].max() + 1])
    ax.set_ylabel(title, fontsize=14)
    ax.xaxis.set_tick_params(labelsize=12)
    ax.set_xlabel('')
    plt.tight_layout()
    plt.savefig(os.path.join(
        save_folder, 'Boxplot_biopsy_type_Pso_{}.pdf'.format(title)), bbox_inches='tight')
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

    # CXCL17
    mask_cxcl17 = df_all['CXCL17'] >= min_counts_threshold
    print('Number of all CXCL17+ spots: ', np.sum(mask_cxcl17))  # 81
    print('Number of Mast cell CXCL17+ spots: ', np.sum(mask_spots & mask_cxcl17))  # 38
    obs_names_cxcl17 = 'CXCL17_others'
    obs_counts_cxcl17 = 'CXCL17_counts'

    # add mast cell annotation to adata
    adata.obs[obs_names_cxcl17] = 'Others'
    adata.obs.loc[mask_spots & mask_cxcl17, obs_names_cxcl17] = 'CXCL17'
    adata.obs[obs_counts_cxcl17] = 0
    adata.obs.loc[mask_spots & mask_cxcl17, obs_counts_cxcl17] = df_all.loc[mask_spots & mask_cxcl17, 'CXCL17']

    df_CXCL17 = pd.DataFrame(data=adata[mask_spots].to_df(layer='counts')['CXCL17'], columns=['CXCL17'])
    df_boxplot_CXCL17 = pd.concat(
        [df_CXCL17, adata.obs.loc[mask_spots, [obs_counts, obs_names, 'biopsy_type', 'specimen']]],
        axis=1)
    df_boxplot_CXCL17 = df_boxplot_CXCL17.loc[df_boxplot_CXCL17['CXCL17'] > 0, :]
    df_boxplot_CXCL17['specimen'] = df_boxplot_CXCL17['specimen'].cat.remove_unused_categories()

    plot_boxplot(df_melted=df_boxplot_CXCL17, y='CXCL17', title='CXCL17 transcript level',
                 save_folder=save_folder)

    df_boxplot = adata.obs.loc[mask_spots, [obs_counts, obs_names, 'biopsy_type', 'specimen']]
    df_boxplot['specimen'] = df_boxplot['specimen'].cat.remove_unused_categories()

    # p-value annotation legend:
    # ns: 5.00e-02 < p <= 1.00e+00
    # *: 1.00e-02 < p <= 5.00e-02
    # **: 1.00e-03 < p <= 1.00e-02
    # ***: 1.00e-04 < p <= 1.00e-03
    # ****: p <= 1.00e-04
    #
    # LESIONAL v.s. NON LESIONAL: Mann-Whitney-Wilcoxon test greater with Bonferroni correction,
    # P_val=6.421e-06 U_stat=1.645e+05
    plot_boxplot(df_melted=df_boxplot, y='Mast cell_counts', title='TPSB2 transcript level',
                 save_folder=save_folder)

    # Save infos for paper
    dict_total_num_spots_specimen = dict(
        zip(adata.obs.groupby('specimen').size().index.to_list(),
            adata.obs.groupby('specimen').size().values))
    df_counts = adata.obs[[obs_counts, 'specimen']].groupby('specimen').agg('sum')
    mapping_type_specimen = dict(zip(adata.obs['specimen'], adata.obs['biopsy_type']))
    mapping_patient_specimen = dict(zip(adata.obs['specimen'], adata.obs['patient']))
    df_counts['Mast cell spots'] = adata.obs.loc[
        adata.obs[obs_names] == 'Mast cell', [obs_names, 'specimen']].groupby('specimen').value_counts().values
    df_counts['biopsy_type'] = list(map(mapping_type_specimen.get, list(df_counts.index)))
    df_counts['total number of spots'] = list(map(dict_total_num_spots_specimen.get, list(df_counts.index)))
    df_counts['patient'] = list(map(mapping_patient_specimen.get, list(df_counts.index)))
    df_counts.to_excel(os.path.join(save_folder, 'Manuscript_Mastcell_information.xlsx'))


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    do_paired = 'paired'
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/MendenLab_projects/psoriasis__Graham_Mackay/output",
        "Figure_2v_TPSB2".format(do_paired), str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')
    # Keep Pso samples
    spatial_adata = spatial_adata[spatial_adata.obs['DISEASE'] == 'Pso'].copy()

    main(adata=spatial_adata, save_folder=savepath, paired=do_paired)
