from python_scripts.utils import rename_observables

import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages, FigureCanvasPdf

import os
from datetime import date


def plot_he(adata, gene, df, ind, obs_counts, pages, gene_symbol):
    specimen = 'V19523-003-V4_1'
    gene_name = gene.split('_')[0]
    sample_adata = adata[(adata.obs['specimen'] == specimen) & (adata.obs[gene] == gene_name)].copy()
    specimen_adata = adata[(adata.obs['specimen'] == specimen)].copy()

    if sample_adata.shape[0] != 0:
        print(specimen)
        sample = sample_adata.obs['sample'].cat.categories[0]
        biopsy_type = sample_adata.obs['biopsy_type'].cat.categories[0]
        disease = sample_adata.obs['DISEASE'].cat.categories[0]

        # Save number of gene+ spots and also add total number of counts per samples
        df.loc[ind, :] = [sample, disease, gene, sample_adata.obs[obs_counts].sum(),
                          sample_adata.shape[0]]
        ind = ind + 1

        # Option 3:
        max_val = sample_adata.obs[obs_counts].max()
        norm = mpl.colors.Normalize(vmin=sample_adata.obs[obs_counts].min(), vmax=max_val)
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.viridis)

        if max_val < 10:
            ticks = np.arange(1, max_val + 1, 1)
        elif max_val < 50:
            ticks = np.linspace(1, max_val, 4, endpoint=True)
            ticks = ticks.round()
        else:
            ticks = np.linspace(1, max_val, 8, endpoint=True)
            ticks = ticks.round()

        fig, ax = plt.subplots()
        # plot whole specimen window
        sc.pl.spatial(adata=specimen_adata,
                      color='biopsy_type', library_id=sample, ax=ax, title='', size=0,
                      palette=['darkorange', 'grey'], show=False, legend_loc='upper right')
        sc.pl.spatial(adata=sample_adata, color=obs_counts, library_id=sample, ax=ax, title='',
                      cmap=cmap.cmap, show=False)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title("{}\n{}".format(biopsy_type, gene_name))
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)

        cbar_ax = fig.axes[-1]
        cbar_ax.get_yaxis().labelpad = 15
        cbar_ax.get_yaxis().spacing = 'proportional'
        cbar_ax.get_yaxis().format = '%1i'
        cbar_ax.set_yticks(ticks)
        cbar_ax.set_ylabel('{} transcript level / spot'.format(gene_symbol), rotation=90, fontsize=12)
        cbar_ax.tick_params(labelsize=10)

        canvas = FigureCanvasPdf(fig)
        canvas.print_figure(pages)
        plt.close(fig=fig)

    return df, ind


def plot_celltype_markers_per_specimen(adata, obs_names_tmp, obs_counts_tmp, save_folder, gene_symbol='TPSB2'):
    # remove unused categories
    adata.obs['specimen'] = adata.obs['specimen'].cat.remove_unused_categories()

    df = pd.DataFrame(columns=['Sample', 'Disease', 'Gene', 'Total No. Counts', 'No. Spots'])

    # Plot H&E images
    with PdfPages(os.path.join(save_folder, 'H&E_{}.pdf'.format(gene_symbol))) as pages:
        ind = 0
        if isinstance(obs_names_tmp, list):
            for gene in obs_names_tmp:
                df_tmp, ind = plot_he(adata, gene=gene, df=df, ind=ind, obs_counts=obs_counts_tmp[ind],
                                      pages=pages, gene_symbol=gene_symbol)
                df = pd.concat([df, df_tmp], axis=1)
        else:
            df, _ = plot_he(adata, gene=obs_names_tmp, df=df, ind=ind, obs_counts=obs_counts_tmp,
                            pages=pages, gene_symbol=gene_symbol)

    df.to_excel(os.path.join(save_folder, 'Information_{}.xlsx'.format(gene_symbol)))


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

    # remove junction and merge dermis depth into dermis
    adata = rename_observables.remove_junction(adata=adata)

    plot_celltype_markers_per_specimen(
        adata=adata, obs_names_tmp=obs_names_cxcl17, obs_counts_tmp=obs_counts_cxcl17,
        save_folder=save_folder, gene_symbol='CXCL17')

    plot_celltype_markers_per_specimen(
        adata=adata, obs_names_tmp=obs_names, obs_counts_tmp=obs_counts,
        save_folder=save_folder, gene_symbol='TPSB2')


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    do_paired = 'paired'
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/MendenLab_projects/psoriasis__Graham_Mackay/output",
        "SFig_2A_counts_TPSB2_CXCL17".format(do_paired), str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')
    # Keep Pso samples
    spatial_adata = spatial_adata[spatial_adata.obs['DISEASE'] == 'Pso'].copy()

    main(adata=spatial_adata, save_folder=savepath, paired=do_paired)
