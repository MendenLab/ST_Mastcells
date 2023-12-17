from python_scripts.utils import rename_observables

import scanpy as sc
import numpy as np
import pandas as pd
from copy import copy
from operator import itemgetter

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages, FigureCanvasPdf
from statannot import add_stat_annotation

import os
from datetime import date


def plot_umap(adata_tmp, obs_counts, gene, pages, gene_symbol):
    # make red colormap
    reds = copy(mpl.cm.Reds)
    reds.set_under("lightgray")

    if gene != 'cell_types':
        gene_name = gene.split('_')[0]
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        sc.pl.umap(adata_tmp[adata_tmp.obs[gene] == 'Others'], color=gene, frameon=True, ncols=1,
                   palette=['lightgrey'],
                   ax=ax[0], show=False, title='', s=5)
        sc.pl.umap(adata_tmp[adata_tmp.obs[gene] == gene_name], color=gene, frameon=True, ncols=1,
                   palette=['orangered'],
                   ax=ax[0], show=False, title='', s=20)
        ax[0].set_ylabel('UMAP2', fontsize=18)
        ax[0].set_xlabel('UMAP1', fontsize=18)

        sc.pl.umap(adata_tmp[adata_tmp.obs[gene] == 'Others'],
                   color=gene, frameon=True, ncols=1, palette=['lightgrey'],
                   ax=ax[1], show=False, legend_loc='right upper', size=5, title='')
        scpl = ax[1].scatter(
            adata_tmp[adata_tmp.obs[gene] == gene_name].obsm['X_umap'][:, 0],
            adata_tmp[adata_tmp.obs[gene] == gene_name].obsm['X_umap'][:, 1],
            cmap=reds, s=2, c=adata_tmp[adata_tmp.obs[gene] == gene_name].obs[obs_counts])
        ax[1].set_ylabel('UMAP2', fontsize=18)
        ax[1].set_xlabel('UMAP1', fontsize=18)

        cb = fig.colorbar(scpl, orientation='vertical')
        cb.ax.set_ylabel('{} transcript level / spot'.format(gene_symbol), rotation=90, labelpad=15, fontsize=15)
    else:
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        sc.pl.umap(adata_tmp, color=gene, frameon=True, ncols=1,
                   ax=ax, show=False, size=5)

        ax.set_ylabel('UMAP2', fontsize=18)
        ax.set_xlabel('UMAP1', fontsize=18)
    sns.despine(fig=fig, ax=ax, trim=False)
    plt.tight_layout()

    # fig.savefig(os.path.join(save_folder, '{}_{}_UMAP.pdf'.format(disease, gene.split('_')[0])),
    #             bbox_inches='tight')

    canvas = FigureCanvasPdf(fig)
    canvas.print_figure(pages)

    plt.close(fig=fig)


def plot_he(adata, gene, df, ind, obs_counts, pages, gene_symbol):
    for specimen in adata.obs['specimen'].cat.categories:
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
            # ax.imshow(sample_adata.uns['spatial'][sample]['images']['hires'])
            # ax.invert_yaxis()
            # scpl = ax.scatter(
            #     sample_adata.obsm[
            #         'spatial'][:, 0] * sample_adata.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef'],
            #     sample_adata.obsm[
            #         'spatial'][:, 1] * sample_adata.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef'],
            #     cmap=cmap.cmap, s=2, norm=norm, c=sample_adata.obs[obs_counts])

            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title("{}\n{}".format(biopsy_type, gene_name))
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_visible(False)

            # cb = fig.colorbar(scpl, orientation='vertical', ticks=ticks)
            # cb.ax.set_ylabel('{} transcript level / spot'.format(gene_symbol), rotation=90)
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


def plots(adata_tmp, disease, obs_names_tmp, obs_counts_tmp, save_folder, gene_symbol='TPSB2'):
    with PdfPages(os.path.join(save_folder, 'UMAP_{}_{}.pdf'.format(gene_symbol, disease))) as pages:
        if isinstance(obs_names_tmp, list):
            for ind, gene in enumerate(obs_names_tmp):
                plot_umap(adata_tmp, obs_counts=obs_counts_tmp[ind], gene=gene, pages=pages, gene_symbol=gene_symbol)
        else:
            plot_umap(adata_tmp, obs_counts=obs_counts_tmp, gene=obs_names_tmp, pages=pages, gene_symbol=gene_symbol)


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

    # remove junction and merge dermis depth into dermis
    adata = rename_observables.remove_junction(adata=adata)

    # Select from first patient the second lesion replicate
    specimen = 'V19523-003-V4_1'
    adata_specimen = adata[adata.obs['specimen'] == specimen].copy()
    sample = adata_specimen.obs['sample'].unique()[0]
    biopsy_type = adata_specimen.obs['biopsy_type'].unique()[0]
    patient = adata_specimen.obs['patient'].unique()[0]

    # Mast cells biopsy type
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 10))
    fig.suptitle("patient: {}\n{}".format(patient, biopsy_type))
    sc.pl.spatial(adata=adata_specimen, color='tissue_layer', library_id=sample, ax=ax[0], title='',
                  show=False)
    ax[0].set_xlabel('')
    ax[0].set_ylabel('')
    ax[0].invert_yaxis()
    ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), title='', fancybox=True, shadow=True,
                 ncol=2, prop={'size': 12}, frameon=False, title_fontsize=14)

    sc.pl.spatial(adata=adata_specimen, color=obs_names, library_id=sample, ax=ax[1], title='',
                  palette=['orangered', 'grey'], show=False)
    ax[1].set_xlabel('')
    ax[1].set_ylabel('')
    ax[1].invert_yaxis()
    ax[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), title='', fancybox=True, shadow=True,
                 ncol=1, prop={'size': 12}, frameon=False, title_fontsize=14)

    fig.savefig(os.path.join(save_folder, 'H&E_Mastcells_{}.pdf'.format('tissuelayers_biopsy_type')),
                bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    # create saving folder in current project path
    today = date.today()
    do_paired = 'paired'
    savepath = os.path.join(
        "/Volumes/CH__data/Projects/MendenLab_projects/psoriasis__Graham_Mackay/output",
        "SFig_3A_TPSB2".format(do_paired), str(today))
    os.makedirs(savepath, exist_ok=True)

    spatial_adata = sc.read(
        '/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08/st_QC_normed_BC_project_PsoADLP.h5')
    # Keep Pso samples
    spatial_adata = spatial_adata[spatial_adata.obs['DISEASE'] == 'Pso'].copy()

    main(adata=spatial_adata, save_folder=savepath, paired=do_paired)
