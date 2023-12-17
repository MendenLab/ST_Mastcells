#!/usr/bin/env python
"""Main script for calculating the (Weighted) Spatial Pearson Correlation for different methods
    File name: main.py
    Author: Christina Hillig
    Credits: Christina Hillig, Ali Farnoud-Niedermayr, Michael Menden
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

from python_scripts.spatial_correlation import density_clustering
from python_scripts.utils import gene_lists

import os
from datetime import date
import scanpy as sc
import anndata
import pandas as pd
from operator import itemgetter


def main(save_folder: str, adata: anndata.AnnData, radius: [int, list], cond_genes: list, genes_resps: dict,
         find_associated_genes: bool = False, corr_method: str = 'pearson', get_plots: bool = False):
    """ Call conditional-based density clustering function

    Parameters
    ----------
    save_folder : str
        path to output directry
    adata : anndata.AnnData
        adata object containing either raw counts or both raw and normed counts.
        Later is necessary if find_associated_genes is True
    radius : int, list
        radius/distance from center spot to surrounding nearest neighbor spots
        e.g.  1 or list [1, 2, 3, ..]
    cond_genes : list
        name of genes to investigate in spots
    genes_resps : dict
        name of cond_genes associated genes such as responder
    find_associated_genes : bool
        if you want to find other than known associated genes of cond_genes
        if True, please provided filtered, normed adata object with sizefactors in .obs as input
    corr_method : str
        which weighted correlation method to use: spearman or pearson (default)
    get_plots : bool
        create evaluation plots

    Returns
    -------

    """
    # parameter
    tissue_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']
    epidermis_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']

    # # Run conditional density clustering
    # adata, counts_dict, df_stats_responders_in_vs_outlesion_sdc, \
    # df_stats_cytokines_responders_in_sdc, df_radius_vs_spearman = density_clustering.main(
    #     adata=adata, save_folder=save_folder, tissue_types=tissue_layers, epidermis_layers=epidermis_layers,
    #     radii=radius, get_plots=get_plots, corr_method=corr_method, conditional_genes=cond_genes,
    #     conditionalgenes_responders=genes_resps, find_associated_genes=find_associated_genes)
    # Run conditional density clustering
    adata, counts_dict, df_radius_vs_spearman = density_clustering.main(
        adata=adata, save_folder=save_folder, tissue_types=tissue_layers, epidermis_layers=epidermis_layers,
        radii=radius, get_plots=get_plots, corr_method=corr_method, conditional_genes=cond_genes,
        conditionalgenes_responders=genes_resps, find_associated_genes=find_associated_genes)

    # return adata, counts_dict, df_stats_responders_in_vs_outlesion_sdc, \
    #        df_stats_cytokines_responders_in_sdc, df_radius_vs_spearman
    return adata, counts_dict, df_radius_vs_spearman


if __name__ == '__main__':
    today = date.today()
    # replace os.environ['PYTHONPATH'].split(os.pathsep)[0] with sys.path[2] -> can be run then in terminal
    path = "/Volumes/CH__data/Projects/MendenLab_projects/psoriasis__Graham_Mackay"
    # create saving folder in current project path
    task = 'TPSB2_vs_CXCL17'
    savepath = os.path.join(path, "output", task, str(today))
    os.makedirs(savepath, exist_ok=True)

    # Load Raw anndata --> used for publication figure 4
    unpp_st_adata = sc.read(
        os.path.join('/Volumes/CH__data/Projects/data/annData_objects/spatial/2022-04-08',
                     "Spatial Transcriptomics_unpp_cleaned_LPADPso.h5"))
    # save .uns
    dict_spatial = unpp_st_adata.uns['spatial']
    # Subset to only lesion Pso samples
    unpp_st_adata = unpp_st_adata[unpp_st_adata.obs['DISEASE'] == 'Pso'].copy()
    # Read out paired samples
    unpp_st_adata = unpp_st_adata[unpp_st_adata.obs['project'].isin(['P15509', 'P16357'])].copy()
    # store spatial information
    unpp_st_adata.obs['sample'] = unpp_st_adata.obs['sample'].cat.remove_unused_categories()

    dict_spatial_wanted = dict((k, dict_spatial[k]) for k in list(unpp_st_adata.obs['sample'].cat.categories) if k in dict_spatial)
    unpp_st_adata.uns['spatial'] = dict_spatial_wanted

    # 1. Get cytokines and responders
    conditional_genes, conditionalgenes_responders = gene_lists.TPSB2_vs_CXCL17()

    radii = [0]

    # unpp_st_adata, final_counts_dict, final_df_stats_responders_in_vs_outlesion_sdc, \
    # final_df_stats_cytokines_responders_in_sdc, final_df_radius_vs_spearman = main(
    #     save_folder=savepath, adata=unpp_st_adata, cond_genes=conditional_genes,
    #     genes_resps=conditionalgenes_responders, radius=radii, get_plots=False, find_associated_genes=False,
    #     corr_method='spearman')

    unpp_st_adata, final_counts_dict, final_df_radius_vs_spearman = main(
        save_folder=savepath, adata=unpp_st_adata, cond_genes=conditional_genes,
        genes_resps=conditionalgenes_responders, radius=radii, get_plots=True, find_associated_genes=False,
        corr_method='spearman')

    writer = pd.ExcelWriter(os.path.join(savepath, 'Figure_Cocorrelation.xlsx'), engine='xlsxwriter')
    if not isinstance(final_df_radius_vs_spearman, list):
        max_corr_value = final_df_radius_vs_spearman[
            final_df_radius_vs_spearman.columns[
                final_df_radius_vs_spearman.columns.str.contains('correlation')]].idxmax(axis=0)
        # Read out counts and correlation values
        if len(conditional_genes) > 1:
            list_fig5eh_counts = list(itemgetter(*max_corr_value)(final_counts_dict))
            df_fig5eh_counts = pd.concat(list_fig5eh_counts, axis=0, keys=list(max_corr_value.index))
        else:
            list_fig5eh_counts = final_counts_dict[max_corr_value.values[0]]
            df_fig5eh_counts = pd.DataFrame(list_fig5eh_counts)
        df_fig5eh_correlation = final_df_radius_vs_spearman.iloc[max_corr_value]
        final_df_radius_vs_spearman.to_excel(writer, sheet_name='radius_vs_correlation')
        df_fig5eh_correlation.to_excel(writer, sheet_name='correlation')
    else:
        list_fig5eh_counts = final_counts_dict[list(final_counts_dict.keys())[0]]
        df_fig5eh_counts = pd.DataFrame(list_fig5eh_counts)

    df_fig5eh_counts.to_excel(writer, sheet_name='counts')
    writer.close()


