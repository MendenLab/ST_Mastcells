from python_scripts.spatial_correlation import tools, helper_functions as ht


def data_preparation(adata, tissue_layers, epidermis_layers, conditional_genes, conditionalgenes_responders):
    """Prepare data before applying conditional density clustering algorithm

    Parameters
    ----------
    adata : annData
    tissue_layers : list
    conditional_genes : list
    epidermis_layers : str, list
    conditionalgenes_responders : dict

    Returns
    -------

    """
    # save .uns
    dict_spatial = adata.uns['spatial']

    # Subset adata to tissue_types of interest: upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS'
    if tissue_layers:
        bool_col = adata.obs[tissue_layers] == 1
        merged = bool_col.sum(axis=1)
        adata = adata[merged == 1].copy()
        # Rename tissue region 'JUNCTION' to basal EPIDERMIS because some spots got both labels
        adata = ht.interface_to_epidermis(adata, tissue_layers='tissue_layer', epidermis_layers=epidermis_layers)

    # Get counts of cyotkines and their responders in the EPIDERMIS
    # - distance between spots: 100 µm, spot diameter: 55 µm
    # - better: use index array
    resp_label = []
    for cyto in conditional_genes:
        adata = tools.add_columns_genes(adata=adata, genes=cyto, label=cyto, count_threshold=1)

    for ind, cyto_resps in enumerate(conditionalgenes_responders.keys()):
        resp_label.append("_".join([cyto_resps, "Responders"]))
        adata = tools.add_columns_genes(adata=adata, genes=conditionalgenes_responders[cyto_resps],
                                        label=resp_label[ind], count_threshold=1)

        # store spatial information
        adata.obs['sample'] = adata.obs['sample'].cat.remove_unused_categories()

        dict_spatial_wanted = dict(
            (k, dict_spatial[k]) for k in list(adata.obs['sample'].cat.categories) if k in dict_spatial)
        adata.uns['spatial'] = dict_spatial_wanted

    return adata
