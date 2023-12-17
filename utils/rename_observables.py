import pandas as pd
import numpy as np
from operator import itemgetter


def interface_to_epidermis(adata, tissue_layers_obsname, skin_layers):
    # Rename tissue region 'INTERFACE' to upper, middle or basal EPIDERMIS because some spots got both labels
    m_interface = adata.obs['JUNCTION'] == 1
    if isinstance(skin_layers, list):
        df_temp = adata.obs[skin_layers][m_interface]
    else:
        df_temp = adata.obs[[skin_layers]][m_interface]
    df_temp = df_temp.loc[:, df_temp.columns].replace(1, pd.Series(df_temp.columns, df_temp.columns))
    df_temp['JUNCTION'] = '0'
    for col in df_temp.columns[:-1]:
        df_temp['JUNCTION'] += df_temp[col].astype(str)
    df_temp["JUNCTION"] = df_temp.JUNCTION.str.replace('0', '')
    adata.obs[tissue_layers_obsname][m_interface] = df_temp["JUNCTION"]
    adata.obs[tissue_layers_obsname] = adata.obs[tissue_layers_obsname].cat.remove_unused_categories()

    return adata


def remove_junction(adata):
    # remove JUNCTION
    mask = adata.obs["spot_type"] == "JUNCTION"
    mask_middle_epidermis = adata.obs["middle EPIDERMIS"] == 1
    mask_basal_epidermis = adata.obs["basal EPIDERMIS"] == 1
    mask_upper_epidermis = adata.obs["upper EPIDERMIS"] == 1
    mask_dermis = adata.obs["DERMIS"] == 1
    adata.obs.loc[mask & mask_middle_epidermis, "spot_type"] = "middle EPIDERMIS"
    adata.obs.loc[mask & mask_basal_epidermis, "spot_type"] = "basal EPIDERMIS"
    adata.obs.loc[mask & mask_upper_epidermis, "spot_type"] = "upper EPIDERMIS"
    adata.obs.loc[mask & mask_dermis, "spot_type"] = "DERMIS"
    adata.obs["spot_type"] = adata.obs["spot_type"].cat.remove_unused_categories()

    adata.obs["spot_type"] = adata.obs['spot_type'].cat.reorder_categories(
        ["upper EPIDERMIS", "middle EPIDERMIS", "basal EPIDERMIS", "DERMIS", 'MUSCLE', "VESSEL", "HAIR FOLLICLE",
         "SEBACEOUS GLAND", "SWEAT GLAND"])

    colors_spottype = dict(zip(adata.obs['spot_type'].cat.categories.to_list(), [
        "#1f77b4", "#ff7f0e", "#279e68", '#e377c2', '#8c564b', '#aa40fc', '#b5bd61', '#17becf', '#aec7e8']))

    adata.uns['spot_type_colors'] = itemgetter(*adata.obs.loc[
                   mask, 'spot_type'].cat.categories.to_list())(colors_spottype)

    # remove JUNCTION in tissue_layer
    adata.obs["tissue_layer"] = adata.obs["tissue_layer"].astype(str)
    mask = adata.obs["tissue_layer"] == "JUNCTION"
    mask_middle_epidermis = adata.obs["middle EPIDERMIS"] == 1
    mask_basal_epidermis = adata.obs["basal EPIDERMIS"] == 1
    mask_upper_epidermis = adata.obs["upper EPIDERMIS"] == 1
    mask_dermis = adata.obs["DERMIS"] == 1
    mask_dermis_depth1 = adata.obs["DERdepth1"] == 1
    mask_dermis_depth2 = adata.obs["DERdepth2"] == 1
    adata.obs.loc[mask & mask_middle_epidermis, "tissue_layer"] = "middle EPIDERMIS"
    adata.obs.loc[mask & mask_basal_epidermis, "tissue_layer"] = "basal EPIDERMIS"
    adata.obs.loc[mask & mask_upper_epidermis, "tissue_layer"] = "upper EPIDERMIS"
    adata.obs.loc[mask & mask_dermis_depth1, "tissue_layer"] = "DERdepth1"
    adata.obs.loc[mask & mask_dermis_depth2, "tissue_layer"] = "DERdepth2"
    adata.obs.loc[mask & mask_dermis, "tissue_layer"] = "DERdepth1"

    # merge dermis depths into dermis
    adata.obs.loc[mask_dermis, "tissue_layer"] = 'DERMIS'
    adata.obs["tissue_layer"] = adata.obs["tissue_layer"].astype('category')
    adata.obs["tissue_layer"] = adata.obs["tissue_layer"].cat.reorder_categories(
        ["upper EPIDERMIS", "middle EPIDERMIS", "basal EPIDERMIS", 'DERMIS'])

    adata.uns['tissue_layer_colors'] = np.asarray(['limegreen', 'mediumseagreen', 'darkgreen', '#e377c2'])

    adata.uns['biopsy_type_colors'] = np.asarray(['lightcoral', 'teal'])

    return adata
