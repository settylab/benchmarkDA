import anndata
import numpy as np
import pandas as pd

def csv2anndata(data_path: str, obs_path: str):
    # read the csv files
    X_pca = pd.read_csv(data_path, index_col=0)
    obs = pd.read_csv(obs_path, index_col=0)
    assert (X_pca.index == obs.index).all(), \
        "the index of data and meta data is different from each other"
    # create the anndata
    adata = anndata.AnnData(X_pca, obs=obs, dtype=np.float64)
    adata.obs.index.name = "cell"
    return adata
    

def calculate_outcome(bm: pd.DataFrame):
    # calculate the metrics, TP, FP, TN, FN and so on.
    TP = ((bm['true'] == bm['pred']) & (bm['pred'] != 'NotDA')).sum()
    FP = ((bm['true'] != bm['pred']) & (bm['pred'] != 'NotDA')).sum()
    FN = ((bm['true'] != bm['pred']) & (bm['pred'] == 'NotDA')).sum()
    TN = ((bm['true'] == bm['pred']) & (bm['pred'] == 'NotDA')).sum()
    
    metrics_dic = {
        'TP': TP, 'FP': FP,
        'FN': FN, 'TN': TN,
        'TPR': [TP / (TP + FN)],
        'FPR': [FP / (FP + TN)],
        'TNR': [TN / (TN + FP)],
        'FNR': [FN / (FN + TP)],
        'FDR': [FP / (TP + FP)],
        'Precision': [TP / (TP + FP)],
        'Power': [1 - FN / (FN + TP)],
        'Accuracy': [(TP + TN) / (TP + TN + FP + FN)]
    }
    
    return pd.DataFrame(metrics_dic, index=['metric'])


def out2bm(out, adata, runtime=None):
    # prepare the data frame as output
    bm = pd.DataFrame()
    bm['true_prob'] = adata.obs['Condition2_prob']
    bm['true'] = adata.obs['true_labels']
    bm['pred'] = out
    bm['method'] = 'meld'
    if runtime:
        bm['runtime'] = runtime
    
    return bm


def c_index(continuous_result, continuous_ground_truth):
    mask = continuous_ground_truth[None, :] < continuous_ground_truth[:, None]
    result = continuous_result[None, :] < continuous_result[:, None]
    cindex = np.sum(result[mask]) / np.sum(mask)
    return cindex
