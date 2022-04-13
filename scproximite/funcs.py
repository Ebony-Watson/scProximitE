###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

"""
Utility functions for evalualtion framework for proximity metric performance in scRNA-seq

Author: Ebony Watson
Date: 13-April-2022
"""

import genieclust  # Only required for PSI
import matplotlib.pyplot as plt
from natsort import natsorted  # Only required for Agglomerative clustering
import numpy as np
import os
import otscomics  # Only required for Optimal Transport
import pandas as pd
import pickle
import scanpy as sc
import scipy as sp
from scipy.stats import kruskal, rankdata, tiecorrect  # Only required for Kruskal test
import scikit_posthocs as sph  # Only required for Dunn test
import seaborn as sns
from sklearn import metrics as ms
from sklearn.cluster import AgglomerativeClustering
import torch  # Only required for Optimal Transport
from tqdm import tqdm  # Only required for Optimal Transport

sns.set_theme(style="white")  # Only required for Dunn vis


def proximity_matrix(data, metric, anndata, save_rep="prox_matrix"):
    """
    
    Builds square similarity (distance) matrix and adds it to anndata object

    Parameters
    ----------

     data : 2D np.ndarray
         cell x gene Data  (i.e. `anndata.layers["normalised"] or 'anndata.X')`)
     metric : str
         proximity metric with which to calculate the matrix
     anndata : anndata.AnnData
         AnnData object to add resulting matrix to (in `.uns` slot)
     save_rep : str, optional (default="prox_matrix")
         name of `.obsp` key to save similarity matrix to within anndata
     
     Returns
     -------
     
     `anndata` is edited in place, adding similarity matrix to `anndata.obsp[save_rep]`
     
     """

    data = sp.sparse.csr_matrix(data)

    if metric == 'spearman':
        prox_mat = sp.stats.spearmanr(data.todense(), axis=1)
        anndata.obsp[save_rep] = 1 - prox_mat.correlation

    elif metric == 'OT':
        data = data.T
        data = pd.DataFrame(data.todense())
        data = data.iloc[np.argsort(data.std(1))[::-1][:10_000]]  # Select 10k most varying features
        data = np.array(data)
        data_norm = data / data.sum(0)
        data_norm += 1e-6
        data_norm /= data_norm.sum(0)
        C = otscomics.cost_matrix(data, 'cosine')  # Compute cost matrix with cosine, removed: n_iter=5,
        D_ot = otscomics.OT_distance_matrix(data_norm, C, eps=.1, n_batches=50, dtype=torch.float32)
        D_ot[D_ot < 0] = 0  # Enforce positivity, just in case
        anndata.obsp[save_rep] = D_ot

    elif metric in ['dice', 'hamming', 'jaccard', 'kulsinski', 'rogerstanimoto', 'sokalmichener', 'sokalsneath',
                    'yule']:
        anndata.layers['bool_X'] = data.todense().astype(bool)
        anndata.obsp[save_rep] = ms.pairwise_distances(anndata.layers['bool_X'], metric=metric)

    elif metric == 'pearson':
        anndata.obsp[save_rep] = ms.pairwise_distances(data.todense(), metric='correlation')

    else:
        anndata.obsp[save_rep] = ms.pairwise_distances(data.todense(), metric=metric)

    return anndata


def merge_metrics(mat, anndata, metric):
    """
    Performs

    Parameters
    ----------

    Returns
    -------

    """
    mat = np.nan_to_num(mat, copy=False)  # replace any nan/inf vals with 0

    if metric == 'phi_s':
        anndata.obsp[f'{metric}_matrix'] = abs(mat)  # take absolute values of matrix & store in anndata object

    elif metric == 'zi_kendall':
        mat = (mat - np.min(mat)) / (np.max(mat) - np.min(mat))  # convert between 0 and 1
        np.fill_diagonal(mat, 0)
        anndata.obsp[f'{metric}_matrix'] = 1 - mat  # subtract values from 1 & store in anndata object

    else:
        anndata.obsp[f'{metric}_matrix'] = 1 - mat  # subtract values & store in anndata object
    print(f'Outcome for evaluation for {metric} as a true distance matrix is',
          sp.spatial.distance.is_valid_dm(anndata.obsp[f'{metric}_matrix'], tol=0.000000000001))

    return anndata


def clustering(anndata, graph=None, prox_matrix=None, save_rep='clusters', method='leiden', optimisation_labs=None,
               linkage='average', rs=1):
    """
            Performs clustering of k-nearest neighbor graph (sparse connectivity) and adds to anndata object

            Parameters
            ----------
            graph: knn graph
            anndata : anndata.AnnData
                AnnData object to add resulting graph to (in `.uns` slot)
                Sparse connectivity graph (i.e. `anndata.obsp["euclidean_KNNG"])`)
            save_rep : str, optional (default="clusters")
                name of `.obs` key to save cluster IDs to within anndata
            method : str, optional (default='leiden')
            optimisation_labs : pandas series, optional (default = "None")
                Ground truth cell-type labels (i.e. anndata.obs['cell_line'])
                If provided, resolution parameter of the clustering method is tuned until number of clusters derived
                matches ground truth number.
            rs: int, optional (default=1)
            Returns
            -------
            `anndata` is edited in place, adding cluster IDs to `anndata.obs[save_rep]`
            """
    if method == 'leiden':
        c = 0
        if optimisation_labs is not None:
            maxn = 2
            minn = 0
            x = 1
            prev_x = x
            sc.tl.leiden(adata=anndata, adjacency=graph, resolution=x, key_added=save_rep, random_state=rs,
                         directed=False, use_weights=False)
            N = len(set(optimisation_labs))
            while len(set(anndata.obs[save_rep])) != N and c < 5000:
                if len(set(anndata.obs[save_rep])) < N:
                    minn = x
                    x = (maxn + x) / 2
                if len(set(anndata.obs[save_rep])) > N:
                    maxn = x
                    x = (minn + x) / 2
                if prev_x == x:
                    if len(set(anndata.obs[save_rep])) < N:
                        x = x * np.random.uniform(low=1.5, high=4.0)
                    elif len(set(anndata.obs[save_rep])) > N:
                        x = x / np.random.uniform(low=1.5, high=4.0)
                sc.tl.leiden(adata=anndata, adjacency=graph, resolution=x, key_added=save_rep, random_state=rs,
                             directed=False, use_weights=False)
                c += 1
                prev_x = x
        else:
            sc.tl.leiden(adata=anndata, adjacency=graph, key_added=save_rep, directed=False, random_state=rs,
                         use_weights=False)
    elif method == 'agglomerative':

        clust = AgglomerativeClustering(n_clusters=len(set(optimisation_labs)), affinity='precomputed',
                                                        linkage=linkage).fit(X=prox_matrix)
        cluster_results = pd.Categorical(
            values=np.array(clust.labels_).astype('U'),
            categories=natsorted(np.unique(np.array(clust.labels_)).astype('U')))
        anndata.obs[save_rep] = cluster_results
    else:
        print("Select between leiden and agglomerative clustering")


def cluster_stats(data, labels_true, labels_pred, prox_metric, eval_list):
    labels_true = labels_true.cat.codes.astype(str).astype('category')
    """
    Performs 

    Parameters
    ----------

    Returns
    -------
    dict: dictionary
        dictionary object containing the results for each clustering evaluation metric:

            ARI :
            Adjusted Rand Index - RI computes similarity by considering all pairs of samples & counting pairs that are
            assigned in the same or different clusters in the predicted and true clusterings. Adjusted Rand Index has
            been corrected for chance. Ranges between 0 (random labeling) to 1 (identical). Takes the true and predicted
            labels as input.
            #https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html

            NMI:
            Normalised Mutual Information - Measures the similarity between two labels of the same Data, normalised to
            a 0 (no mutual information) to 1 (perfect correlation) scale. Takes the true and predicted
            labels as input.
            #https://scikit-learn.org/stable/modules/generated/sklearn.metrics.mutual_info_score.html#sklearn.metrics.mutual_info_score
            #https://scikit-learn.org/stable/modules/generated/sklearn.metrics.normalized_mutual_info_score.html

            AMI:
            Adjusted Mutual Information - Measures the similarity between two labels of the same Data, adjusted for
            chance to ensure number of clusters does not skew result. Upper limit is 1 (identical labelling) and whilst
            random labelling has an expected value of 0, may become negative. Takes the true and predicted
            labels as input, as well as method for calculating of the normalizer in denominator - default is arithmetic.
            #https://scikit-learn.org/stable/modules/generated/sklearn.metrics.mutual_info_score.html#sklearn.metrics.mutual_info_score
            #https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_mutual_info_score.html#sklearn.metrics.adjusted_mutual_info_score

            FM:
            Fowlkes-Mallows score -  Measures similarity of two clusterings of a set of points, defined as the
            geometric mean between of the precision and recall. Ranges between 0 and 1, with a higher value indicating
            greater similarity between two clusters. Takes the true and predicted labels as input.
            #https://scikit-learn.org/stable/modules/generated/sklearn.metrics.fowlkes_mallows_score.html

            SIL:
            Silhouette coefficient - Compute the mean Silhouette Coefficient of all samples, where Silhouette
            Coefficient is a measure of how well samples are clustered with samples that are similar to themselves.
            Score is bounded between -1 for incorrect clustering and +1 for highly dense clustering (i.e higher
            score = better defined clusters). Values near 0 indicate overlapping clusters. Takes the pairwise distance
            matrix (when metric='precomputed') and a set of labels as input.
            #https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_score.html#sklearn.metrics.silhouette_score

            VRC:
            Variance Ratio Criterion (also called Calinski-Harabasz Index) - The ratio of the sum of between-clusters
            dispersion and of within-cluster dispersion for all clusters, where dispersion is defined as the sum of
            distances squared. A higher score indicates clusters are dense and well separated. Takes a A list of
            n_features-dimensional Data points and a set of labels as input.

            DBS:
            Davies-Bouldin Index - Signifies the average ‘similarity’ between clusters, where the similarity is a
            measure that compares distance between clusters with the size of the clusters themselves. A lower score
            indicates better separation between clusters. LIMITED TO EUCLIDEAN SPACE.

            Jac:
            Jaccards Index - The size of the intersection divided by the size of the union of two label sets.
    """
    ARI = ms.adjusted_rand_score(labels_true=labels_true, labels_pred=labels_pred)
    AMI = ms.adjusted_mutual_info_score(labels_true=labels_true, labels_pred=labels_pred)
    SIL = ms.silhouette_score(X=data.obsp[f'{prox_metric}_matrix'], labels=labels_pred, metric='precomputed',
                              sample_size=None, random_state=0)  # sample size parameter? random state?
    VRC = ms.calinski_harabasz_score(X=data.X.todense(), labels=labels_pred)
    DBS = ms.davies_bouldin_score(X=data.X.todense(), labels=labels_pred)
    HS = ms.homogeneity_score(labels_true=labels_true, labels_pred=labels_pred)
    FM = ms.fowlkes_mallows_score(labels_true=labels_true, labels_pred=labels_pred,
                                  sparse=False)  # Should sparse be true?
    MCC = ms.matthews_corrcoef(y_true=labels_true, y_pred=labels_pred)
    CS = ms.completeness_score(labels_true=labels_true, labels_pred=labels_pred)
    # CGM = ms.cluster.contingency_matrix(labels_true=labels_true, labels_pred=labels_pred)
    CFM = ms.cluster.pair_confusion_matrix(labels_true=labels_true, labels_pred=labels_pred)
    Jac = (CFM[1, 1] / (CFM[1, 1] + CFM[0, 1] + CFM[1, 0]))
    NMI = ms.normalized_mutual_info_score(labels_true=labels_true,
                                          labels_pred=labels_pred)  # has different averaging methods
    PSI = genieclust.compare_partitions.pair_sets_index(labels_true, labels_pred)  # has different averaging methods
    ev_dict = {
        'ARI': ARI,
        'AMI': AMI,
        'SIL': SIL,
        'VRC': VRC,
        'DBS': DBS,
        'Jac': Jac,
        'HS': HS,
        'MCC': MCC,
        'CS': CS,
        'FM': FM,
        # 'CGM': CGM,
        # 'CFM': CFM,
        'NMI': NMI,
        'PSI': PSI
    }
    eval_dict = ev_dict.fromkeys(eval_list)
    for e in eval_dict.keys():
        eval_dict[e] = ev_dict.get(e)
    return eval_dict


def reform(data, levels=None):
    """
    Performs

    Parameters
    ----------

    Returns
    -------

    """
    if levels is None:
        levels = ['Eval', 'Metric', 'k_val']
    for level1 in data:
        for level2 in data[level1]:
            for level3 in data[level1][level2]:
                values = np.array(data[level1][level2][level3])
                try:
                    values[values == -100] = np.nan
                except:
                    print("Warning: If you get this there may be an issue in your calculations. "
                          "Please contact the authors :) ")
                    print(level1, level2, level3)
                    print(values)
                data[level1][level2][level3] = values
    reform_dict = {(level1_key, level2_key, level3_key): values
                   for level1_key, level2_dict in data.items()
                   for level2_key, level3_dict in level2_dict.items()
                   for level3_key, values in level3_dict.items()}
    results = pd.DataFrame(reform_dict).T
    results.index.set_names(levels, inplace=True)
    return results


def merge_evaluations(anndata, metrics, k_vals, eval_metrics, clustering_method, data_dir, random_seeds,
                      save_as_csv=False):
    """
    Performs 

    Parameters
    ----------

    Returns
    -------

    """
    all_results = dict.fromkeys(eval_metrics, )
    cols = [str(k) + 'K_stats' for k in k_vals]
    for eval_m in eval_metrics:
        all_results[eval_m] = dict.fromkeys(metrics, )
        for metric in metrics:
            df = pd.DataFrame(columns=cols, index=random_seeds)
            for seed in random_seeds:
                dat = sc.read_h5ad(os.path.join(data_dir, anndata + f'_Clust_{seed}.h5ad'))
                results = pd.DataFrame.from_dict({(i, j): dat.uns[i][j]
                                                  for i in dat.uns.keys()
                                                  for j in dat.uns[i].keys()}, orient='index')
                names = ['k_value', 'Metric']
                results.index.set_names(names, inplace=True)
                for k in df.columns:
                    if (len(set(dat.obs[f'{metric}_{clustering_method}{k[0:-7]}k']))) == (len(set(dat.obs['labels_truth']))):
                        df.at[seed, k] = results.loc[pd.IndexSlice[f'{k}', f'{metric}'], f'{eval_m}']
                    else:
                        df.at[seed, k] = np.nan
            all_results[eval_m][metric] = df
    filename = os.path.join(data_dir, anndata + '_Evaluations')
    outfile = open(filename, 'wb')
    pickle.dump(all_results, outfile)
    outfile.close()
    if save_as_csv is True:
        csv_results = reform(all_results)
        csv_results.to_csv(os.path.join(data_dir, f'{anndata}_Evaluations_Complete.csv'))
    return all_results


def kruskal_test(data, anndata, metrics, eval_metrics, k_vals, random_seeds, data_dir):
    """
    Performs 

    Parameters
    ----------

    Returns
    -------

    """
    Kruskal_results = pd.DataFrame(columns=list(metrics), index=eval_metrics)
    data = data.assign(kruskal=np.nan)
    for eval_m in Kruskal_results.index.values:
        sub = data.loc[pd.IndexSlice[eval_m]]
        for metric in Kruskal_results.columns:
            K_list = []
            k_stats = [str(k) + 'K_stats' for k in k_vals]
            for i, k_val in enumerate(k_stats):
                kdat = sub.xs((f"{metric}", f"{k_val}"))
                kdat = kdat.dropna()
                if len(kdat) < len(random_seeds):
                    print(f'{k_val} for {metric} in {eval_m} contains only {len(kdat)} samples')
                K_list.append(kdat)
            args = list(map(np.asarray, K_list))
            alldata = np.concatenate(args)
            ranked = rankdata(alldata)
            ties = tiecorrect(ranked)
            if ties == 0:
                krusk = 1.0
            else:
                krusk = kruskal(*K_list)
                krusk = krusk[1]
            Kruskal_results.at[eval_m, metric] = krusk
            data.loc[(eval_m, metric), 'kruskal'] = krusk
    Kruskal_results.to_csv(os.path.join(data_dir, f'{anndata}_Kruskal_results.csv'))
    print(Kruskal_results)
    return data


def dunn_test(data, anndata, metrics, eval_metrics, k_vals, random_seeds, data_dir, printout=True, visualise=False,
              save=False):
    """
    Performs 

    Parameters
    ----------

    Returns
    -------

    """
    Dunn_dict = dict.fromkeys(eval_metrics, )
    for eval_m in Dunn_dict.keys():
        sub = data.loc[pd.IndexSlice[eval_m]]
        mets = sub.index.get_level_values(0).unique()
        met_dict = dict.fromkeys(mets, )
        for metric in met_dict.keys():
            df = sub.xs((f"{metric}"))
            df = df.reset_index()
            df_melt = pd.melt(df, id_vars=['k_val'], value_vars=list(random_seeds))
            df_melt['k_val'] = pd.Categorical(df_melt['k_val'], categories=df_melt['k_val'].unique(), ordered=True)
            met_dict[metric] = sph.posthoc_dunn(df_melt, val_col='value', group_col='k_val', p_adjust='holm')
            if printout is True:
                print(f'Dunn test results for {metric} in {anndata}: {met_dict[metric]}')
            if visualise is True:
                mask = np.triu(np.ones_like(met_dict[metric], dtype=bool))
                f, ax = plt.subplots(figsize=(11, 9))
                cmap = sns.diverging_palette(230, 20, as_cmap=True)
                ax = sns.heatmap(met_dict[metric], mask=mask, cmap=cmap, vmin=0.0, vmax=0.1, center=0.05,
                                 square=True, linewidths=.5, cbar_kws={"shrink": .5}, annot=True, annot_kws={"size": 7})
                ax.set_title(f'Dunn Post-Hoc for {metric} in {anndata}')
                plt.show()
                if save is True:
                    plt.savefig(f'{data_dir}{anndata}_{metric}_DunnResults.png', bbox_inches='tight')
        Dunn_dict[eval_m] = met_dict
    result = reform(Dunn_dict)
    result.to_csv(f'{data_dir}{anndata}_Dunn_Results.csv', header=True, index=True)
    return Dunn_dict


def load_pickle(data, data_dir, mode='mean', metrics=None, eval_metrics=None, k_vals=None):
    """
    Performs

    Parameters
    ----------

    Returns
    -------

    """
    infile = open(os.path.join(data_dir, data + '_Evaluations'), 'rb')
    dat = pickle.load(infile)
    infile.close()
    results = reform(dat)
    dataframe = None
    if mode == 'full':
        dataframe = results
    elif mode == 'mean':
        if eval_metrics is None:
            eval_metrics = results.index.get_level_values('Eval').unique()
        if metrics is None:
            metrics = results.index.get_level_values('Metric').unique()
        if k_vals is None:
            k_vals = results.index.get_level_values('k_val').unique()
        all_results = dict.fromkeys(k_vals, )
        for k in k_vals:
            eval_df = pd.DataFrame(columns=metrics)
            for metric in metrics:
                for eval_m in eval_metrics:
                    series = results.loc[pd.IndexSlice[f'{eval_m}', f'{metric}', f'{k}']]
                    series = series.dropna()
                    eval_df.at[eval_m, metric] = series.mean()
            all_results[k] = eval_df
        final_df = pd.DataFrame.from_dict({(i, j): all_results[i][j]
                                           for i in all_results.keys()
                                           for j in all_results[i].keys()}, orient='index')
        names = ['k_val', 'Metric']
        final_df.index.set_names(names, inplace=True)
        dataframe = final_df
    return dataframe


def vis_dataframe(datasets, datasets_dict, metrics, k_vals, eval_metrics):
    """
    Performs 

    Parameters
    ----------

    Returns
    -------

    """
    metrics_df = []
    scores_df = []
    kval_df = []
    eval_df = []
    datasetname_df = []
    structure_df = []  # Discrete or cont
    rarity_df = []  # Rare/ultra rare etc
    condition_df = []
    property_df = []

    measures = eval_metrics
    k_test = k_vals
    inc_metrics = metrics
    data_dict = {key: datasets_dict[key] for key in datasets}
    for dataset in data_dict.keys():
        for measure in measures:
            dat = data_dict[dataset][measure].reset_index().pivot(columns='k_val',
                                                                  index='Metric', values=measure)
            for k in k_test:
                for i, metric_label in enumerate(dat.index):
                    scores_df.append(dat[k][i])
                    metrics_df.append(metric_label)  # save the label
                    kval_df.append(k)
                    eval_df.append(measure)
                    datasetname_df.append(dataset)
                    dataset_list = dataset.split('_')
                    structure_df.append(dataset_list[0])  # Separate out the dataset
                    rarity_df.append(dataset_list[1])
                    condition_df.append(dataset_list[0] + ' ' + dataset_list[1])
                    prop_label = 'Base' if len(dataset_list) < 3 else dataset_list[2]
                    property_df.append(prop_label)

    vis_df = pd.DataFrame()
    vis_df['Metric'] = metrics_df
    vis_df['Performance'] = scores_df
    vis_df['k_value'] = kval_df
    vis_df['Evaluation_Metric'] = eval_df
    vis_df['Dataset'] = datasetname_df
    vis_df['Data_Structure'] = structure_df
    vis_df['Population_Balance'] = rarity_df
    vis_df['Condition'] = condition_df
    vis_df['Property'] = property_df
    return vis_df
