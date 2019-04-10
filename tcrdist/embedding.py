import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import palettable

from sklearn.decomposition import PCA, KernelPCA
from sklearn.manifold import Isomap, LocallyLinearEmbedding, TSNE, MDS

"""Embedding of pairwise-distance matrices."""

__all__ = [ 'embedDistanceMatrix',
            'plotEmbedding',
            'clusteredScatter']

def embedDistanceMatrix(dmatDf, method='kpca', n_components=2, **kwargs):
    """Two-dimensional embedding of sequence distances in dmatDf,
    returning Nx2 x,y-coords: tsne, isomap, pca, mds, kpca, sklearn-tsne"""
    if isinstance(dmatDf, pd.DataFrame):
        dmat = dmatDf.values
    else:
        dmat = dmatDf

    if method == 'isomap':
        isoObj = Isomap(n_neighbors=10, n_components=n_components)
        xy = isoObj.fit_transform(dmat)
    elif method == 'mds':
        mds = MDS(n_components=n_components,
                  max_iter=3000,
                  eps=1e-9,
                  random_state=15,
                  dissimilarity="precomputed",
                  n_jobs=1)
        xy = mds.fit(dmat).embedding_
        rot = PCA(n_components=n_components)
        xy = rot.fit_transform(xy)
    elif method == 'pca':
        pcaObj = PCA(n_components=None)
        xy = pcaObj.fit_transform(dmat)[:, :n_components]
    elif method == 'kpca':
        pcaObj = KernelPCA(n_components=dmat.shape[0], kernel='precomputed', eigen_solver='dense')
        try:
            gram = dist2kernel(dmat)
        except:
            print('Could not convert dmat to kernel for KernelPCA; using 1 - dmat/dmat.max() instead')
            gram = 1 - dmat / dmat.max()
        xy = pcaObj.fit_transform(gram)[:, :n_components]
    elif method == 'lle':
        lle = LocallyLinearEmbedding(n_neighbors=30, n_components=n_components, method='standard')
        xy = lle.fit_transform(dist)
    elif method == 'sklearn-tsne':
        tsneObj = TSNE(n_components=n_components, metric='precomputed', random_state=0, perplexity=kwargs['perplexity'])
        xy = tsneObj.fit_transform(dmat)
    else:
        print('Method unknown: %s' % method)
        return

    assert xy.shape[0] == dmatDf.shape[0]
    xyDf = pd.DataFrame(xy[:, :n_components], index=dmatDf.index, columns=np.arange(n_components))
    if method == 'kpca':
        """Not sure how negative eigenvalues should be handled here, but they are usually
        small so it shouldn't make a big difference"""
        xyDf.explained_variance_ = pcaObj.lambdas_[:n_components]/pcaObj.lambdas_[pcaObj.lambdas_>0].sum()
    return xyDf

def plotEmbedding(dmatDf,
                  xyDf=None,
                  labels=None,
                  method='kpca',
                  plotLabels=False,
                  plotDims=[0, 1],
                  weights=None,
                  txtSize='large',
                  alpha=0.8,
                  sz=50,
                  mxSz=500,
                  markers=None,
                  plotLegend=True,
                  colors=None,
                  markerLabels=None,
                  continuousLabel=False):
    """Two-dimensional plot of embedded distance matrix, colored by labels"""
    
    if labels is None:
        labels = np.zeros(dmatDf.shape[0])

    assert dmatDf.shape[0] == dmatDf.shape[1]
    assert labels.shape[0] == dmatDf.shape[0]

    uLabels = freqSort(labels)
    
    if xyDf is None:
        xyDf = embedDistanceMatrix(dmatDf, method=method, n_components=np.max(plotDims) + 1)
    
    clusteredScatter(xyDf,
                     labels=labels,
                     plotDims=plotDims,
                     weights=weights,
                     alpha=alpha,
                     sz=sz,
                     mxSz=mxSz,
                     markerLabels=markerLabels,
                     markers=markers,
                     continuousLabel=continuousLabel,
                     colors=colors)
   
    if plotLabels:
        annotationParams = dict(xytext=(0, 5), textcoords='offset points', size=txtSize)
        for coli, col in enumerate(dmatDf.columns):
            if plotLabels:
                axh.annotate(col, xy=(xyDf.loc[col, plotDims[0]], xyDf.loc[col, plotDims[1]]), **annotationParams)

    if len(uLabels) > 1 and plotLegend:
        if labels is None and markerLabels is None:
            pass
        else:
            if labels is None:
                legTit = markerLabels.name
            elif markerLabels is None:
                legTit = abels.name
            else:
                legTit = '%s | %s' % (labels.name, markerLabels.name)
            plt.legend(loc=0, title=legTit)
    if hasattr(xyDf, 'explained_variance_'):
        plt.xlabel('KPCA %1.0f (%1.0f%% variance explained)' % (plotDims[0]+1, 100*xyDf.explained_variance_[plotDims[0]]))
        plt.ylabel('KPCA %1.0f (%1.0f%% variance explained)' % (plotDims[1]+1, 100*xyDf.explained_variance_[plotDims[1]]))
    else:
        plt.xlabel('KPCA %1.0f' % (plotDims[0]+1))
        plt.ylabel('KPCA %1.0f' % (plotDims[1]+1))
    plt.show()
    return xyDf

def clusteredScatter(xyDf,
                     labels=None,
                     plotDims=[0, 1],
                     weights=None,
                     alpha=0.8,
                     sz=50,
                     mxSz=500,
                     markerLabels=None,
                     markers=None,
                     colors=None,
                     continuousLabel=False):
    """Produce a scatter plot with axes, shaded by values in labels and with specified markers

    Parameters
    ----------
    xyDf : pd.DataFrame
        One column for each plot dimension specified
    labels : pd.Series
        Holds categorical or continuous metadata used for coloring
    plotDims : list len 2
        Specifies the columns in xyDf for plotting on X and Y axes
    weights : pd.Series
        Relative weights that are mapped for sizing each symbol
    alpha : float
        Transparency of each point
    sz : float
        Size of each symbol or minimum size of a symbol if there are weights
    mxSz : float
        Maximum size of a symbol if there are weights
    markerLabels : pd.Series
        Holds categorical labels used for plotting different symbols
    markers : list
        List of marker symbols to use. Defaults to: "ov8sp*hDPX"
    colors : list
        List of colors to use for categorical labels or a colormap
        for a continuousLabel. Defaults to colors from Set1 or the YlOrRd colormap
    continuousLabel : bool
        Indicates whether labels are categorical or continuous"""
    if weights is None:
        sVec = sz * pd.Series(np.ones(xyDf.shape[0]), index=xyDf.index)
    else:
        sVec = weights * mxSz + sz
    
    if labels is None:
        labels = pd.Series(np.ones(xyDf.shape[0]), index=xyDf.index)
        useColors = False
    else:
        useColors = True

    if continuousLabel:
        if not colors is None:
            cmap = colors
        else:
            cmap = palettable.colorbrewer.sequential.YlOrRd_9.mpl_colormap
    else:
        cmap = None
        uLabels = freqSort(labels)
        
        if colors is None:
            nColors = min(max(len(uLabels), 3), 9)
            colors = palettable.colorbrewer.get_map('Set1', 'Qualitative', nColors).mpl_colors
        elif isinstance(colors, pd.Series):
            colors = colors[uLabels].values

    if markerLabels is None:
        markerLabels = pd.Series(np.ones(xyDf.shape[0]), index=xyDf.index)
        useMarkers = False
    else:
        useMarkers = True

    uMLabels = freqSort(markerLabels)
    
    if markers is None:
        nMarkers = len(uMLabels)
        markers = ['o', 'v', '8', 's', 'p', '*', 'h', 'D', 'P', 'X'][:nMarkers]
    elif isinstance(markers, pd.Series):
        markers = markers[uMLabels].values

    figh = plt.gcf()
    plt.clf()
    axh = figh.add_axes([0.05, 0.05, 0.9, 0.9])
    axh.patch.set_facecolor('white')
    # axh.axis('off')
    figh.patch.set_facecolor('white')

    if continuousLabel:
        for mi, m in enumerate(uMLabels):
            ind = markerLabels == m
            if useMarkers:
                labS = '%s (N=%d)' % (m, ind.sum())
            else:
                labS = None
            
            plt.scatter(xyDf.loc[ind, plotDims[0]],
                        xyDf.loc[ind, plotDims[1]],
                        marker=markers[mi],
                        s=sVec.loc[ind],
                        alpha=alpha,
                        c=labels.loc[ind],
                        label=labS,
                        cmap=cmap)
    else:
        for vi, v in enumerate(uLabels):
            for mi, m in enumerate(uMLabels):
                ind = (labels == v) & (markerLabels == m)
                if useMarkers and useColors:
                    labS = '%s|%s (N=%d)' % (v, m, ind.sum())
                elif useColors:
                    labS = '%s (N=%d)' % (v, ind.sum())
                elif useMarkers:
                    labS = '%s (N=%d)' % (m, ind.sum())
                else:
                    labS = None
                
                plt.scatter(xyDf.loc[ind, plotDims[0]],
                            xyDf.loc[ind, plotDims[1]],
                            marker=markers[mi],
                            s=sVec.loc[ind],
                            alpha=alpha,
                            c=[colors[vi % len(colors)], ] * ind.sum(),
                            label=labS,
                            cmap=cmap)
        axh.set_xticks(())
    axh.set_yticks(())
    plt.show()

def dist2kernel(dmat):
    """Convert a distance matrix into a similarity kernel
    for KernelPCA or kernel regression methods.

    Implementation of D2K in MiRKAT, Zhao et al., 2015

    Parameters
    ----------
    dmat : ndarray shape (n,n)
        Pairwise-distance matrix.

    Returns
    -------
    kernel : ndarray shape (n,n)"""

    n = dmat.shape[0]
    I = np.identity(n)
    """m = I - dot(1,1')/n"""
    m = I - np.ones((n, n))/float(n)
    kern = -0.5 * np.linalg.multi_dot((m, dmat**2, m))

    if isinstance(dmat, pd.DataFrame):
        return pd.DataFrame(kern, index=dmat.index, columns=dmat.columns)
    else:
        return kern

def freqSort(labels):
    freq = {}
    for lab in labels:
        try:
            freq[lab] += 1
        except KeyError:
            freq[lab] = 1

    uLabels = sorted(np.unique(labels), key=freq.get, reverse=True)
    return uLabels