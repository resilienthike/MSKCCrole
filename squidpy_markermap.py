import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
import squidpy as sq

sc.logging.print_header()
print(f"squidpy=={sq.__version__}")

from markermap.vae_models import MarkerMap, train_model
from markermap.utils import (
    new_model_metrics,
    plot_confusion_matrix,
    split_data,
    plot_umap_embedding
)

from markermap import vae_models
print(vae_models.__file__)
from markermap import utils
print(utils.__file__)

def dataset_squidpy(batch_size=64):
  #get data
  adata = sq.datasets.visium_hne_adata()
  group_by = 'annotation'
  adata.obs[group_by] = adata.obs['cluster']
  # adata.obs[group_by] = adata.obsm['spatial']
  print(adata)
  print(adata.X.shape)
  print(adata.obs['annotation'].shape)
  print(np.unique(adata.obs['annotation']))
  adata.X = adata.X.toarray()

  # we will use 70% training data, 10% vaidation data during the training of the marker map, then 20% for testing
  train_indices, val_indices, test_indices = split_data(
    adata.X,
    adata.obs[group_by],
    [0.7, 0.1, 0.2],
  )
  train_val_indices = np.concatenate([train_indices, val_indices])

  train_dataloader, val_dataloader = MarkerMap.prepareData(
    adata,
    train_indices,
    val_indices,
    group_by,
    None, #layer, just use adata.X
    batch_size=batch_size,
  )

  return adata, train_dataloader, val_dataloader, train_indices, val_indices, test_indices, train_val_indices

def squidpy_markermap(
  params={
    'z_size': 512,
    'hidden_layer_size': 2048,
    'k': 100,
  },
  adata=None,
  train_dataloader=None,
  val_dataloader=None,
  verbose=False,
):
  print(params)
  supervised_marker_map = MarkerMap(
    adata.X.shape[1],
    params['hidden_layer_size'],
    params['z_size'],
    len(adata.obs['annotation'].unique()),
    params['k'],
    loss_tradeoff=0,
  )

  train_model(supervised_marker_map, train_dataloader, val_dataloader, verbose=verbose)
  return supervised_marker_map


def eval_squidpy(supervised_marker_map, adata, train_val_indices, test_indices):
  misclass_rate, test_rep, cm, Y_pred = new_model_metrics(
    adata[train_val_indices, :].X,
    adata[train_val_indices, :].obs['annotation'],
    adata[test_indices, :].X,
    adata[test_indices, :].obs['annotation'],
    markers = supervised_marker_map.markers().clone().cpu().detach().numpy(),
    X = adata.X
  )

  print('Misclass_rate:', misclass_rate)
  print('f1-score:', test_rep['weighted avg']['f1-score'])
  print('Number of genes for training:', adata.X.shape[1])
  plot_confusion_matrix(cm, adata.obs['annotation'].unique())
  return Y_pred

# def get_results

if __name__ == "__main__":
  print('hi')