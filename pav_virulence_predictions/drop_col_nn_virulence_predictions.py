import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import LabelEncoder
from sklearn.manifold import TSNE
from collections import Counter
import math
import random
from sklearn.model_selection import KFold
from sklearn.metrics import confusion_matrix
import torch.optim.lr_scheduler as lr_scheduler
from sklearn.metrics import roc_auc_score
import sys

file_path = sys.argv[1]
ogs_to_drop_file = sys.argv[2]

# read in virulence phenotypes data
df_phenotypes = pd.read_csv(file_path + "virulence_phenotypes.txt", sep='\t')
df_phenotypes.columns = ['isolate', 'group', 'rice_type', 'symptom_score', 'variety']
df_phenotypes = df_phenotypes.drop('rice_type', axis=1)
df_phenotypes = df_phenotypes.drop('group', axis=1)
# take median score of biological replicates
df_phenotypes = df_phenotypes.groupby(['isolate','variety'])['symptom_score'].median().reset_index()
# round down some of the samples with only two replicates
df_phenotypes = df_phenotypes.applymap(lambda x: math.floor(x) if isinstance(x, float) else x)
# subtract one from all to zero index ordinal scale
df_phenotypes = df_phenotypes.applymap(lambda x: x - 1 if isinstance(x, int) else x)
# pivot table and rename indices
df_phenotypes = df_phenotypes.pivot(index=['isolate'], columns='variety', values='symptom_score').reset_index()
df_phenotypes = df_phenotypes.rename_axis(None, axis=1)
df_phenotypes.set_index('isolate', inplace=True)
# this isolate was not properly sequenced so remove it from the data
df_phenotypes = df_phenotypes[df_phenotypes.index != 'IN0092']
# make sure to sort alphabetically to match with groups array
df_phenotypes = df_phenotypes.sort_index(axis=0)
# grab list of rice varieties
varieties = df_phenotypes.columns
# scale and center
df_phenotypes_untransformed = df_phenotypes
df_phenotypes = StandardScaler().fit_transform(df_phenotypes)
## read in PAV matrix for all isolates
df_pav = pd.read_csv(file_path + "pav_matrix.txt", sep=' ')
df_pav.columns = [elem.split('_')[0] for elem in df_pav.columns ]
df_pav = df_pav.transpose()
## this is the outgroup for orthogrouping so we don't need it
df_pav = df_pav[df_pav.index != 'NI907']
## drop fully conserved genes
df_pav = df_pav.loc[:, df_pav.var() != 0]
## normalize data and convert to tensor
scaler = StandardScaler()
X = scaler.fit_transform(df_pav)
X = torch.tensor(X, dtype=torch.float32)

class My_Network(nn.Module):
    def __init__(self, num_layers, hidden_size, input_dim, output_dim, dropout=False, dropout_perc=0.5):
        super().__init__()
        self.input_layer  = nn.Sequential(
            nn.Linear(input_dim, hidden_size),
            nn.BatchNorm1d(hidden_size),
        )
        self.output_layer = nn.Sequential(
            nn.Linear(hidden_size, output_dim),
            nn.Sigmoid()
        )
        self.hidden_layers = nn.Sequential()
        for i in range(num_layers):
          self.hidden_layers.add_module(f"hidden_layer_{i}", nn.Linear(hidden_size, hidden_size))
          if dropout:
            self.hidden_layers.add_module(f"dropout_{i}", nn.Dropout(dropout_perc))
          self.hidden_layers.add_module(f"activation_{i}", nn.ReLU())

    def forward(self, X):
      out = self.input_layer(X)
      out = self.hidden_layers(out)
      logits = self.output_layer(out)
      return logits

def train_nn(model, X_train, y_train, X_test, y_test, epochs=15, batch_size=32, lr=1e-3, verbose=False):
    """
    Q:  write the training loop following the schema shown above.

    Inputs
    - model: the model to be trained - a PyTorch nn.Module class object
    - X_train, y_train, X_val, y_val: training and validation data
    - epochs: num epochs, or the number of times we want to run through the entire training data
    - batch_size: number of data points per batch
    - lr: learning rate
    - optimizer: optimizer used

    Outputs
    - losses: a list of losses
    - accuracies: a list of validation accuracies
    - train_accs: a list of training accuracies
    """

    N, D = X_train.shape

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    loss_fn = nn.BCELoss()

    losses = []

    # scheduler = lr_scheduler.LinearLR(optimizer, start_factor=1.0, end_factor=0.01, total_iters=epochs)

    for epoch in range(epochs):
      logits = model(X_train)

      loss = loss_fn(logits, y_train)

      optimizer.zero_grad()
      loss.backward()
      optimizer.step()

      training_loss = loss

      testing_loss = loss_fn(model(X_test), y_test)

      # print epoch, loss, and current test accuracy
      if verbose:
        print(f"Epoch {epoch}:\t training loss {training_loss} & testing loss {testing_loss}")

    return losses

## first convert the phenotypes to disease and no disease
def convert_to_binary(x):
    if x >= 0 and x <= 2:
        return 0
    elif x >= 3 and x <= 5:
        return 1
    else:
        return x
    

df_phenotypes_binary = df_phenotypes_untransformed.applymap(convert_to_binary)


layers = 1
nodes = 50
epochs = 500
lr = 1e-4
dropout = False
dropout_perc = 0.5
y = torch.tensor(df_phenotypes_binary.values, dtype=torch.float32)
average_fprs = []
average_tprs = []
average_f1s = []
for k in [0,1,2]:
  kfold_fprs = []
  kfold_tprs = []
  kfold_f1s = []
  kf = KFold(n_splits=5, shuffle=True)
  for train_index, test_index in kf.split(X):
    X_train, X_test, y_train, y_test = X[train_index], X[test_index], y[train_index], y[test_index]
    my_model = My_Network(1,50,X_train.shape[1], y.shape[1], dropout=False)
    losses = train_nn(my_model, X_train, y_train, X_test, y_test, epochs=epochs, lr=lr, verbose=False)
    ## evaluate ##
    my_model.eval()
    preds = my_model(X_test)
    tns = 0
    fps = 0
    fns = 0
    tps = 0
    for i in range(len(preds)):
      tn, fp, fn, tp = confusion_matrix(y_test[i], torch.round(preds[i]).detach().numpy()).ravel()
      tns += tn
      fps += fp
      fns += fn
      tps += tp
    fpr = fps / (fps + tns)  # false positive rate
    tpr = tps / (tps + fns)  # true positive rate
    kfold_fprs.append(fpr)
    kfold_tprs.append(tpr)
    precision = tps / (tps+fps)
    recall = tps / (tps+fns)
    f1 = 2*tps/(2*tps+fps+fns)
    kfold_f1s.append(f1)
  average_tpr = sum(kfold_tprs)/len(kfold_tprs)
  average_fpr = sum(kfold_fprs)/len(kfold_fprs)
  average_f1 = sum(kfold_f1s)/len(kfold_f1s)
  average_tprs.append(average_tpr)
  average_fprs.append(average_fpr)
  average_f1s.append(average_f1)
final_fpr = sum(average_fprs)/len(average_fprs)
final_tpr = sum(average_tprs)/len(average_tprs)
final_f1 = sum(average_f1s)/len(average_f1s)

# set original f1
full_model_f1 = final_f1

with open(ogs_to_drop_file, 'r') as f:
    ogs_to_drop = [line.strip() for line in f.readlines()]

layers = 1
nodes = 50
epochs = 500
lr = 1e-4
dropout = False
dropout_perc = 0.5
y = torch.tensor(df_phenotypes_binary.values, dtype=torch.float32)
for og in ogs_to_drop:
  df_pav_subset = df_pav.drop(og, axis=1).copy(deep=True)
  ## normalize data and convert to tensor as before
  scaler = StandardScaler()
  X = scaler.fit_transform(df_pav)
  X = torch.tensor(X, dtype=torch.float32)
  average_f1s = []
  for k in [0,1,2]:
    kfold_f1s = []
    kf = KFold(n_splits=5, shuffle=True)
    for train_index, test_index in kf.split(X):
      X_train, X_test, y_train, y_test = X[train_index], X[test_index], y[train_index], y[test_index]
      my_model = My_Network(1,50,X_train.shape[1], y.shape[1], dropout=False)
      losses = train_nn(my_model, X_train, y_train, X_test, y_test, epochs=epochs, lr=lr, verbose=False)
      ## evaluate ##
      my_model.eval()
      preds = my_model(X_test)
      tns = 0
      fps = 0
      fns = 0
      tps = 0
      for i in range(len(preds)):
        tn, fp, fn, tp = confusion_matrix(y_test[i], torch.round(preds[i]).detach().numpy()).ravel()
        tns += tn
        fps += fp
        fns += fn
        tps += tp
      fpr = fps / (fps + tns)  # false positive rate
      tpr = tps / (tps + fns)  # true positive rate
      precision = tps / (tps+fps)
      recall = tps / (tps+fns)
      f1 = 2*tps/(2*tps+fps+fns)
      kfold_f1s.append(f1)
    average_f1 = sum(kfold_f1s)/len(kfold_f1s)
    average_f1s.append(average_f1)
  subset_f1 = sum(average_f1s)/len(average_f1s)
  diff_f1 = subset_f1 - full_model_f1
  print(og, diff_f1)