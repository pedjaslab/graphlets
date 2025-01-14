'''
Code to organize svml output from graphlet counting into a usable sparse matrix file
Ross stewart
Dec 2024
'''

import os
import glob
import pandas as pd
import numpy as np
import sys
from scipy.sparse import coo_matrix, csr_matrix, save_npz

assert len(sys.argv) == 3, f"Usage: python {sys.argv[0]} <path_to_svml_dir> <path_to_results_dir>"


'''
filepaths
'''

svml_results_dir = sys.argv[1] # location of graphlet counting output
results_dir = sys.argv[2] # directory to store results
pdb_ids_out_f = f"{results_dir}/pdb_ids.npy" # pdb id corresponding to every data point, useful if you are stratifying by protein. leave empty if you don't wish to save
features_out_f = f"{results_dir}/features.npz" # feature matrix file name
labels_out_f = f"{results_dir}/labels.npy" # label file name
os.makedirs(results_dir,exist_ok=True)


'''
read graphlet features from each `.svml` file
'''

df = pd.DataFrame()

residue_idx = 0 # the row of the sparse matrix
labels = []
pdb_ids = []
sparse_feats = {}
for svml in glob.glob(f'{svml_results_dir}/*.svml'):
    pdb_id = svml.split('/')[-1].split('.')[0]
    with open(svml,'r') as f:
        for line in f:
            vals = line.strip().split()
            label = int(vals[0])
            idx = int(vals[-1].replace('#',''))
            feats_and_counts  = vals[1:-1]
            for feat_and_count in feats_and_counts:
                feat,count = feat_and_count.split(':')
                # assert int(feat) not in sparse_feats
                sparse_feats[(residue_idx,int(feat))] = float(count)
            residue_idx += 1
            labels.append(label)
            pdb_ids.append(pdb_id)

if len(pdb_ids_out_f) != 0:
    np.save(pdb_ids_out_f, pdb_ids)
    print(f'pdb ids saved in {pdb_ids_out_f} with shape {np.array(pdb_ids).shape}')

np.save(labels_out_f, labels)
print(f'labels saved in {labels_out_f} with shape {np.array(labels).shape}')


'''
convert into `scipy.sparse` format
'''

rows,cols,vals = [],[],[]
for (row,col), val in sparse_feats.items():
    rows.append(row)
    cols.append(col)
    vals.append(val)

# assume the highest col index appears in the data
n_rows = max(rows) + 1
n_cols = max(cols) + 1

sparse_mat = coo_matrix((vals, (rows, cols)), shape=(n_rows, n_cols)).tocsr()
# print(sparse_mat)


'''
optional - shift feature indices if they do not start at 0
'''

# feature indices sometimes start at 6 million or something, shift them to start at 0 and increment by 1
def shift_feature_indices(X):
    unique_features = np.unique(X.indices)
    feature_mapping = {old_index: new_index for new_index, old_index in enumerate(unique_features)}
    new_indices = np.array([feature_mapping[old_index] for old_index in X.indices])
    new_csr_mat = csr_matrix((X.data, new_indices, X.indptr), shape=(X.shape[0], len(unique_features)))
    
    return new_csr_mat
    
sparse_mat = shift_feature_indices(sparse_mat)


'''
save feature file
'''

save_npz(features_out_f, sparse_mat)
print(f'sparse matrix saved as {features_out_f} with shape {sparse_mat.shape}')













