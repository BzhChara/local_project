import pandas as pd
import numpy as np
import esm
import torch.nn as nn
import torch
import argparse
from sklearn.feature_selection import VarianceThreshold
# protein
model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()

batch_converter = alphabet.get_batch_converter()
model.eval()
#number
protein = pd.read_csv('protein.csv')
seq = protein.iloc[:, 0].astype(str).tolist()

seq = ["".join([char.upper() for char in s if not char.isdigit()]) for s in seq]
#print(seq)

chunk_size = 1
my_data =[]

for i in range(0, len(seq), chunk_size):
    sequences = seq[i:i + chunk_size]
    data_list = []
    for j in range(0,len(sequences)):
        data_list.append((1, sequences[j][:len(sequences[j])]))

    batch_labels, batch_strs, batch_tokens = batch_converter(data_list)#padding
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)#real lens

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]
    sequence_representations = []
    for i, tokens_len in enumerate(batch_lens):
        sequence_representations.append(token_representations[i, 1 : tokens_len - 1].mean(0))
    my_data.append(sequence_representations)
protein_all_data =[]
for i in my_data:
    for j in i:
        j =j.tolist()
        protein_all_data.append(j)
protein_data = np.array(protein_all_data)
protein_feature = pd.DataFrame(protein_data)
protein_feature.to_csv('protein_features.csv', index=False)