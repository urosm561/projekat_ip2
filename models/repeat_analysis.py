import pandas as pd
import numpy as np
import itertools

from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth, association_rules
from sklearn.decomposition import PCA
from collections import Counter
from matplotlib import pyplot as plt

from data_split import prepare_data

all_kmers = ["".join(p) for p in itertools.product(list("ACDEFGHIKLMNPQRSTVWY"), repeat = 2)]

def extract_kmers(seq, k = 2):
    kmers = [seq[i : i + k] for i in range(len(seq) - k + 1)] if len(seq) >= k else seq
    counts = Counter(kmers)
    vec = np.array([counts[kmer] for kmer in all_kmers], dtype = float)

    return vec

SEED = 561
data_path = "proteins/data/virus_protein_repeats.csv"

train_data, valid_data, test_data = prepare_data(data = data_path, seed = SEED, train_size = 0.8, test_size = 0.2)
print(train_data)