import numpy as np
import pandas as pd

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AMINO_ACID_TO_INDEX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

PROTVEC_DF = pd.read_csv("utils/protVec_100d_3grams.csv", sep = "\\t", header = None)
PROTVEC_DF = PROTVEC_DF.applymap(lambda x: str(x).strip('"') if isinstance(x, str) else x)
PROTVEC_DF.columns = ["3gram"] + [f"dim_{i}" for i in range(100)]
PROTVEC_DICT = {row["3gram"]: row[1:].values.astype(float) for _, row in PROTVEC_DF.iterrows()}

def get_kmers(seq, k = 3):
    return [seq[i : i + k] for i in range(len(seq) - k + 1)]

def one_hot_encode_amino_acid(amino_acid):
    vec = np.zeros(len(AMINO_ACIDS))
    if amino_acid in AMINO_ACID_TO_INDEX:
        vec[AMINO_ACID_TO_INDEX[amino_acid]] = 1
    return vec

def one_hot_encode_repeat(repeat, max_k = 7):
    encoding = []
    for i in range(max_k):
        if i < len(repeat):
            encoding.append(one_hot_encode_amino_acid(repeat[i]))
        else:
            encoding.append(np.zeros(len(AMINO_ACIDS)))
    return np.concatenate(encoding)

def protvec_encode_repeat(repeat, k = 3):
    kmers = get_kmers(repeat, k)
    vectors = []

    for kmer in kmers:
        vec = PROTVEC_DICT.get(kmer)
        if vec is not None:
            vectors.append(vec)
    
    if vectors:
        return np.mean(vectors, axis = 0)
    
    else:
        return np.zeros(100)
