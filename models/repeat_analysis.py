import pandas as pd
import numpy as np
import itertools

from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth, association_rules
from sklearn.decomposition import PCA
from collections import Counter
from matplotlib import pyplot as plt

from preprocess_data import preprocess_data, PositionalEncoder, PAD

def ngrams(s: str, n: int):
    return {s[i : i + n] for i in range(len(s) - n + 1)} if len(s) >= n else set()

def build_transactions(data: pd.DataFrame, width: int = 7):
    pe = PositionalEncoder(width = width)
    pos_data = pe.transform(data[["repeat"]])

    baskets = []
    for (idx, row), (_, orig) in zip(pos_data.iterrows(), data.iterrows()):
        pos_items = [f"P{i+1}={row[f'p{i+1}']}" for i in range(width) if row[f"p{i+1}"] != PAD]
        raw = str(orig["repeat"])[-width:]
        true_len = len(raw)

        items = pos_items[:]

        if true_len <= 3: items.append("LEN=<=3")
        elif true_len <=5: items.append("LEN=4-5")
        else: items.append("LEN=6-7")

        for bg in ngrams(raw, 2): items.append(f"BG={bg}")
        # for tg in ngrams(raw, 3): items.append(f"TG={tg}")

        for ch in set(raw):
            items.append(f"AA={ch}")

        items.append(f"VIRUS={orig['virus_name']}")

        baskets.append(sorted(set(items)))
    
    te = TransactionEncoder()
    item_set = te.fit_transform(baskets)
    item_set_df = pd.DataFrame(data = item_set, columns = te.columns_)

    return item_set_df

def mine_rules(data: pd.DataFrame, width = 7, min_confidence: float = 0.8, min_lift: float = 1.3, min_support: float = 0.01, min_antecedents = 2, max_len = 5):
    item_set_df = build_transactions(data, width)

    sets = fpgrowth(item_set_df, min_support = min_support, use_colnames = True, max_len = max_len, verbose = 0)
    rules = association_rules(sets, metric = "confidence", min_threshold = min_confidence)

    rules = rules[rules["consequents"].apply(lambda s: len(s) == 1 and next(iter(s)).startswith("VIRUS="))]
    rules = rules[rules["antecedents"].apply(len) >= min_antecedents]
    rules = rules[rules["lift"] >= min_lift].sort_values(["lift", "confidence", "support"], ascending = False)

    return rules

if __name__ == "__main__":
    SEED = 561
    data_path = "proteins/data/virus_protein_repeats.csv"

    train_data, valid_data, test_data = preprocess_data(data = data_path, seed = SEED, train_size = 0.8, test_size = 0.2)
    print(train_data)
    print(mine_rules(train_data))