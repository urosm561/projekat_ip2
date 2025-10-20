import pandas as pd

from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth, association_rules
from typing import Tuple

from preprocess_data import preprocess_data, PositionalEncoder, Metrics, PAD, calculate_metrics, split_into_proteins, rebalance_data, log_training

def ngrams(s: str, n: int):
    return {s[i : i + n] for i in range(len(s) - n + 1)} if len(s) >= n else set()

def build_baskets(data: pd.DataFrame, max_aa_length: int = 7, include_virus: bool = True):
    pe = PositionalEncoder(width = max_aa_length)
    pos_data = pe.transform(data[["repeat"]])

    baskets = []
    for (idx, row), (_, orig) in zip(pos_data.iterrows(), data.iterrows()):
        pos_items = [f"P{i+1}={row[f'p{i+1}']}" for i in range(max_aa_length) if row[f"p{i+1}"] != PAD]
        raw = str(orig["repeat"])[-max_aa_length:]
        true_len = len(raw)

        items = pos_items[:]

        if true_len <= 3: items.append("LEN=<=3")
        elif true_len <=5: items.append("LEN=4-5")
        else: items.append("LEN=6-7")

        for bg in ngrams(raw, 2): items.append(f"BG={bg}")
        # for tg in ngrams(raw, 3): items.append(f"TG={tg}")

        for ch in set(raw):
            items.append(f"AA={ch}")

        if include_virus:
            items.append(f"VIRUS={orig['virus_name']}")

        baskets.append(sorted(set(items)))

    return baskets

def build_transactions(data: pd.DataFrame, max_aa_length: int = 7, include_virus: bool = True):
    baskets = build_baskets(data, max_aa_length, include_virus)
    
    te = TransactionEncoder()
    item_set = te.fit_transform(baskets)
    item_set_df = pd.DataFrame(data = item_set, columns = te.columns_)

    return item_set_df

def train_ar_classifier(
        train_data: pd.DataFrame, 
        test_data: pd.DataFrame,
        seed: int,
        max_aa_length: int = 7, 
        min_confidence: float = 0.7, 
        min_lift: float = 1.3, 
        min_support: float = 0.001, 
        min_antecedents: int = 2, 
        max_antecedents: int = 6,
        resample_strategy: str = None,
        results_folder: str = "results"
) -> Tuple[pd.DataFrame, Metrics]:
    train_data = rebalance_data(train_data, seed = seed, strategy = resample_strategy)
    item_set_df = build_transactions(train_data, max_aa_length)
    print(train_data)

    print("a")
    sets = fpgrowth(item_set_df, min_support = min_support, use_colnames = True, max_len = max_antecedents, verbose = 0)
    print("b")
    rules = association_rules(sets, metric = "confidence", min_threshold = min_confidence)
    print("c")
    
    rules = rules[rules["consequents"].apply(lambda s: len(s) == 1 and next(iter(s)).startswith("VIRUS="))]
    print("d")
    rules = rules[rules["antecedents"].apply(len) >= min_antecedents]
    print("e")
    rules = rules[rules["lift"] >= min_lift].sort_values(["lift", "confidence", "support"], ascending = False)
    print("f")

    print(rules["consequents"].value_counts())

    rule_tuples = []
    for _, r in rules.iterrows():
        ante = set(r["antecedents"])
        cons_item = next(iter(r["consequents"]))
        label = cons_item.split("=", 1)[1]
        rule_tuples.append((
            ante, label, float(r["confidence"]), float(r["lift"]), float(r["support"]), len(ante)
        ))

    test_baskets = build_baskets(test_data, max_aa_length, include_virus = False)

    majority = train_data["virus_name"].mode().iloc[0]

    y_test = test_data["virus_name"].astype(str).tolist()
    y_pred = []

    for items in test_baskets:
        candidates = [rt for rt in rule_tuples if rt[0].issubset(items)]
        if candidates:
            candidates.sort(key = lambda t: (t[2], t[3], t[4], t[5]), reverse = True)
            y_pred.append(candidates[0][1])
        else:
            y_pred.append(majority)

    metrics = calculate_metrics(y_test, y_pred)

    params = {
        "max_aa_length": max_aa_length, 
        "min_confidence": min_confidence, 
        "min_lift": min_lift, 
        "min_support": min_support, 
        "min_antecedents": min_antecedents, 
        "max_antecedents": max_antecedents,
        "resample_strategy": resample_strategy
    }

    log_training(results_folder, "ar", pd.DataFrame(), params, metrics, seed)

    return rules, metrics

if __name__ == "__main__":
    seed = 561
    data_path = "proteins/data/virus_protein_repeats.csv"
    train_size = 0.8
    valid_size = 0
    test_size = 0.2

    spike_data, nucleocapsid_data = split_into_proteins(data_path)

    train_spike, valid_spike, test_spike = preprocess_data(spike_data, seed, train_size, valid_size, test_size)
    train_nucleocapsid, valid_nucleocapsid, test_nucleocapsid = preprocess_data(nucleocapsid_data, seed, train_size, valid_size, test_size)
    
    rules, metrics = train_ar_classifier(train_spike, test_spike, seed = 561, results_folder = "results/spike", resample_strategy = "under")
    rules, metrics = train_ar_classifier(train_nucleocapsid, test_nucleocapsid, seed = 561, results_folder = "results/nucleocapsid", resample_strategy = "under")
    