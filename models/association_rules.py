import pandas as pd

from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth, association_rules
from typing import Tuple

from preprocess_data import preprocess_data, PositionalEncoder, Metrics, PAD, calculate_metrics, split_into_proteins, rebalance_data, log_training

def ngrams(s: str, n: int):
    return {s[i : i + n] for i in range(len(s) - n + 1)} if len(s) >= n else set()

def build_baskets(data: str, include_virus: bool = True):
    # pe = PositionalEncoder(width = max_aa_length)
    # pos_data = pe.transform(data[["repeat"]])
    data = pd.read_csv(data)
    max_aa_length = data.shape[1] - 2

    baskets = []
    for id, row in data.iterrows():
        pos_items = [f"P{i + 1}={row[f'p{i + 1}']}" for i in range(max_aa_length) if row[f"p{i + 1}"] != PAD]
        raw = "".join([row[f"p{i + 1}"] for i in range(max_aa_length) if row[f"p{i + 1}"] != PAD])
        true_len = len(raw)
        items = pos_items[:]

        if true_len <= 3: items.append("LEN=<=3")
        elif true_len <=5: items.append("LEN=4-5")
        else: items.append("LEN=6-7")

        for bg in ngrams(raw, 2): items.append(f"BG={bg}")
        # for tg in ngrams(raw, 3): items.append(f"TG={tg}")

        for ch in set(raw):
            items.append(f"AA={ch}")

        items.append(f"RT={row['repeat_type']}")

        if include_virus:
            items.append(f"VIRUS={row['virus_name']}")

        baskets.append(sorted(set(items)))

    print(baskets)

    # for (idx, row), (_, orig) in zip(pos_data.iterrows(), data.iterrows()):
    #     pos_items = [f"P{i+1}={row[f'p{i+1}']}" for i in range(max_aa_length) if row[f"p{i+1}"] != PAD]
    #     raw = str(orig["repeat"])[-max_aa_length:]
    #     true_len = len(raw)

    #     items = pos_items[:]

    #     if true_len <= 3: items.append("LEN=<=3")
    #     elif true_len <=5: items.append("LEN=4-5")
    #     else: items.append("LEN=6-7")

    #     for bg in ngrams(raw, 2): items.append(f"BG={bg}")
    #     # for tg in ngrams(raw, 3): items.append(f"TG={tg}")

    #     for ch in set(raw):
    #         items.append(f"AA={ch}")

    #     if include_virus:
    #         items.append(f"VIRUS={orig['virus_name']}")

    #     baskets.append(sorted(set(items)))

    return baskets

def build_transactions(data: str, include_virus: bool = True):
    baskets = build_baskets(data, include_virus)
    
    te = TransactionEncoder()
    item_set = te.fit_transform(baskets)
    item_set_df = pd.DataFrame(data = item_set, columns = te.columns_)

    return item_set_df

def train_ar_classifier(
        train_data: str, 
        test_data: str,
        min_confidence: float = 0.7, 
        min_lift: float = 1.3, 
        min_support: float = 0.001, 
        min_antecedents: int = 2, 
        max_antecedents: int = 6,
        results_folder: str = "results"
) -> Tuple[pd.DataFrame, Metrics, Metrics]:
    resample_strategy = train_data.split("_")[-4]
    seed = train_data.split("_")[-3]
    item_set_df = build_transactions(train_data, include_virus = True)
    train_df = pd.read_csv(train_data)
    test_df = pd.read_csv(test_data)
    max_aa_length = train_df.shape[1] - 2

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

    rule_tuples = []
    for _, r in rules.iterrows():
        ante = set(r["antecedents"])
        cons_item = next(iter(r["consequents"]))
        label = cons_item.split("=", 1)[1]
        rule_tuples.append((
            ante, label, float(r["confidence"]), float(r["lift"]), float(r["support"]), len(ante)
        ))

    test_baskets = build_baskets(test_data, include_virus = False)

    majority = train_df["virus_name"].mode().iloc[0]

    y_test = test_df["virus_name"].astype(str).tolist()
    y_pred = []

    for items in test_baskets:
        candidates = [rt for rt in rule_tuples if rt[0].issubset(items)]
        if candidates:
            candidates.sort(key = lambda t: (t[2], t[3], t[4], t[5]), reverse = True)
            y_pred.append(candidates[0][1])
        else:
            y_pred.append(majority)

    metrics_test = calculate_metrics(y_test, y_pred)

    train_baskets = build_baskets(train_data, include_virus = False)

    y_train = train_df["virus_name"].astype(str).tolist()
    y_pred = []

    for items in train_baskets:
        candidates = [rt for rt in rule_tuples if rt[0].issubset(items)]
        if candidates:
            candidates.sort(key = lambda t: (t[2], t[3], t[4], t[5]), reverse = True)
            y_pred.append(candidates[0][1])
        else:
            y_pred.append(majority)

    metrics_train = calculate_metrics(y_train, y_pred)

    params = {
        "min_confidence": min_confidence, 
        "min_lift": min_lift, 
        "min_support": min_support, 
        "min_antecedents": min_antecedents, 
        "max_antecedents": max_antecedents,
        "max_aa_length": max_aa_length, 
        "resample_strategy": resample_strategy
    }

    log_training(results_folder, "ar", pd.DataFrame(), params, metrics_test, metrics_train, seed)

    return rules, metrics_test, metrics_train

if __name__ == "__main__":
    seed = 561
    data_path = "proteins/data/virus_protein_repeats.csv"
    train_size = 0.9
    valid_size = 0
    test_size = 0.1

    valid_under_spike = pd.read_csv("proteins/data/spike_valid_561_7_final.csv")
    valid_under_nucleocapsid = pd.read_csv("proteins/data/nucleocapsid_valid_561_7_final.csv")

    from imblearn.under_sampling import RandomUnderSampler

    Xs, ys = valid_under_spike.drop(columns = ["virus_name"]), valid_under_spike["virus_name"]
    Xn, yn = valid_under_nucleocapsid.drop(columns = ["virus_name"]), valid_under_nucleocapsid["virus_name"]

    rus = RandomUnderSampler(random_state = seed)
    Xs, ys = rus.fit_resample(Xs, ys)
    rus = RandomUnderSampler(random_state = seed)
    Xn, yn = rus.fit_resample(Xn, yn)

    valid_under_spike = pd.concat([Xs, ys], axis = 1)
    valid_under_nucleocapsid = pd.concat([Xn, yn], axis = 1)

    valid_under_spike.to_csv(f"proteins/data/spike_valid_under_{seed}_7_final.csv", index = False)
    valid_under_nucleocapsid.to_csv(f"proteins/data/nucleocapsid_valid_under_{seed}_7_final.csv", index = False)

    pd.concat([
        pd.read_csv(f"proteins/data/nucleocapsid_valid_under_{seed}_7_final.csv"),
        pd.read_csv(f"proteins/data/nucleocapsid_train_under_{seed}_7_final.csv")
    ], ignore_index = True).to_csv(f"proteins/data/nucleocapsid_train_valid_under_{seed}_7_final.csv", index = False)
    pd.concat([
        pd.read_csv(f"proteins/data/spike_valid_under_{seed}_7_final.csv"),
        pd.read_csv(f"proteins/data/spike_train_under_{seed}_7_final.csv")
    ], ignore_index = True).to_csv(f"proteins/data/spike_train_valid_under_{seed}_7_final.csv", index = False)
    
    _, _, _ = train_ar_classifier(f"proteins/data/spike_train_valid_under_{seed}_7_final.csv", f"proteins/data/spike_test_{seed}_7_final.csv", results_folder = "results/spike")
    _, _, _ = train_ar_classifier(f"proteins/data/nucleocapsid_train_valid_under_{seed}_7_final.csv", f"proteins/data/nucleocapsid_test_{seed}_7_final.csv", results_folder = "results/nucleocapsid")
    