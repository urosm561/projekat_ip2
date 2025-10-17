import pandas as pd

from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth, association_rules
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
from matplotlib import pyplot as plt

from preprocess_data import preprocess_data, PositionalEncoder, PAD, Metrics, VIRUS_NAMES, calculate_metrics

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
        max_antecedents: int = 6
):
    item_set_df = build_transactions(train_data, max_aa_length)

    sets = fpgrowth(item_set_df, min_support = min_support, use_colnames = True, max_len = max_antecedents, verbose = 0)
    rules = association_rules(sets, metric = "confidence", min_threshold = min_confidence)

    rules = rules[rules["consequents"].apply(lambda s: len(s) == 1 and next(iter(s)).startswith("VIRUS="))]
    rules = rules[rules["antecedents"].apply(len) >= min_antecedents]
    rules = rules[rules["lift"] >= min_lift].sort_values(["lift", "confidence", "support"], ascending = False)

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

    return rules, metrics

if __name__ == "__main__":
    SEED = 561
    data_path = "proteins/data/virus_protein_repeats.csv"

    train_data, valid_data, test_data = preprocess_data(data = data_path, seed = SEED, train_size = 0.8, test_size = 0.2)
    print(train_data)
    rules, metrics = train_ar_classifier(train_data, test_data, seed = 561)

    print(rules)
    print(metrics)