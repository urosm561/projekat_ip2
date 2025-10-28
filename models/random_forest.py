import pandas as pd
import os

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import ParameterGrid
from typing import Tuple, Optional

from preprocess_data import encode, preprocess_data, Metrics, PositionalEncoder, ALPHABET, VIRUS_NAMES, calculate_metrics, rebalance_data, split_into_proteins, log_training

def train_rf_classifier(
        train_data: str, 
        test_data: str, 
        seed: int, 
        n_estimators: int = 200,
        criterion: str = "gini"
) -> Tuple[Pipeline, Metrics, Metrics]:
    """Train the random forest classifier and then evaluate it.

    Args:
        train_data (pd.DataFrame): Data to train the rf on.
        test_data (pd.DataFrame): Data to evaluate the rf on.
        seed (int): Control randomness.
        max_aa_length (int, optional): Maximum allowed number of amino acids in a repeat. Defaults to 7.

    Returns:
        Tuple[Pipeline, Metrics]: Trained model and all the relevant metrics.
    """
    train_data = pd.read_csv(train_data)
    test_data = pd.read_csv(test_data)

    X_train = train_data.drop(columns = ["virus_name"])
    y_train = train_data["virus_name"]
    X_test = test_data.drop(columns = ["virus_name"])
    y_test = test_data["virus_name"]

    max_aa_length = X_train.shape[1] - 1
    
    ohe_aa = OneHotEncoder(
        categories = [ALPHABET] * max_aa_length,
        handle_unknown = "error",
        sparse_output = False
    )
    ohe_i = OneHotEncoder(
        categories = [["in", "dn"]],
        handle_unknown = "error",
        sparse_output = False
    )

    prep = ColumnTransformer(
        transformers = [
            ("onehot_aa", ohe_aa, [f"p{i + 1}" for i in range(max_aa_length)]),
            ("onehot_i", ohe_i, ["repeat_type"])
        ],
        remainder = "drop"
    )

    rf_classifier = RandomForestClassifier(
        n_estimators = n_estimators,
        max_depth = None,
        random_state = seed,
        criterion = criterion,
        class_weight = "balanced_subsample",
        n_jobs = -1
    )

    pipe = Pipeline([
        ("prep", prep),
        ("rf", rf_classifier)
    ])

    pipe.fit(X_train, y_train)
    metrics_test = None
    metrics_train = None

    if not test_data.empty:
        y_pred = pipe.predict(X_test)
        metrics_test = calculate_metrics(y_test, y_pred)

    y_pred = pipe.predict(X_train)
    metrics_train = calculate_metrics(y_train, y_pred)

    return pipe, metrics_test, metrics_train

def choose_best_rf(
        data_type,
        seed, 
        max_aa_length, 
        resample_strategy, 
        n_estimators, 
        criterion, 
        results_folder = "results"
):
    grid = ParameterGrid({
        "n_estimators": n_estimators,
        "criterion": criterion
    })
    
    best_val_f1 = -1.0
    best_val_metrics = None
    best_val_params = None
    best_max_aa_length = None
    best_resample_strategy = None
    rows = []

    data_folder = "proteins/data"

    for l in max_aa_length:
        for r in resample_strategy:
            tmp_train = f"{data_type}_train_{r}_{seed}_{l}_final.csv"
            tmp_valid = f"{data_type}_valid_{seed}_{l}_final.csv"
            tmp_test = f"{data_type}_test_{seed}_{l}_final.csv"
            tmp_train = os.path.join(data_folder, tmp_train)
            tmp_valid = os.path.join(data_folder, tmp_valid)
            tmp_test = os.path.join(data_folder, tmp_test)

            if not os.path.exists(tmp_train) or not os.path.exists(tmp_valid) or not os.path.exists(tmp_test):
                main_path = os.path.join(data_folder, f"{data_type}.csv")
                tr, v, te = preprocess_data(main_path, seed = seed)
                tr1 = rebalance_data(os.path.join(data_folder, f"{data_type}_train_{seed}.csv"), seed = seed, strategy = "over")
                tr2 = rebalance_data(os.path.join(data_folder, f"{data_type}_train_{seed}.csv"), seed = seed, strategy = "under")
                tr3 = rebalance_data(os.path.join(data_folder, f"{data_type}_train_{seed}.csv"), seed = seed, strategy = None)
                tr1_path = os.path.join(data_folder, f"{data_type}_train_over_{seed}.csv")
                tr2_path = os.path.join(data_folder, f"{data_type}_train_under_{seed}.csv")
                tr3_path = os.path.join(data_folder, f"{data_type}_train_None_{seed}.csv")
                v_path = os.path.join(data_folder, f"{data_type}_valid_{seed}.csv")
                te_path = os.path.join(data_folder, f"{data_type}_test_{seed}.csv")
                paths = [tr1_path, tr2_path, tr3_path, v_path, te_path]
                for path in paths:
                    encode(path, max_aa_length = l)

    for l in max_aa_length:
        for r in resample_strategy:
            for params in grid:
                train_data = f"{data_type}_train_{r}_{seed}_{l}_final.csv"
                valid_data = f"{data_type}_valid_{seed}_{l}_final.csv"
                test_data = f"{data_type}_test_{seed}_{l}_final.csv"
                train_data = os.path.join(data_folder, train_data)
                valid_data = os.path.join(data_folder, valid_data)
                test_data = os.path.join(data_folder, test_data)

                _, metrics, _ = train_rf_classifier(train_data, valid_data, seed, **params)
                f1 = metrics.f1
                
                rows.append({
                    **params,
                    "max_aa_length": l,
                    "resample_strategy": r,
                    "accuracy": metrics.accuracy,
                    "precision": metrics.precision,
                    "recall": metrics.recall,
                    "f1": metrics.f1
                })

                if f1 > best_val_f1:
                    best_val_f1 = f1
                    best_val_params = params
                    best_max_aa_length = l
                    best_resample_strategy = r
                    best_val_metrics = metrics

    leaderboard = pd.DataFrame(rows).sort_values("f1", ascending = False).reset_index(drop = True)

    train_df = pd.read_csv(train_data)
    valid_df = pd.read_csv(valid_data)
    test_df = pd.read_csv(test_data)
    combined_train_valid = pd.concat([train_df, valid_df], ignore_index = True)
    tmp_path = os.path.join(data_folder, "combined_train_valid.csv")
    combined_train_valid.to_csv(tmp_path, index = False)
    combined_train_valid_test = pd.concat([train_df, valid_df, test_df], ignore_index = True)
    tmp_path2 = os.path.join(data_folder, "combined_train_valid_test.csv")
    combined_train_valid_test.to_csv(tmp_path2, index = False)
    empty_df = pd.DataFrame(columns = ["repeat", "repeat_type", "virus_name"])
    empty_path = os.path.join(data_folder, "empty.csv")
    empty_df.to_csv(empty_path)
    _, metrics_test, metrics_train = train_rf_classifier(tmp_path, test_data, seed, **best_val_params)
    model, _, _ = train_rf_classifier(tmp_path2, empty_path, seed, **best_val_params)
    os.remove(tmp_path)
    os.remove(tmp_path2)
    os.remove(empty_path)

    best_val_params["max_aa_length"] = best_max_aa_length
    best_val_params["resample_strategy"] = best_resample_strategy
    log_training(results_folder, "rf", leaderboard, best_val_params, metrics_test, metrics_train, seed)

    return model, metrics_test, best_val_params

if __name__ == "__main__":
    data_path = "proteins/data/virus_protein_repeats.csv"
    seed = 561
    train_size = 0.7
    valid_size = 0.2
    test_size = 0.1

    model_nucleocapsid, metrics_nucleocapsid, params_nucleocapsid = choose_best_rf(
        "nucleocapsid", seed, 
        max_aa_length = [5, 6, 7],
        resample_strategy = ["over", "under", None],
        n_estimators = [100, 200, 400],
        criterion = ["gini", "entropy"],
        results_folder = "results/nucleocapsid" 
    )

    model_spike, metrics_spike, params_spike = choose_best_rf(
        "spike", seed, 
        max_aa_length = [5, 6, 7],
        resample_strategy = ["over", "under", None],
        n_estimators = [100, 200, 400],
        criterion = ["gini", "entropy"],
        results_folder = "results/spike" 
    )
    