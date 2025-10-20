import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import ParameterGrid
from typing import Tuple, Optional

from preprocess_data import preprocess_data, Metrics, PositionalEncoder, ALPHABET, VIRUS_NAMES, calculate_metrics, rebalance_data, split_into_proteins, log_training

def train_rf_classifier(
        train_data: pd.DataFrame, 
        test_data: pd.DataFrame, 
        seed: int, 
        max_aa_length: int = 7,
        resample_strategy: Optional[str] = None,
        n_estimators: int = 200,
        criterion: str = "gini"
) -> Tuple[Pipeline, Metrics]:
    """Train the random forest classifier and then evaluate it.

    Args:
        train_data (pd.DataFrame): Data to train the rf on.
        test_data (pd.DataFrame): Data to evaluate the rf on.
        seed (int): Control randomness.
        max_aa_length (int, optional): Maximum allowed number of amino acids in a repeat. Defaults to 7.

    Returns:
        Tuple[Pipeline, Metrics]: Trained model and all the relevant metrics.
    """
    train_data = rebalance_data(train_data, seed = seed, strategy = resample_strategy)

    X_train = train_data[["repeat"]]
    y_train = train_data["virus_name"]
    X_test = test_data[["repeat"]]
    y_test = test_data["virus_name"]

    pos = PositionalEncoder(width = max_aa_length)
    ohe = OneHotEncoder(
        categories = [ALPHABET] * max_aa_length,
        handle_unknown = "error",
        sparse_output = False
    )

    prep = Pipeline([
        ("pos", pos),
        ("ohe", ColumnTransformer(
            transformers = [("onehot", ohe, [f"p{i + 1}" for i in range(max_aa_length)])],
            remainder = "drop"
        ))
    ])

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
    metrics = None

    if not test_data.empty:
        y_pred = pipe.predict(X_test)
        metrics = calculate_metrics(y_test, y_pred)

    return pipe, metrics

def choose_best_rf(
        train_data, 
        valid_data, 
        test_data, 
        seed, 
        max_aa_length, 
        resample_strategy, 
        n_estimators, 
        criterion, 
        results_folder = "results"
):
    grid = ParameterGrid({
        "max_aa_length": max_aa_length,
        "resample_strategy": resample_strategy,
        "n_estimators": n_estimators,
        "criterion": criterion
    })
    
    best_val_f1 = -1.0
    best_val_metrics = None
    best_val_params = None
    rows = []

    for params in grid:
        _, metrics = train_rf_classifier(train_data, valid_data, seed, **params)
        f1 = metrics.f1
        
        rows.append({
            **params,
            "accuracy": metrics.accuracy,
            "precision": metrics.precision,
            "recall": metrics.recall,
            "f1": metrics.f1
        })

        if f1 > best_val_f1:
            best_val_f1 = f1
            best_val_params = params
            best_val_metrics = metrics

    leaderboard = pd.DataFrame(rows).sort_values("f1", ascending = False).reset_index(drop = True)

    _, metrics = train_rf_classifier(pd.concat([train_data, valid_data], ignore_index = True), test_data, seed, **best_val_params)
    model, _ = train_rf_classifier(pd.concat([train_data, valid_data, test_data], ignore_index = True), pd.DataFrame(columns = ["repeat", "virus_name"]), seed, **best_val_params)

    log_training(results_folder, "rf", leaderboard, best_val_params, metrics, seed)

    return model, metrics, best_val_params

if __name__ == "__main__":
    data_path = "proteins/data/virus_protein_repeats.csv"
    seed = 561
    train_size = 0.7
    valid_size = 0.2
    test_size = 0.1

    spike_data, nucleocapsid_data = split_into_proteins(data_path)
    
    train_spike, valid_spike, test_spike = preprocess_data(spike_data, seed, train_size, valid_size, test_size)
    train_nucleocapsid, valid_nucleocapsid, test_nucleocapsid = preprocess_data(nucleocapsid_data, seed, train_size, valid_size, test_size)

    model_nucleocapsid, metrics_nucleocapsid, params_nucleocapsid = choose_best_rf(
        train_nucleocapsid, valid_nucleocapsid, test_nucleocapsid, seed, 
        max_aa_length = [5, 6, 7],
        resample_strategy = ["over", "under", None],
        n_estimators = [100, 200, 400],
        criterion = ["gini", "entropy"],
        results_folder = "results/nucleocapsid" 
    )

    model_spike, metrics_spike, params_spike = choose_best_rf(
        train_spike, valid_spike, test_spike, seed, 
        max_aa_length = [5, 6, 7],
        resample_strategy = ["over", "under", None],
        n_estimators = [100, 200, 400],
        criterion = ["gini", "entropy"],
        results_folder = "results/spike" 
    )
    