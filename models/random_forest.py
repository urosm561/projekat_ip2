import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
from sklearn.preprocessing import OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import ParameterGrid
from typing import Tuple, Optional
import seaborn as sns
import matplotlib.pyplot as plt

from preprocess_data import preprocess_data, Metrics, PositionalEncoder, ALPHABET, VIRUS_NAMES, calculate_metrics, rebalance_data

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

def choose_best_rf(train_data, valid_data, test_data, seed, max_aa_length, resample_strategy, n_estimators, criterion):
    grid = ParameterGrid({
        "max_aa_length": max_aa_length,
        "resample_strategy": resample_strategy,
        "n_estimators": n_estimators,
        "criterion": criterion
    })
    
    best_val_f1 = 0

    for params in grid:
        model, metrics = train_rf_classifier(train_data, valid_data, seed, **params)
        f1 = metrics.f1

        if f1 > best_val_f1:
            best_val_f1 = f1
            best_val_params = params

    _, metrics = train_rf_classifier(pd.concat([train_data, valid_data], ignore_index = True), test_data, seed, **best_val_params)
    model, _ = train_rf_classifier(pd.concat([train_data, valid_data, test_data], ignore_index = True), pd.DataFrame(columns = ["repeat", "virus_name"]), seed, **best_val_params)

    return model, metrics, best_val_params

if __name__ == "__main__":
    data_path = "proteins/data/virus_protein_repeats.csv"
    seed = 561
    train_size = 0.7
    valid_size = 0.2
    test_size = 0.1
    columns_to_drop = ["left_end", "right_end", "left_start", "right_start", "sequence_id", "length", "protein_name", "repeat_type"]
    train, valid, test = preprocess_data(data_path, seed, train_size, valid_size, test_size)

    model, metrics, params = choose_best_rf(
        train, valid, test, seed, 
        max_aa_length = [5],
        resample_strategy = ["over"],
        n_estimators = [100, 200],
        criterion = ["gini", "entropy"]    
    )

    print(metrics)
    print(params)

    sns.heatmap(metrics.confusion_matrix, annot=True, fmt='d', cmap='inferno', 
                xticklabels = VIRUS_NAMES, yticklabels = VIRUS_NAMES)
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.show()
    