import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
from sklearn.preprocessing import OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer

from preprocess_data import preprocess_data, Metrics, PositionalEncoder, AMINO_ACIDS, PAD, CATS_21

def rf_classify(train_data: pd.DataFrame, test_data: pd.DataFrame, seed: int, max_aa_length: int = 7) -> Metrics:
    X_train = train_data[["repeat"]]
    y_train = train_data["virus_name"]
    X_test = test_data[["repeat"]]
    y_test = test_data["virus_name"]

    pos = PositionalEncoder(width = max_aa_length)
    ohe = OneHotEncoder(
        categories = [CATS_21] * max_aa_length,
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
        n_estimators = 200,
        max_depth = None,
        random_state = seed,
        criterion = "gini"
    )

    pipe = Pipeline([
        ("prep", prep),
        ("rf", rf_classifier)
    ])

    pipe.fit(X_train, y_train)
    y_pred = pipe.predict(X_test)

    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, average = "weighted")
    recall = recall_score(y_test, y_pred, average = "weighted")
    f1 = f1_score(y_test, y_pred, average = "weighted")
    cm = confusion_matrix(y_test, y_pred)

    return Metrics(accuracy, precision, recall, f1, cm)

if __name__ == "__main__":
    data_path = "proteins/data/virus_protein_repeats.csv"
    seed = 561
    train_size = 0.8
    valid_size = 0
    test_size = 0.2
    columns_to_drop = ["left_end", "right_end", "left_start", "right_start", "sequence_id", "length", "protein_name", "repeat_type"]
    train, valid, test = preprocess_data(data_path, seed, train_size, valid_size, test_size)

    metrics = rf_classify(train, test, seed)

    import seaborn as sns
    import matplotlib.pyplot as plt

    sns.heatmap(metrics.confusion_matrix, annot=True, fmt='d', cmap='inferno')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.show()
    