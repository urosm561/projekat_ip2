import numpy as np
import pandas as pd
import torch.nn as nn

import torch
import os

from torch.utils.data import TensorDataset, DataLoader
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.model_selection import ParameterGrid
from typing import List, Tuple, Dict, Optional

from preprocess_data import encode, preprocess_data, PositionalEncoder, ALPHABET, Metrics, VIRUS_NAMES, calculate_metrics, rebalance_data, log_training, split_into_proteins

class SimpleMLP(nn.Module):
    def __init__(
            self, 
            in_dim: int, 
            num_classes: int, 
            n1: int = 256,
            n2: int = 128,
            n3: int = 64,
            p_drop: float = 0.1
    ):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, n1),
            nn.ReLU(),
            nn.Dropout(p_drop),
            nn.Linear(n1, n2),
            nn.ReLU(),
            nn.Dropout(p_drop),
            nn.Linear(n2, n3),
            nn.ReLU(),
            nn.Dropout(p_drop),
            nn.Linear(n3, num_classes)
        )

    def forward(
            self, 
            x: torch.Tensor
    ) -> torch.Tensor:
        return self.net(x)

def class_weights_from_train(
        y: np.ndarray, 
        VIRUS_NAMES: List[str]
) -> torch.Tensor:
    """Calculate the class weights, because the classes are heavily imbalanced."""
    counts = pd.Series(y).value_counts().reindex(VIRUS_NAMES).fillna(0).astype(int)
    inv = 1.0 / counts.replace(0, np.nan)
    inv = inv / inv.max()
    inv = inv.fillna(0.0)
    return torch.tensor(inv.values, dtype=torch.float32)

def build_features(
        train_data: str, 
        test_data: str
):
    """
    Prepare the data by positionally encoding amino acids up to the `max_aa_length`,
    and then one-hot encoding these columns. We get `max_aa_length` x 21 columns.
    """
    train_data = pd.read_csv(train_data)
    test_data = pd.read_csv(test_data)

    X_train = train_data.drop(columns = ["virus_name"])
    y_train = train_data["virus_name"]
    X_test = test_data.drop(columns = ["virus_name"])
    y_test = test_data["virus_name"]

    max_aa_length = X_train.shape[1] - 1
    print(X_train, max_aa_length)

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

    prep.fit(X_train, y_train)
    X_train = prep.transform(X_train)

    if not test_data.empty:
        X_test = prep.transform(X_test)
    else:
        X_test = np.empty((0, X_train.shape[1]), dtype = X_train.dtype)

    return X_train.astype(np.float32), X_test.astype(np.float32)

def encode_labels(
        y_train, 
        y_test
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Encode the labels in the simplest possible way.
    """
    class2idx = {c: i for i, c in enumerate(VIRUS_NAMES)}
    y_train = y_train.astype(str).map(class2idx).values
    y_test = y_test.astype(str).map(class2idx).values
    return y_train, y_test

def train_nn_classifier(
    train_data: str,
    test_data: str,
    seed: int,
    batch_size: int = 1024,
    lr: float = 1e-3,
    epochs: int = 20,
    p_drop: float = 0.1,
    n1: int = 256,
    n2: int = 128, 
    n3: int = 64,
    device: str = "cuda" if torch.cuda.is_available() else "cpu",
    patience = 5
) -> Tuple[nn.Module, Metrics, Metrics]:
    """
    Train the neural network and then evaluate it.

    Args:
        train_data (pd.DataFrame): Data to train the nn on.
        test_data (pd.DataFrame): Data to evaluate the nn on.
        seed (int): Control randomness.
        max_aa_length (int, optional): Maximum allowed number of amino acids in a repeat. Defaults to 7.
        batch_size (int, optional): Batch size for the DataLoaders. Defaults to 1024.
        lr (float, optional): Learning rate for the Adam optimizer. Defaults to 1e-3.
        epochs (int, optional): Number of epochs for training. Defaults to 20.
        p_drop (float, optional): Dropout parameter. Defaults to 0.1.
        device (str, optional): torch device.

    Returns:
        Tuple[nn.Module, Metrics]: Trained network and all the relevant metrics.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)
    train_df = pd.read_csv(train_data)
    test_df = pd.read_csv(test_data)

    y_train = train_df["virus_name"]
    y_test = test_df["virus_name"]
    
    X_train, X_test = build_features(train_data, test_data)
    y_train, y_test = encode_labels(y_train, y_test)

    in_dim = X_train.shape[1]
    num_classes = len(VIRUS_NAMES)

    Xtr = torch.from_numpy(X_train)
    ytr = torch.from_numpy(y_train).long()
    Xte = torch.from_numpy(X_test)
    yte = torch.from_numpy(y_test).long()

    train_loader = DataLoader(TensorDataset(Xtr, ytr), batch_size=batch_size, shuffle=True)
    test_loader  = DataLoader(TensorDataset(Xte, yte), batch_size=batch_size, shuffle=False)

    model = SimpleMLP(in_dim=in_dim, num_classes=num_classes, p_drop=p_drop, n1=n1, n2=n2, n3=n3).to(device)
    weights = class_weights_from_train(train_df["virus_name"].astype(str).values, VIRUS_NAMES).to(device)
    criterion = nn.CrossEntropyLoss(weight = weights)
    opt = torch.optim.AdamW(model.parameters(), lr = lr)

    best_val_loss = float("inf")
    epochs_since_improvement = 0

    metrics_test = None

    model.train()
    for epoch in range(epochs):
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            logits = model(xb)
            loss = criterion(logits, yb)

            opt.zero_grad()
            loss.backward()
            opt.step()

        model.eval()
        preds, trues = [], []
        val_loss = 0
        val_samples = 0

        if len(test_loader) > 0:

            with torch.no_grad():
                for xb, yb in test_loader:
                    xb = xb.to(device)
                    yb = yb.to(device)
                    logits = model(xb)
                    loss = criterion(logits, yb).sum()
                    val_loss += loss.item()

                    pred = logits.argmax(dim=1).cpu().numpy()
                    preds.append(pred)
                    trues.append(yb.numpy())
                    val_samples += yb.size(0)

            y_pred = np.concatenate(preds)
            y_test = np.concatenate(trues)

            y_pred = [VIRUS_NAMES[i] for i in y_pred.tolist()]
            y_test = [VIRUS_NAMES[i] for i in y_test.tolist()]

            metrics_test = calculate_metrics(y_test, y_pred)

            val_loss /= val_samples

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                epochs_since_improvement = 0
            else:
                epochs_since_improvement += 1

                if epochs_since_improvement == patience:
                    break

    model.eval()
    preds_train, trues_train = [], []
    with torch.no_grad():
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            logits = model(xb)
            preds_train.append(logits.argmax(dim=1).cpu().numpy())
            trues_train.append(yb.cpu().numpy())

    y_pred_train = np.concatenate(preds_train)
    y_true_train = np.concatenate(trues_train)
    y_pred_train = [VIRUS_NAMES[i] for i in y_pred_train.tolist()]
    y_true_train = [VIRUS_NAMES[i] for i in y_true_train.tolist()]
    metrics_train = calculate_metrics(y_true_train, y_pred_train)

    return model, metrics_test, metrics_train

def choose_best_nn(
        data_type, 
        seed, 
        max_aa_length, 
        resample_strategy,    
        p_drop,
        n1,
        n2, 
        n3,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
        batch_size: int = 1024,
        lr: float = 1e-3,
        epochs: int = 20,
        patience = 5,
        results_folder = "results"
):
    grid = ParameterGrid({
        "p_drop": p_drop,
        "n1": n1,
        "n2": n2,
        "n3": n3
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

                _, metrics, _ = train_nn_classifier(train_data, valid_data, seed, **params, 
                                                device = device, batch_size = batch_size, 
                                                lr = lr, epochs = epochs, patience = patience)
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

    _, metrics_test, metrics_train = train_nn_classifier(tmp_path, test_data, seed, **best_val_params, 
                                     device = device, batch_size = batch_size, 
                                     lr = lr, epochs = epochs, patience = patience)
    
    model, _, _ = train_nn_classifier(tmp_path2, empty_path, seed, **best_val_params, 
                                     device = device, batch_size = batch_size, 
                                     lr = lr, epochs = epochs, patience = patience)
    
    os.remove(tmp_path)
    os.remove(tmp_path2)
    os.remove(empty_path)

    best_val_params["max_aa_length"] = best_max_aa_length
    best_val_params["resample_strategy"] = best_resample_strategy
    log_training(results_folder, "nn", leaderboard, best_val_params, metrics_test, metrics_train, seed)

    return model, metrics, best_val_params

if __name__ == "__main__":
    data_path = "proteins/data/virus_protein_repeats.csv"
    seed = 561
    train_size = 0.7
    valid_size = 0.2
    test_size = 0.1

    model_nucleocapsid, metrics_nucleocapsid, params_nucleocapsid = choose_best_nn(
        "nucleocapsid", seed, 
        max_aa_length = [5, 6, 7],
        resample_strategy = ["over", None],
        p_drop = [0.1],
        n1 = [256, 128],
        n2 = [128, 64],
        n3 = [64, 32],
        results_folder = "results/nucleocapsid" 
    )

    model_spike, metrics_spike, params_spike = choose_best_nn(
        "spike", seed, 
        max_aa_length = [5, 6, 7],
        resample_strategy = ["over", None],
        p_drop = [0.1],
        n1 = [256, 128],
        n2 = [128, 64],
        n3 = [64, 32],
        results_folder = "results/spike" 
    )
