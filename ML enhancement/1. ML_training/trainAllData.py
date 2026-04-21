# -*- coding: utf-8 -*-
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional, Dict

import numpy as np
import tensorflow as tf
tf.get_logger().setLevel("ERROR")

from tensorflow import keras
from tensorflow.keras import layers

from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_curve, average_precision_score

# -------------------------
# Features
# -------------------------
kFeatureNames = [
    "is_found_by_MST","is_found_by_RC","s_frac","dist_in_avg_norm","dist_2min_out_avg_norm",
    "dist_2min_in_avg_norm","rc_in_avg_norm","rc_2min_out_avg_norm","rc_2min_in_avg_norm",
    "lp_in_avg","lp_2max_out_avg","lp_2max_in_avg","dist_knn_out_frac_avg","lp_knn_out_frac_avg",
    "rc_knn_out_frac_avg","dist_first_in_avg_norm","lp_first_in_avg","rc_first_in_avg_norm",
    "dist_boundary_gap_avg_norm","lp_boundary_gap_avg","rc_boundary_gap_avg_norm","mst_dist_avg_norm",
    "mst_1mx_avg","mst_rc_avg_norm","greedy_cycle_dist_avg_norm","greedy_cycle_1mx_avg",
    "greedy_cycle_rc_avg_norm",
]
active_features = list(range(27))

# -------------------------
# Data IO
# -------------------------
@dataclass
class Sample:
    x: List[float]
    label: int
    verts: str = ""
    id: int = -1

def read_file(path: Path) -> List[Sample]:
    out: List[Sample] = []
    kNumFeatures = 27
    kFirstFeatureCol = 3
    kLabelCol = kFirstFeatureCol + kNumFeatures  # 30
    kMinCols = kLabelCol + 1                     # 31

    with path.open("r", newline="") as f:
        reader = csv.reader(f)
        first = True
        for row in reader:
            if not row:
                continue
            if first:
                first = False
                continue
            if len(row) < kMinCols:
                continue

            try:
                sid = int(row[0])
                verts = row[2]
                feats = [float(row[kFirstFeatureCol + i]) for i in range(kNumFeatures)]
                label = int(row[kLabelCol])
            except ValueError:
                continue

            out.append(Sample(x=feats, label=label, verts=verts, id=sid))
    return out

def load_folder(folder: str, id_offset: int = 0) -> List[Sample]:
    folder_path = Path(folder)
    all_samples: List[Sample] = []
    file_count = 0

    for p in folder_path.iterdir():
        if not p.is_file():
            continue
        if p.suffix.lower() not in (".txt", ".csv"):
            continue

        tmp = read_file(p)
        if tmp:
            file_count += 1
            if id_offset != 0:
                for s in tmp:
                    s.id += id_offset
            all_samples.extend(tmp)

    print(f"Finished reading files: {file_count} | total samples = {len(all_samples)}")
    return all_samples

def to_numpy(data: List[Sample], feat_idx: List[int]) -> Tuple[np.ndarray, np.ndarray]:
    X = np.array([[s.x[j] for j in feat_idx] for s in data], dtype=np.float32)
    y = np.array([s.label for s in data], dtype=np.int32)
    return X, y

# -------------------------
# Model
# -------------------------
def make_mlp_model(in_dim: int, hidden: List[int]) -> keras.Model:
    inputs = keras.Input(shape=(in_dim,))
    x = inputs
    for h in hidden:
        x = layers.Dense(h, activation="tanh")(x)
    outputs = layers.Dense(1, activation="sigmoid")(x)
    return keras.Model(inputs, outputs, name="MLPNet")

# -------------------------
# PR / threshold helpers
# -------------------------
def _sorted_pr(precision: np.ndarray, recall: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    precision = np.asarray(precision, dtype=np.float64)
    recall = np.asarray(recall, dtype=np.float64)
    order = np.argsort(recall)  # ascending recall
    return precision[order], recall[order]

def _precision_envelope(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=np.float64)
    return np.maximum.accumulate(p[::-1])[::-1]  # non-increasing w.r.t recall

def precision_at_recall(precision: np.ndarray, recall: np.ndarray, r0: float) -> float:
    p, r = _sorted_pr(precision, recall)
    p = _precision_envelope(p)
    return float(np.interp(float(r0), r, p))

def recall_at_precision(precision: np.ndarray, recall: np.ndarray, p0: float) -> float:
    p, r = _sorted_pr(precision, recall)
    p = _precision_envelope(p)
    mask = p >= float(p0)
    return float(r[mask].max()) if np.any(mask) else 0.0

def thresholds_for_targets(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    recall_targets=(0.99, 0.999),
    precision_targets=(1/7, 1/9),
) -> Dict[str, Optional[float]]:
    """
    Finds thresholds based on direct threshold scanning.
    - For recall targets: returns the HIGHEST threshold whose recall >= target
      (i.e., most selective while still meeting the recall requirement).
    - For precision targets: returns the HIGHEST threshold whose precision >= target!!
    """
    y_true = np.asarray(y_true, dtype=np.int32)
    y_prob = np.asarray(y_prob, dtype=np.float64)

    # Candidate thresholds: use unique probs (plus endpoints)
    thr = np.unique(y_prob)
    if thr.size == 0:
        return {}

    # Evaluate from high->low threshold (more selective to less selective)
    thr = thr[::-1]

    P = (y_true == 1).sum()
    out: Dict[str, Optional[float]] = {}

    # Pre-sort probs descending for faster cumulative computation
    order = np.argsort(-y_prob)
    y_sorted = y_true[order]
    prob_sorted = y_prob[order]

    # As threshold decreases, predicted positives accumulate
    tp_cum = np.cumsum(y_sorted == 1)
    fp_cum = np.cumsum(y_sorted == 0)

    # For any threshold t, predicted positives are those with prob >= t.
    # We can map each unique threshold to the last index where prob >= t.
    # Using searchsorted on descending arrays:
    # Convert to ascending for searchsorted convenience:
    prob_asc = prob_sorted[::-1]
    tp_asc = tp_cum[::-1]
    fp_asc = fp_cum[::-1]

    def metrics_at_threshold(t: float) -> Tuple[float, float]:
        # index of first element in ascending prob_asc that is >= t
        idx = np.searchsorted(prob_asc, t, side="left")
        if idx >= prob_asc.size:
            # no predicted positives
            tp = 0
            fp = 0
        else:
            # predicted positives are indices idx..end in ascending array
            tp = int(tp_asc[idx])
            fp = int(fp_asc[idx])

        recall = (tp / P) if P > 0 else 0.0
        precision = (tp / (tp + fp)) if (tp + fp) > 0 else 1.0
        return precision, recall

    # Recall targets
    for r0 in recall_targets:
        best_t = None
        for t in thr:
            p, r = metrics_at_threshold(float(t))
            if r >= float(r0):
                best_t = float(t)
                break  # because thr is high->low: first match is highest t
        out[f"threshold_recall_{str(r0).replace('.','p')}"] = best_t

    # Precision targets
    for p0 in precision_targets:
        best_t = None
        for t in thr:
            p, r = metrics_at_threshold(float(t))
            if p >= float(p0):
                best_t = float(t)
                break
        out[f"threshold_precision_{str(p0).replace('.','p')}"] = best_t

    return out

# -------------------------
# Train + Evaluate + Save
# -------------------------
def main():
    # Folders
    F_SMALL = r"/data/storage_1_5_tb/TSP/FEATURES_RUN/PORTGEN/FEATURES/PORTGEN_SMALL"
    F_MID   = r"/data/storage_1_5_tb/TSP/FEATURES_RUN/PORTGEN/FEATURES/PORTGEN_MID"
    F_PORTC = r"/data/storage_1_5_tb/TSP/FEATURES_RUN/PORTCGEN/FEATURES/PORTCGEN"

    # Reproducibility
    seed = 12345
    tf.random.set_seed(seed)
    np.random.seed(seed)

    # Load + merge
    samples_small = load_folder(F_SMALL, id_offset=0)
    samples_mid   = load_folder(F_MID,   id_offset=1_000_000_000)
    samples_portc = load_folder(F_PORTC, id_offset=2_000_000_000)
    all_samples = samples_small + samples_mid + samples_portc

    print(
        f"Merged dataset: small={len(samples_small)} | "
        f"mid={len(samples_mid)} | portc={len(samples_portc)} | "
        f"total={len(all_samples)}"
    )
    if not all_samples:
        raise SystemExit("No samples loaded from any folder.")

    # Convert to arrays
    X_all, y_all = to_numpy(all_samples, active_features)

    # Shuffle + stratified 80/20 split
    Xtr, Xte, ytr, yte = train_test_split(
        X_all, y_all,
        test_size=0.20,
        random_state=seed,
        shuffle=True,
        stratify=y_all,
    )

    print(f"Split: train={len(ytr)} test={len(yte)} | pos(train)={(ytr==1).sum()} pos(test)={(yte==1).sum()}")

    in_dim = X_all.shape[1]

    # Architecture
    hidden = [256, 128, 64, 32]
    model = make_mlp_model(in_dim, hidden)

    # Hyperparams
    lr = 1e-3
    epochs = 10
    batch_size = 32768

    # class_weight on TRAIN only (avoid leakage)
    n0 = int((ytr == 0).sum())
    n1 = int((ytr == 1).sum())
    pos_weight_val = (n0 / n1) if n1 > 0 else 1.0
    class_weight = {0: 1.0, 1: float(pos_weight_val)}
    print(f"class_weight(train) = {class_weight} (n0={n0}, n1={n1})")

    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=lr),
        loss=keras.losses.BinaryCrossentropy(from_logits=False),
    )

    model.fit(
        Xtr, ytr,
        epochs=epochs,
        batch_size=batch_size,
        class_weight=class_weight,
        verbose=2,
        shuffle=True,
    )

    # Evaluate on test
    y_prob = model.predict(Xte, batch_size=batch_size, verbose=0).ravel().astype(np.float64)

    ap = float(average_precision_score(yte, y_prob))
    precision, recall, thresholds = precision_recall_curve(yte, y_prob)

    p_at_r099  = precision_at_recall(precision, recall, 0.99)
    p_at_r0999 = precision_at_recall(precision, recall, 0.999)

    r_at_p1_7 = recall_at_precision(precision, recall, 1/7)
    r_at_p1_9 = recall_at_precision(precision, recall, 1/9)

    thr_dict = thresholds_for_targets(
        y_true=yte,
        y_prob=y_prob,
        recall_targets=(0.99, 0.999),
        precision_targets=(1/7, 1/9),
    )

    # Print metrics
    print("\n================ TEST METRICS ================")
    print(f"AP (Average Precision / full AUPRC): {ap:.6f}")
    print(f"Precision @ Recall=0.99:            {p_at_r099:.6f}")
    print(f"Precision @ Recall=0.999:           {p_at_r0999:.6f}")
    print(f"Recall @ Precision=1/7:             {r_at_p1_7:.6f}")
    print(f"Recall @ Precision=1/9:             {r_at_p1_9:.6f}")

    print("\n================ THRESHOLDS (on TEST) ================")
    print(f"Highest threshold with Recall >= 0.99:   {thr_dict.get('threshold_recall_0p99')}")
    print(f"Highest threshold with Recall >= 0.999:  {thr_dict.get('threshold_recall_0p999')}")
    print(f"Highest threshold with Precision >= 1/7: {thr_dict.get('threshold_precision_0p1428')}")
    print(f"Highest threshold with Precision >= 1/9: {thr_dict.get('threshold_precision_0p1111')}")

    # Save model + metadata + PR arrays
    out_dir = Path("trained_model")
    out_dir.mkdir(parents=True, exist_ok=True)

    model_path = out_dir / "mlp_256_128_64_32.keras"
    model.save(model_path)

    np.save(out_dir / "active_features.npy", np.array(active_features, dtype=np.int32))
    np.save(out_dir / "feature_names.npy", np.array(kFeatureNames, dtype=object))

    # Save PR curve arrays:
    # sklearn: precision & recall length = len(thresholds)+1
    # We'll also store aligned per-threshold precision/recall (dropping the first point)
    pr_path = out_dir / "test_pr_curve_data.npz"
    np.savez_compressed(
        pr_path,
        y_test=yte.astype(np.int32),
        y_prob=y_prob.astype(np.float64),
        precision=precision.astype(np.float64),
        recall=recall.astype(np.float64),
        thresholds=thresholds.astype(np.float64),
        # aligned with thresholds (each threshold corresponds to precision[i+1], recall[i+1])
        precision_at_threshold=precision[1:].astype(np.float64),
        recall_at_threshold=recall[1:].astype(np.float64),
        ap=np.array([ap], dtype=np.float64),
        p_at_r099=np.array([p_at_r099], dtype=np.float64),
        p_at_r0999=np.array([p_at_r0999], dtype=np.float64),
        r_at_p1_7=np.array([r_at_p1_7], dtype=np.float64),
        r_at_p1_9=np.array([r_at_p1_9], dtype=np.float64),
        threshold_recall_0p99=np.array([thr_dict.get("threshold_recall_0p99")], dtype=object),
        threshold_recall_0p999=np.array([thr_dict.get("threshold_recall_0p999")], dtype=object),
        threshold_precision_1_over_7=np.array([thr_dict.get("threshold_precision_0p1428")], dtype=object),
        threshold_precision_1_over_9=np.array([thr_dict.get("threshold_precision_0p1111")], dtype=object),
    )

    print("\n================ SAVED ================")
    print(f"Model:         {model_path.resolve()}")
    print(f"PR data (npz): {pr_path.resolve()}")
    print(f"Active feats:  {(out_dir / 'active_features.npy').resolve()}")

if __name__ == "__main__":
    main()