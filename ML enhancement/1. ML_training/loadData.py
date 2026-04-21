import numpy as np
from pathlib import Path
from tensorflow import keras

def load_trained_model(model_dir="trained_model"):
    model_dir = Path(model_dir)
    model = keras.models.load_model(model_dir / "mlp_256_128_64_32.keras")
    active_features = np.load(model_dir / "active_features.npy").tolist()
    return model, active_features

def predict_proba(model, X_raw_27, active_features, batch_size=32768):
    """
    X_raw_27: numpy array shape (N, 27) in the SAME feature order as training.
    returns probabilities shape (N,)
    """
    X = X_raw_27[:, active_features].astype(np.float32, copy=False)
    return model.predict(X, batch_size=batch_size, verbose=0).ravel()

# Example:
# model, active_features = load_trained_model()
# probs = predict_proba(model, X_raw_27, active_features)