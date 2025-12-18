import numpy as np

def compute_base_uncertainty(accuracy: float, precision: float, other: float = 0.0) -> float:

    base_unc = np.sqrt(accuracy**2 + precision**2 + other**2)

    return base_unc



