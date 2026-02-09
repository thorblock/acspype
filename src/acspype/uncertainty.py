import numpy as np

def compute_base_uncertainty(accuracy: float, precision: float, other: float = 0.0) -> float:
    """Compute the base uncertainty associated with a sensor."""
    base_unc = np.sqrt(accuracy**2 + precision**2 + other**2)
    return base_unc



