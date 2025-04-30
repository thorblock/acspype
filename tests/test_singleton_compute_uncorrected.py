import numpy as np
from acspype.processing import compute_uncorrected

TEST_DATA = {'a_signal': (1000,2000,3000,4000,5000),
             'a_reference': (100,200,300,400,500)}

def test():
    assert isinstance(TEST_DATA['a_signal'], tuple)
    assert isinstance(compute_uncorrected(TEST_DATA['a_signal'], TEST_DATA['a_reference']),np.ndarray)


