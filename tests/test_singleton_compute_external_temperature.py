from acspype.processing import compute_external_temperature

TEST_DATA = 43210

def test():
    assert isinstance(TEST_DATA, int)    
    assert isinstance(compute_external_temperature(TEST_DATA), float)


