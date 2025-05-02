from acspype.tscor import ACSTS4CorReader

TEST_FILEPATH = '../dev_tools/tscor/TS4.cor'

def test():
    tscor = ACSTS4CorReader(TEST_FILEPATH)

    assert len(tscor.wavelengths) == len(tscor.psi_t)
    assert len(tscor.wavelengths) == len(tscor.psi_s_c)
    assert len(tscor.wavelengths) == len(tscor.psi_s_a)

