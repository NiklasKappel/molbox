import numpy as np
from molbox.box import Box


class TestBox:
    def test_init(self):
        corner_lo = [0, 0, 0]
        corner_hi = [1, 1, 1]
        box = Box(corner_lo, corner_hi)
        assert np.allclose(box.corner_lo, corner_lo)
        assert np.allclose(box.corner_hi, corner_hi)
