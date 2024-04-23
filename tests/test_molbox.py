import numpy as np
from molbox.box import Box
from molbox.molbox import MolBox


class TestMolBox:
    def test_init_side_effects(self):
        corner_lo = [0, 0, 0]
        corner_hi = [1, 1, 1]
        box = Box(corner_lo, corner_hi)
        molbox = MolBox(box=box)
        box.corner_lo = corner_hi
        assert np.allclose(molbox.box.corner_lo, corner_lo)
