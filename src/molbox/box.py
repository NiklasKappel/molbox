from typing import Any

import numpy as np


class Box:
    def __init__(self, corner_lo: Any, corner_hi: Any) -> None:
        self.corner_lo = np.array(corner_lo)
        self.corner_hi = np.array(corner_hi)
