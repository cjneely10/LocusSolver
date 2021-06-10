from src.filters import FilterResult
from src.util.superlocus import SuperLocus


class Tier:
    def __init__(self, tier: int):
        self._tier = tier

    def __call__(self, super_locus: SuperLocus) -> FilterResult:
        pass
