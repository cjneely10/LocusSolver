from typing import Optional, Tuple, List

from src.util.feature import Feature

# name: List of applicable features as locus
# ex: genemark: [(1,50), (60, 100)]
FilterResult = List[Tuple[str, List[Feature]]]
