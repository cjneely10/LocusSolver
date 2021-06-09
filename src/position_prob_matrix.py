from typing import List, Sequence, Dict

import numpy as np


class PPM:
    def __init__(self, data: List[str], mapping: Sequence[str] = ("A", "C", "G", "T", "N")):
        self._mapping = {val: i for i, val in enumerate(mapping)}
        self._matrix = np.ones((len(mapping), len(data[0])))
        self._to_matrix(data)

    def _to_matrix(self, data: List[str]):
        for sequence_string in data:
            for i, char in enumerate(sequence_string):
                self._matrix[self._mapping[char], i] += 1
        self._matrix /= np.sum(self._matrix, axis=0)

    @property
    def matrix(self) -> np.ndarray:
        return self._matrix

    @property
    def ids(self) -> Dict[str, int]:
        return self._mapping

    def score(self, query_string: str) -> float:
        total = 1.0
        for i, char in enumerate(query_string):
            if i > self.matrix.shape[1]:
                break
            total *= self._matrix[self.ids[char], i]
        return total
