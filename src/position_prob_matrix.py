from typing import List, Sequence, Dict, Tuple, Generator

import numpy as np


class PositionProbabilityMatrix:
    def __init__(self,
                 seq_length: int,
                 data: Generator[str, None, None],
                 mapping: Sequence[str] = ("A", "C", "G", "T", "N")):
        self._mapping = {val: i for i, val in enumerate(mapping)}
        self._matrix = np.ones((len(mapping), seq_length), dtype=float)
        self._to_matrix(data)

    def _to_matrix(self, data: Generator[str, None, None]):
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
        if len(query_string) != self._matrix.shape[1]:
            raise ValueError(f"Query (current size={len(query_string)}) must be size {self._matrix.shape[1]}")
        total = 1.0
        for i, char in enumerate(query_string):
            total *= self._matrix[self._mapping[char], i]
        return total
