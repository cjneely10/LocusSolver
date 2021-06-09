from typing import List, Tuple

from src.superlocus import SuperLocus, Feature


class SuperLocusList:
    def __init__(self, data: List[Tuple[str, Feature]]):
        self._data = self._sort(data)

    @property
    def sorted(self):
        return self._data

    @staticmethod
    def _sort(data: List[Tuple[str, Feature]]) -> List[SuperLocus]:
        data.sort(key=lambda feature_tuple: feature_tuple[1].start)
        stack = [SuperLocus(data[0][1], data[0][0])]
        for (identifier, feature) in data[1:]:
            top = stack[-1]
            if not top.overlaps(feature):
                stack.append(SuperLocus(feature, identifier))
            else:
                top.add_feature(feature, identifier)
        return stack
