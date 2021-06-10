import unittest
from pathlib import Path

from src.filters.simple_filters import longest_cds
from src.util.annotation import Annotation
from src.locus_filter import LocusFilter


class TestTransitionModel(unittest.TestCase):
    def test_load(self):
        model1 = Annotation(
            Path("data/sample.fna"),
            Path("data/sample.gmes.gff3"),
            "genemark"
        )
        model2 = Annotation(
            Path("data/sample.fna"),
            Path("data/sample.aug.gff3"),
            "augustus"
        )
        locus_filter = LocusFilter(annotation_models=[model1, model2])
        print(next(locus_filter.filter(longest_cds)))


if __name__ == '__main__':
    unittest.main()
