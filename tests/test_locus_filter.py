import unittest
from pathlib import Path

from BCBio import GFF

from src.filters.simple_filters import *
from src.util.locus_filter import LocusFilter
from src.util.annotation import Annotation


class TestLocusFilter(unittest.TestCase):
    def test_load(self):
        model1 = Annotation(
            Path("data/sample.gmes.gff3"),
            "1"
        )
        model2 = Annotation(
            Path("data/sample.aug.gff3"),
            "2"
        )
        locus_filter = LocusFilter(annotation_models=[model1, model2])
        with open("data/sample.merged.gff3", "w") as merged_file:
            GFF.write(locus_filter.filter(Tier(1, {1: "1", 2: "2"})), merged_file)


if __name__ == '__main__':
    unittest.main()
