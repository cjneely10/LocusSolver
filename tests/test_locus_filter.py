import unittest
from pathlib import Path

from BCBio import GFF

from src.filters.simple_filters import *
from src.locus_filter import LocusFilter
from src.util.annotation import Annotation


class TestLocusFilter(unittest.TestCase):
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
        with open("data/sample.merged.gff3", "w") as merged_file:
            GFF.write(locus_filter.filter(PriorityFilter({1: "genemark", 2: "augustus"})), merged_file)


if __name__ == '__main__':
    unittest.main()
