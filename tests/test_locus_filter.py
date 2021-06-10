import unittest
from pathlib import Path

from src.util.annotation import Annotation
from src.locus_filter import LocusFilter


class TestTransitionModel(unittest.TestCase):
    def test_load(self):
        model1 = Annotation(
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.mask.Repeats.fna"),
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.AbinitioGeneMark.gff3"),
            "genemark"
        )
        model2 = Annotation(
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.mask.Repeats.fna"),
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.augustus.gff3"),
            "augustus"
        )
        locus_filter = LocusFilter(annotation_models=[model1, model2])
        print(next(locus_filter.filter()))


if __name__ == '__main__':
    unittest.main()
