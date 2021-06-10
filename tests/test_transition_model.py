import unittest
from pathlib import Path

from src.models.annotation_model import AnnotationModel
from src.models.locus_filter import LocusFilter


class TestTransitionModel(unittest.TestCase):
    def test_load(self):
        model1 = AnnotationModel(
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.mask.Repeats.fna"),
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.AbinitioGeneMark.gff3"),
            "genemark"
        )
        model2 = AnnotationModel(
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.mask.Repeats.fna"),
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.augustus.gff3"),
            "augustus"
        )
        locus_filter = LocusFilter(annotation_models=[model1, model2])


if __name__ == '__main__':
    unittest.main()
