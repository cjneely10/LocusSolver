import unittest
from pathlib import Path

from src.models.annotation_model import AnnotationModel
from src.models.transition_model import TransitionModel


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
        t_model = TransitionModel(annotation_models=[model1, model2])


if __name__ == '__main__':
    unittest.main()
