import unittest
from pathlib import Path

from src.models.annotation_model import AnnotationModel
from src.models.transition_model import TransitionModel


class TestTransitionModel(unittest.TestCase):
    def test_load(self):
        model = AnnotationModel(
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.mask.Repeats.fna"),
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.AbinitioGeneMark.gff3")
        )
        t_model = TransitionModel(annotation_models=[model])


if __name__ == '__main__':
    unittest.main()
