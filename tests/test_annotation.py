import unittest
from pathlib import Path

from src.util.annotation import Annotation


class TestAnnotation(unittest.TestCase):
    def test_load(self):
        model = Annotation(
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.mask.Repeats.fna"),
            Path("data/GCA_000151265.1_Micromonas_pusilla_CCMP1545_v2.0_genomic.AbinitioGeneMark.gff3")
        ).intron_model
        print(model.matrix)


if __name__ == '__main__':
    unittest.main()
