import unittest

from src.models.position_prob_matrix import PositionProbabilityMatrix


class TestPPM(unittest.TestCase):
    @staticmethod
    def generator_fxn():
        test_list = [
            "GAGGTAAAC",
            "TCCGTAAGT",
            "CAGGTTGGA",
            "ACAGTCAGT",
            "TAGGTCATT",
            "TAGGTACTG",
            "ATGGTAACT",
            "CAGGTATAC",
            "TGTGTGAGT",
            "AAGGTAAGT",
        ]
        for val in test_list:
            yield val

    def test_matrix(self):
        print(PositionProbabilityMatrix(9, TestPPM.generator_fxn()).score("GAGGTAAAC"))


if __name__ == '__main__':
    unittest.main()
