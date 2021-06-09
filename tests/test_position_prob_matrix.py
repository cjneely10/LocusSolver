import unittest

from src.position_prob_matrix import PositionProbabilityMatrix


class TestPPM(unittest.TestCase):
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

    def test_matrix(self):
        print(PositionProbabilityMatrix(TestPPM.test_list).score("GAGGTAAAC"))


if __name__ == '__main__':
    unittest.main()
