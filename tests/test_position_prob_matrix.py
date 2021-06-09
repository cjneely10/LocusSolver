import unittest

from src.position_prob_matrix import PPM


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
        print(PPM(TestPPM.test_list).matrix)


if __name__ == '__main__':
    unittest.main()
