import unittest

from src.superlocus import Feature
from src.superlocus_tree import SuperLocusList


class TestSuperLocusList(unittest.TestCase):
    def test_insert(self):
        slt = SuperLocusList([
            ("test", Feature(50, 100, -1)),
            ("test", Feature(1, 40, -1)),
            ("test", Feature(150, 200, -1)),
            ("test", Feature(90, 160, -1)),
        ])
        for super_locus in slt.sorted:
            print(super_locus)

    def test_sort(self):
        self.fail()


if __name__ == '__main__':
    unittest.main()
