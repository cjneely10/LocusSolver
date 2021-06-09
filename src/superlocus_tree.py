from src.superlocus import SuperLocus, Feature


class SuperLocusTree:
    class _TreeNode:
        def __init__(self, super_locus: SuperLocus):
            self.current = super_locus
            self.left = None
            self.right = None

    def __init__(self):
        self.head = None

    def insert(self, feature: Feature):
        pass

    def sort(self):
        pass
