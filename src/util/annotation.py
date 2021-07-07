from pathlib import Path
from typing import Sequence, Optional, Tuple, List, Dict

from BCBio import GFF
from Bio.SeqFeature import SeqFeature

from src.util.superlocus import SuperLocus
from src.util.superlocus_list import SuperLocusList


class Annotation(dict):
    def __init__(self,
                 gff3_file: Path,
                 identifier: str,
                 features: Optional[Sequence[str]] = None,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        if features is None:
            features = ["transcript", "gene", "CDS", "exon"]
        self._gff_file = gff3_file
        self._features = features
        self._identifier = identifier

    def _to_features(self) -> List[Tuple[str, List[Tuple[str, SeqFeature]]]]:
        out = []
        with open(self._gff_file, "r") as gff_ptr:
            for record in GFF.parse(gff_ptr):
                out.append((record.id, [(self._identifier, feature) for feature in record.features]))
        return out

    @staticmethod
    def merge(models: List["Annotation"]) -> Dict[str, List[SuperLocus]]:
        feature_slls = {}
        for model in models:
            for (record_id, feature_list) in model._to_features():
                if record_id not in feature_slls.keys():
                    feature_slls[record_id] = feature_list
                else:
                    feature_slls[record_id].extend(feature_list)
        for record_id, feature_list in feature_slls.items():
            feature_slls[record_id] = SuperLocusList(feature_list).sorted
        return feature_slls
