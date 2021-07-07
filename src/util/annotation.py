from pathlib import Path
from typing import Sequence, Optional, Tuple, List, Dict

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

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
        self._load()

    @staticmethod
    def load_genome(genome_file: Path):
        genome_ptr = open(genome_file, "r")
        Annotation.genome_dict = SeqIO.to_dict(SeqIO.parse(genome_ptr, "fasta"))
        genome_ptr.close()

    def _to_features(self) -> List[Tuple[str, List[Tuple[str, SeqFeature]]]]:
        out = []
        gene_list: List[SeqFeature]
        for record_id, gene_list in self.items():
            _add = []
            for feature in gene_list:
                _add.append((self._identifier, feature))
            out.append((record_id, _add))
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

    def _load(self):
        gff3_ptr = open(self._gff_file, "r")
        record: SeqRecord
        for record in GFF.parse(gff3_ptr):
            self[record.id] = record.features
        gff3_ptr.close()
