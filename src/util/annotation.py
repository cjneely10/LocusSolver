from pathlib import Path
from typing import Sequence, Optional, Generator, Tuple, List, Dict

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from src.models.position_prob_matrix import PositionProbabilityMatrix
from src.util.feature import Feature
from src.util.superlocus import SuperLocus
from src.util.superlocus_list import SuperLocusList


class Annotation(dict):
    genome_dict: dict = None

    def __init__(self,
                 genome_file: Path,
                 gff3_file: Path,
                 identifier: str,
                 features: Optional[Sequence[str]] = None,
                 front_size: int = 10,
                 end_size: int = 20,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        if features is None:
            features = ["transcript", "gene", "CDS"]
        self._gff_file = gff3_file
        self._features = features
        self._identifier = identifier
        if Annotation.genome_dict is None:
            Annotation.load_genome(genome_file)
        self._intron_model = PositionProbabilityMatrix(front_size + end_size, self._intron_list(front_size, end_size))

    @property
    def intron_model(self) -> PositionProbabilityMatrix:
        return self._intron_model

    @staticmethod
    def load_genome(genome_file: Path):
        genome_ptr = open(genome_file, "r")
        Annotation.genome_dict = SeqIO.to_dict(SeqIO.parse(genome_ptr, "fasta"))
        genome_ptr.close()

    def _to_features(self) -> List[Tuple[str, List[Tuple[str, Feature]]]]:
        out = []
        for record_id, gene_dict in self.items():
            _add = []
            for feature in gene_dict:
                _add.append((self._identifier,
                             Feature(feature.location.start, feature.location.end, feature.strand,
                                     feature.sub_features)))
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

    def _intron_list(self, front_size: int, end_size: int) -> Generator[str, None, None]:
        gff3_ptr = open(self._gff_file, "r")
        record: SeqRecord
        for record in GFF.parse(gff3_ptr,
                                base_dict=Annotation.genome_dict,
                                limit_info=dict(gff_type=self._features)):
            self[record.name] = []
            for feature in record.features:
                sub_features = feature.sub_features
                self[record.name].append(feature)
                if len(sub_features) > 1:
                    intron_indices = []
                    for i in range(len(sub_features) - 1):
                        intron_indices.append((sub_features[i].location.end, sub_features[i + 1].location.start))
                    for index in intron_indices:
                        seq = record.seq[index[0] - 1: index[0] + front_size - 1] + \
                              record.seq[index[1] - end_size + 1: index[1] + 1]
                        if len(seq) < front_size + end_size:
                            continue
                        if feature.strand == 1:
                            yield str(seq).upper()
                        else:
                            yield str(seq.reverse_complement()).upper()
        gff3_ptr.close()
