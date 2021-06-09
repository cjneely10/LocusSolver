from pathlib import Path
from typing import Sequence, Optional, List, Generator

from BCBio import GFF
from Bio import SeqIO

from src.position_prob_matrix import PositionProbabilityMatrix


class Annotation(dict):
    genome_dict: dict = None

    def __init__(self,
                 genome_file: Path,
                 gff3_file: Path,
                 features: Optional[Sequence[str]] = None,
                 front_size: int = 10,
                 end_size: int = 20,
                 *args, **kwargs):
        if features is None:
            features = ["transcript", "exon"]
        self._gff_file = gff3_file
        self._features = features
        if Annotation.genome_dict is None:
            Annotation.load_genome(genome_file)
        super().__init__(*args, **kwargs)
        self._intron_model = PositionProbabilityMatrix(front_size + end_size, self._intron_list(front_size, end_size))

    @property
    def intron_model(self):
        return self._intron_model

    @staticmethod
    def load_genome(genome_file: Path):
        genome_ptr = open(genome_file, "r")
        Annotation.genome_dict = SeqIO.to_dict(SeqIO.parse(genome_ptr, "fasta"))
        genome_ptr.close()

    def _intron_list(self, front_size: int, end_size: int) -> Generator[str, None, None]:
        gff3_ptr = open(self._gff_file, "r")
        for record in GFF.parse(gff3_ptr, base_dict=Annotation.genome_dict, limit_info=dict(gff_type=self._features)):
            for feature in record.features:
                sub_features = feature.sub_features
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
