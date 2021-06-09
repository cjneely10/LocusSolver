from pathlib import Path
from typing import Sequence, Optional, List, Tuple

from BCBio import GFF
from Bio import SeqIO


class Annotation(dict):
    genome_dict: dict = None

    def __init__(self, genome_file: Path, gff3_file: Path, features: Optional[Sequence[str]] = None, *args, **kwargs):
        if features is None:
            features = ["transcript", "exon"]
        if Annotation.genome_dict is None:
            Annotation.load_genome(genome_file)
        super().__init__(*args, **kwargs)
        self._load(gff3_file, features)

    @staticmethod
    def load_genome(genome_file: Path):
        genome_ptr = open(genome_file, "r")
        Annotation.genome_dict = SeqIO.to_dict(SeqIO.parse(genome_ptr, "fasta"))
        genome_ptr.close()

    def _load(self, gff3_file: Path, features: Optional[Sequence[str]]):
        gff3_ptr = open(gff3_file, "r")
        for record in GFF.parse(gff3_ptr, base_dict=Annotation.genome_dict, limit_info=dict(gff_type=features)):
            self[record.name] = record
        gff3_ptr.close()

    def intron_list(self, front_size: int = 10, end_size: int = 20) -> List[str]:
        out = []
        for key, val in self.items():
            for feature in val.features:
                sub_features = feature.sub_features
                if len(sub_features) > 1:
                    intron_indices = []
                    for i in range(len(sub_features) - 1):
                        intron_indices.append((sub_features[i].location.end, sub_features[i + 1].location.start))
                    for index in intron_indices:
                        seq = val.seq[index[0] - 1: index[0] + front_size - 1] + \
                              val.seq[index[1] - end_size + 1: index[1] + 1]
                        if len(seq) < 30:
                            continue
                        if feature.strand == 1:
                            out.append(str(seq).upper())
                        else:
                            out.append(str(seq.reverse_complement()).upper())
        return out
