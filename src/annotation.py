from pathlib import Path
from typing import List, Sequence, Optional

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Annotation(dict):
    def __init__(self, genome_file: Path, gff3_file: Path, features: Optional[Sequence[str]] = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._load(genome_file, gff3_file, features)

    def _load(self, genome_file: Path, gff3_file: Path, features: Optional[Sequence[str]]):
        limit_info = {}
        if features is not None:
            limit_info = {"limit_info": dict(gff_type=features)}
        genome_ptr = open(genome_file, "r")
        genome_dict = SeqIO.to_dict(SeqIO.parse(genome_ptr, "fasta"))
        genome_ptr.close()
        gff3_ptr = open(gff3_file, "r")
        record: SeqRecord
        for record in GFF.parse(gff3_ptr, base_dict=genome_dict, **limit_info):
            self[record.name] = record
        gff3_ptr.close()

    def intron_list(self):
        for key, val in self.items():
            for feature in val.features:
                print(feature.sub_features)
