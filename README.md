# LocusSolver

Accessory script for merging features from GFF3 files. This script is included as part of the YAPIM implementation of EukMetaSanity.

## Installation

`LocusSolver` expects Python3.8.

```shell
pip install git+https://github.com/cjneely10/LocusSolver.git
```

## Usage

```shell
Merge GFF3 files into single annotation set

Usage:
    locus_solver [SWITCHES] 

Meta-switches:
    -h, --help                            Prints this help message and quits
    --help-all                            Print help messages of all subcommands and quit
    -v, --version                         Prints the program's version and quits

Switches:
    -i, --input-gff3 INPUT_GFF3S:str      Path to annotation GFF3 file. Provide in expected priority sort order; may be given
                                          multiple times; required
    -o, --output OUTPUT_PATH:str          Path to output merged results, default stdout
    -t, --tier TIER:int                   Set number of programs that must corroborate an input, default 1

```
Example:

```shell
locus_solver -i data/sample.gmes.gff3 -i data/sample.aug.gff3 -o data/sample.merged.gff3
locus_solver -t 2 -i data/sample.gmes.gff3 -i data/sample.aug.gff3 -o data/sample.merged.gff3
```