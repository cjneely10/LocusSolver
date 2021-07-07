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
    locus_solver [SWITCHES] gff3_files...

Meta-switches:
    -h, --help                        Prints this help message and quits
    --help-all                        Print help messages of all subcommands and quit
    -v, --version                     Prints the program's version and quits

Switches:
    -o, --output OUTPUT_PATH:str      Path to output merged results, default stdout
    -t, --tier TIER:int               Set number of programs that must corroborate an input, default 1


```

Provide GFF3 files in order of precedence (e.g. first GFF3 features have higher priority than those in second file, etc.)

Example:

Ensure that your files are in GFF3 format using GFFread:

```shell
gffread "<file>" -G -o "<file>.gff3"
```

Run `locus_solver` on GFF3 files:

```shell
locus_solver data/sample.gmes.gff3 data/sample.aug.gff3 -o data/sample.merged.test.gff3
locus_solver -t 2 data/sample.gmes.gff3 data/sample.aug.gff3 -o data/sample.merged.test.gff3
```