liftover-utils
==============

Utilities for lifting over genome coordinates in _C. elegans_. Currently supports:

* BCF/VCF
* GFF
* Bed
* Custom (specify chromosome, start position, and optionally end position)

Makes use of `remap_gff_between_releases.pl` by Gary Williams.

	Usage:
	  liftover.py <file> <release1> <release2> (bcf|vcf|gff|bed)
	  liftover.py <file> <release1> <release2> <chrom_col> <start_pos_column> [<end_pos_column>] [options]

	Options:
	  -h --help     Show this screen.
	  --delim=<delim>  File Delimiter; Default is a tab [default: TAB].

## Installation

```
pip install https://github.com/AndersenLab/liftover-utils/archive/v0.1.tar.gz
```
