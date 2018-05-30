# ani-m
Average Nucleotide Identity analysis with MUMmer under the hood

## Usage

```
This script calculates the average nucleotide identity (ANI).
  If multiple queries are given, they will all be compared to
  the reference genome.
  Usage: ani-m.pl reference.fasta query.fasta [query2.fasta...]
  --symmetric        When given, the reference and query will be
                     swapped in a second calculated ANI value.
  --mash-filter 0.9  Use the software mash to skip any ANI
                     calculations when the mash distance
                     is more than this number. If Mash is not
                     found, no filtering will be run.
  --verbose
```

## Requirements

* MUMmer v3.23 or above
* (optional) Mash v2 or above
* Perl
