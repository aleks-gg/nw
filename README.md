## Prepare the environment
Install pixi: ```curl -fsSL https://pixi.sh/install.sh | bash ```  
Install the environment: ```pixi install```

## Run the script
```
usage: nw.py [-h] [--output OUTPUT] [--match-score MATCH_SCORE] [--mismatch-score MISMATCH_SCORE] [--indel-score INDEL_SCORE] fasta_path

Perform global alignment of first two sequences in a fasta file using the Needleman-Wunsch algorithm.

positional arguments:
  fasta_path            Path to the FASTA file containing the sequences

options:
  -h, --help            show this help message and exit
  --output OUTPUT       Path to use for the output file, default = nw_out.fasta
  --match-score MATCH_SCORE
                        Score for a match, default = 1
  --mismatch-score MISMATCH_SCORE
                        Score for a mismatch, default = -1
  --indel-score INDEL_SCORE
                        Score for an indel, default = -1
```
Example:  
```
pixi run python3 nw.py seqs.fasta --output "test.fasta" --match 1 --mismatch -1 --indel -2
```

### Input and output from above example
Script output:
```
Alignment score: -5
```
seqs.fasta:
```
>a, test_a
GTCGACGCA
>b, test_b
GATTACAATTCGA
```
test.fasta:
```
>a, test_a
G--T-CGACGC-A
>b, test_b
GATTACAATTCGA
```
