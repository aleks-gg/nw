import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

class NWAligner:
    def __init__(self, match_score: int = 1, mismatch_score: int = -1, indel_score: int = -1):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.indel_score = indel_score
        
    def align(self, seq_a: str, seq_b: str) -> tuple[str, str, int]:
        cols_num = len(seq_a) + 1
        rows_num = len(seq_b) + 1

        matrix = np.zeros((rows_num, cols_num), dtype=int)
        
        for i in range(1, cols_num):
            matrix[0][i] = i * self.indel_score
        for i in range(1, rows_num):
            matrix[i][0] = i * self.indel_score
        
        for i in range(1, rows_num):
            for j in range(1, cols_num):
                score_top_left = matrix[i-1][j-1] + (self.match_score if seq_a[j-1] == seq_b[i-1] else self.mismatch_score)
                score_top = matrix[i-1][j] + self.indel_score
                score_left = matrix[i][j-1] + self.indel_score
                
                matrix[i][j] = max(score_top_left, score_top, score_left)
        
        seq_a_aligned = ""
        seq_b_aligned = ""  
        col = cols_num - 1
        row = rows_num - 1
        alignment_score = 0
        
        while col > 0 or row > 0:
            current_score = matrix[row][col]
            
            if row > 0 and col > 0 and current_score == matrix[row-1][col-1] + (self.match_score if seq_a[col-1] == seq_b[row-1] else self.mismatch_score): # top-left
                seq_a_aligned += seq_a[col-1]
                seq_b_aligned += seq_b[row-1]
                alignment_score += (self.match_score if seq_a[col-1] == seq_b[row-1] else self.mismatch_score)
                row -= 1
                col -= 1
            elif col > 0 and current_score == matrix[row][col-1] + self.indel_score: # left
                seq_a_aligned += seq_a[col-1]
                seq_b_aligned += '-'
                alignment_score += self.indel_score
                col -= 1
            elif row > 0 and current_score == matrix[row-1][col] + self.indel_score: # top
                seq_a_aligned += '-'
                seq_b_aligned += seq_b[row-1]
                alignment_score += self.indel_score
                row -= 1
        
        seq_a_aligned = seq_a_aligned[::-1]
        seq_b_aligned = seq_b_aligned[::-1]
        
        return seq_a_aligned, seq_b_aligned, alignment_score
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Perform global alignment of first two sequences in a fasta file using the Needleman-Wunsch algorithm.")
    parser.add_argument("fasta_path", type=str, help="Path to the FASTA file containing the sequences")
    parser.add_argument("--output", type=str, default="nw_out.fasta", help="Path to use for the output file")
    parser.add_argument("--match-score", type=int, default=1, help="Score for a match")
    parser.add_argument("--mismatch-score", type=int, default=-1, help="Score for a mismatch")
    parser.add_argument("--indel-score", type=int, default=-1, help="Score for an indel")
    
    args = parser.parse_args()
    
    records = []
    for rec in SeqIO.parse(args.fasta_path, "fasta"):
        records.append(rec)
    assert len(records) > 1
    seq_a = records[0]
    seq_b = records[1]
    
    nw = NWAligner(args.match_score, args.mismatch_score, args.indel_score)
    seq_a_aligned, seq_b_aligned, score = nw.align(seq_a.seq, seq_b.seq)
    
    print("Alignment score:", score)
    
    record_a_out = SeqRecord(Seq(seq_a_aligned), id=seq_a.id, description=seq_a.description)
    record_b_out = SeqRecord(Seq(seq_b_aligned), id=seq_b.id, description=seq_b.description)
    SeqIO.write([record_a_out, record_b_out], args.output, "fasta")
    