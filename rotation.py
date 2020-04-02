from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, SeqRecord
import os
import sys
import textdistance

''' Read biological sequence from FASTA file.'''
def read_sequence(file):
    seqs = []
    for seq_record in SeqIO.parse(file, "fasta"):
        print("================================")
        sys.stdout.write(seq_record.id + "\n")
        sys.stdout.flush()
        sys.stdout.write("Length: {}\n".format(len(seq_record)))
        print("Sequence: {}".format(seq_record.seq))
        print("================================")
        return seq_record

'''Get the rotated sequence corresponding to the given rotated fingerprint.'''
def get_rotation_by_index(ix, seq):
    rot_seq = seq[ix:] + seq[:ix]
    return rot_seq

'''Returns rotated sequence for given file and index'''
def get_rotated_sequence(file,index):
    seq = read_sequence(file)
    output_file = "{}_rot_{}.fasta".format(str(seq.id), str(index))
    rot_seq = get_rotation_by_index(index, seq)
    print("Rotated Sequence:\n {}".format(rot_seq.seq))
    SeqIO.write(rot_seq, output_file, "fasta")

if __name__ == "__main__":
    file = input("File: ")
    index = int(input("Index: "))
    get_rotated_sequence(file,index)
