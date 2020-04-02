from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import os
import sys
import textdistance

def read_sequence(file):
    seqs = []
    for seq_record in SeqIO.parse(file, "fasta"):
        sys.stdout.write(seq_record.id + "\n")
        sys.stdout.flush()
        sys.stdout.write("Length: {}\n".format(len(seq_record)))
        print("Sequence: {}".format(seq_record.seq))
        return seq_record


def write_sequence(seq, filename):
    SeqIO.write(seq, filename, "fasta")
    print("written")

'''
Get the rotated sequence corresponding to the given rotated fingerprint.
'''
def get_rotation_by_index(ix, seq):
    # print("prefix to append: {}".format(seq[:ix]))
    # print("suffix: {}".format(seq[ix:]))
    rot_seq = seq[ix:] + seq[:ix]
    # print("Rotated Seq: {}".format(rot_seq))
    return rot_seq

'''Returns rotated sequence for given file and index'''
def get_rotated_sequence():
    file = input("File: ")
    index = int(input("Index: "))

    seq = read_sequence(file)
    output_file = "{}_rot_{}".format(str(seq.id), str(index))
    rot_seq = get_rotation_by_index(index, seq)
    print("Rotated Sequence:\n {}".format(rot_seq.seq))
    write_sequence(rot_seq, output_file)

get_rotated_sequence()
