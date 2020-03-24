from alignment import Fingerprint, FingerprintException
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from collections import deque
import editdistance
import os
import sys
import textdistance

def read_sequence(file):
    seqs = []
    for seq_record in SeqIO.parse(file, "fasta"):
        sys.stdout.write(seq_record.id + "\n")
        sys.stdout.flush()
        sys.stdout.write("Length: {}\n".format(len(seq_record)))
        # print("Sequence: {}".format(seq_record.seq))
        return seq_record.seq


def write_sequence(seq, filename):
    SeqIO.write(seq, filename, "fasta")

'''
Rotates the fingerprint signature by 1 hash each time.
'''
def rotate_fingerprint(fp):
    return fp[1:] + [fp[0]]


'''
Get the rotated sequence corresponding to the given rotated fingerprint.
'''
def get_rotation_by_index(ix, seq):
    print("prefix to append: {}".format(seq[:ix]))
    print("suffix: {}".format(seq[ix:]))
    rot_seq = seq[ix:] + seq[:ix]
    # print("Rotated Seq: {}".format(rot_seq))
    return rot_seq


'''
Returns the Levenshtein edit distance between two sequences of hashes.
'''
def get_fingerprint_dist(lla, llb):
    print(" ============ROTATING FINGERPRINT============")
    a_hashes = [item[0] for item in lla]
    b_hashes = [item[0] for item in llb]
    # print("A: {}\nB:{}".format(a_hashes, b_hashes))
    print("Difference in Length of Fingerprints: {}".format(len(a_hashes) - len(b_hashes)))
    dist = textdistance.levenshtein.distance(a_hashes, b_hashes)
    print("Start Index in Sequence B: {}".format(llb[0]))
    print("Edit Distance of Fingerprints: {}".format(dist))
    return (dist, llb[0][1])

'''
Finds the rotation of b_mins at which the
edit distance between a_mins and b_mins is minimal.
'''
def get_min_dist_rotation(a_mins, b_mins):
    min_dist_rot = 0
    split_point = 0
    min_dist = max(len(a_mins), len(b_mins))
    rot_b_mins = b_mins
    for i in range(len(a_mins)):
        rotation_dist_and_ix = get_fingerprint_dist(a_mins, rot_b_mins)
        rotation_dist = rotation_dist_and_ix[0]
        rotation_ix = rotation_dist_and_ix[1]
        if min_dist >= rotation_dist:
            min_dist_rot =  i
            min_dist = rotation_dist
            split_point = rotation_ix
        # print("Rotation {}; Dist: {} Index: {}".format(i, rotation_dist, rotation_ix))
        rot_b_mins = rotate_fingerprint(rot_b_mins)
    return (min_dist_rot, min_dist, split_point)


if __name__ == '__main__':
    seq_a_file = input("Sequence A File: ")
    seq_a = read_sequence(seq_a_file)

    seq_b_file = input("Sequence B File: ")
    seq_b = read_sequence(seq_b_file)

    baseline_lev_dist = editdistance.eval(str(seq_a), str(seq_b))

    print("Actual Lev Edit Distance Between Sequences: {}".format(baseline_lev_dist))

    afp = Fingerprint(k=5, window_len=5, base=7, modulo=101)
    a_mins = afp.generate(str(seq_a))


    bfp = Fingerprint(k=5, window_len=5, base=7, modulo=101)
    b_mins = bfp.generate(str(seq_b))

    a_fp_len = len(a_mins)
    b_fp_len = len(b_mins)
    len_diff = a_fp_len - b_fp_len

    print("A fingerprint length: {} B fingerprint length: {}".format(a_fp_len, b_fp_len, len_diff))
    # print("A Fingerprint: {}".format(a_mins))
    # print("B Fingerprint: {}".format(b_mins))

    md = get_min_dist_rotation(a_mins, b_mins)
    print("===========Best Rotation of Fingerprint===========")
    print("ROTATION {} with MIN DISTANCE {}; INDEX {}".format(md[0], md[1], md[2]))

    rotation_split_point = md[2]

    rotated_b_seq = get_rotation_by_index(rotation_split_point, seq_b)
    print("Seq A: {}".format(seq_a))
    print("Seq B: {}".format(seq_b))
    # print("Rot B: {}".format(rotated_b_seq))
    # print("Seq A: {}".format(seq_a))
    # for i in range(len(seq_a)):
        # print("{} {} {}".format())

    adj_dist = textdistance.levenshtein.distance(seq_a, rotated_b_seq)
    print("Adjusted Lev Distance: {}".format(adj_dist))
