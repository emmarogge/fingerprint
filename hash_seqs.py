from alignment import Fingerprint, FingerprintException
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from collections import deque
import editdistance
import matplotlib.pyplot as plt
import os
import sys
import textdistance

save_figure_dir = "/Users/emmarogge/workspace/fingerprint/figures/"

def read_sequence(file):
    seqs = []
    for seq_record in SeqIO.parse(file, "fasta"):
        sys.stdout.write(seq_record.id + "\n")
        sys.stdout.flush()
        sys.stdout.write("Length: {}\n".format(len(seq_record)))
        # print("Sequence: {}".format(seq_record.seq))
        return seq_record


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
def get_min_dist_rotation(a_mins, b_mins, a_name, b_name, k, w):
    min_dist_rot = 0
    split_point = 0
    min_dist = max(len(a_mins), len(b_mins))
    rot_b_mins = b_mins
    rot_to_plot = []
    indices_to_plot = []
    distances_to_plot = []
    for i in range(len(a_mins)):
        rotation_dist_and_ix = get_fingerprint_dist(a_mins, rot_b_mins)
        rotation_dist = rotation_dist_and_ix[0]
        rotation_ix = rotation_dist_and_ix[1]
        if min_dist >= rotation_dist:
            min_dist_rot =  i
            min_dist = rotation_dist
            split_point = rotation_ix
        # print("Rotation {}; Dist: {} Index: {}".format(i, rotation_dist, rotation_ix))
        rot_to_plot.append(i)
        indices_to_plot.append(rotation_ix)
        distances_to_plot.append(rotation_dist)
        rot_b_mins = rotate_fingerprint(rot_b_mins)
    plt.plot(rot_to_plot, distances_to_plot, 'r')
    plt.xlabel("Rotation")
    plt.ylabel("Levenshtein Distance")
    plt.title("k = {}, w = {}".format(k, w))
    fig_title = save_figure_dir + 'graph_' + a_name + '_' + b_name + '_k_' + str(k) + '_w_' + str(w) + '.png'
    plt.savefig(fig_title, dpi = 300)
    return (min_dist_rot, min_dist, split_point)

if __name__ == '__main__':

    # Setting parameters for kgram length & MinHash window size
    k_val = input("Input a kgram length: ")
    w_size = input("Input a window size: ")

    # Uploading biological sequences to be aligned
    seq_a_file = input("Sequence A filepath: ")
    seq_a = read_sequence(seq_a_file)
    a_name = seq_a.id

    seq_b_file = input("Sequence B filepath: ")
    seq_b = read_sequence(seq_b_file)
    b_name = seq_b.id

    baseline_lev_dist = editdistance.eval(str(seq_a.seq), str(seq_b.seq))

    print("Lev Edit Distance Between Sequences for Original Split Point: {}".format(baseline_lev_dist))

    _k = int(k_val)
    _w = int(w_size)
    afp = Fingerprint(k=_k, window_len=_w, base=7, modulo=101)
    a_mins = afp.generate(str(seq_a.seq))


    bfp = Fingerprint(k=_k, window_len=_w, base=7, modulo=101)
    b_mins = bfp.generate(str(seq_b.seq))

    a_fp_len = len(a_mins)
    b_fp_len = len(b_mins)
    len_diff = a_fp_len - b_fp_len

    print("A fingerprint length: {} B fingerprint length: {}".format(a_fp_len, b_fp_len, len_diff))
    # print("A Fingerprint: {}".format(a_mins))
    # print("B Fingerprint: {}".format(b_mins))

    md = get_min_dist_rotation(a_mins, b_mins, a_name, b_name, _k, _w)
    # print("===========Best Rotation of Fingerprint===========")
    # print("ROTATION {} with MIN DISTANCE {}; INDEX {}".format(md[0], md[1], md[2]))

    rotation_split_point = md[2]

    rotated_b_seq = get_rotation_by_index(rotation_split_point, seq_b.seq)
    print("Seq A: {}".format(seq_a))
    print("Seq B: {}".format(seq_b))
    # print("Rot B: {}".format(rotated_b_seq))
    # print("Seq A: {}".format(seq_a))
    # for i in range(len(seq_a)):
        # print("{} {} {}".format())

    adj_dist = textdistance.levenshtein.distance(seq_a.seq, rotated_b_seq)
    print("Adjusted Lev Distance: {}".format(adj_dist))
