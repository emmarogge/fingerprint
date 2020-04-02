from alignment import Fingerprint, FingerprintException
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, SeqRecord
from collections import deque
from difflib import SequenceMatcher
import editdistance
import matplotlib.pyplot as plt
import os
import rotation
import sys
import textdistance

save_figure_dir = "/Users/emmarogge/workspace/fingerprint/figures/"

'''
Rotates the fingerprint signature by 1 hash each time.
'''
def rotate_fingerprint(fp):
    return fp[1:] + [fp[0]]

'''
Returns the Levenshtein edit distance between two sequences of hashes.
'''
def get_fingerprint_dist(lla, llb):
    a_hashes = [item[0] for item in lla]
    b_hashes = [item[0] for item in llb]
    dist = textdistance.levenshtein.distance(a_hashes, b_hashes)
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
        # print("Rotation: {} Distance: {}: Index: {}".format(i, rotation_dist_and_ix[0], rotation_dist_and_ix[1]))
        rotation_dist = rotation_dist_and_ix[0]
        rotation_ix = rotation_dist_and_ix[1]
        if min_dist >= rotation_dist:
            min_dist_rot =  i
            min_dist = rotation_dist
            split_point = rotation_ix
        rot_to_plot.append(i)
        indices_to_plot.append(rotation_ix)
        distances_to_plot.append(rotation_dist)
        rot_b_mins = rotate_fingerprint(rot_b_mins)
    plt.plot(rot_to_plot, distances_to_plot)
    plt.xlabel("Rotation")
    plt.ylabel("Lev. Dist. Between Fingerprints")
    plt.title("{} vs {}, k = {}, w = {}".format(a_name, b_name, k, w))
    fig_title = save_figure_dir + 'graph_' + a_name + '_' + b_name + '_k_' + str(k) + '_w_' + str(w) + '.png'
    plt.savefig(fig_title, dpi = 300)
    return (min_dist_rot, min_dist, split_point)

def refine_best_split_point(seq_a, rot_seq, w, b_len):
    print("================================")
    a_len = len(seq_a.seq)
    a_win = seq_a.seq[:(w)]

    b_len = len(rot_seq.seq)
    b_win = rot_seq.seq[b_len-(w):b_len] + rot_seq.seq[:(w)]

    print("Seq A - {}".format(a_win))
    print("Rot B - {}".format(b_win))

    s = SequenceMatcher(None, a_win, b_win)
    # print("len awin {} len bwin {}".format(len(a_win), len(b_win)))
    match = s.find_longest_match(0, len(a_win), 0, len(b_win))
    print("Match: {} {} {}".format(match[0], match[1], match[2]))
    refined_start_index = None
    # If the start index is less than the refining window size
    if match[1] < w:
        # We rotate counter-clockwise
        refined_start_index = b_len - match[1] - 1
    else:
        # We rotate clockwise
        refined_start_index = match[1]
    return refined_start_index


def caLSH(k_val, w_size, seq_a_file, seq_b_file):

    # Uploading biological sequences to be aligned
    seq_a = rotation.read_sequence(seq_a_file)
    a_name = seq_a.id

    seq_b = rotation.read_sequence(seq_b_file)
    b_name = seq_b.id
    b_len = len(seq_b.seq)
    print("b length: {}".format(b_len))
    # Calculate baseline Levenshtein edit distance between sequences.
    baseline_lev_dist = editdistance.eval(str(seq_a.seq), str(seq_b.seq))
    print("Edit Distance w/ No Rotation: {}".format(baseline_lev_dist))

    # Create fingerprints of original sequences.
    afp = Fingerprint(k=int(k_val), window_len=int(w_size), base=7, modulo=101)
    a_mins = afp.generate(str(seq_a.seq))
    bfp = Fingerprint(k=int(k_val), window_len=int(w_size), base=7, modulo=101)
    b_mins = bfp.generate(str(seq_b.seq))

    # Find rotation of Fingerprint B with minimum edit distance from Fingerprint A.
    md = get_min_dist_rotation(a_mins, b_mins, a_name, b_name, k_val, w_size)
    print("ROTATION {} with MIN DISTANCE {}; INDEX {}".format(md[0], md[1], md[2]))

    rotation_split_point = md[2]

    # Get the rotation of the original sequence corresponding to best fingerprint rotation.
    rotated_b_seq = rotation.get_rotation_by_index(rotation_split_point, seq_b)
    print("Rotation: {}".format(str(rotated_b_seq)))

    start_index = refine_best_split_point(seq_a, rotated_b_seq, w_size, b_len)
    print("Refined Start Index: {}".format(start_index))

    refined_rot_b_seq = rotation.get_rotation_by_index(start_index, rotated_b_seq)
    print("Refined: {}".format(refined_rot_b_seq.seq))

    # Compute Levenshtein edit distance of sequence A and rotated sequence B.
    adj_dist = textdistance.levenshtein.distance(seq_a.seq, rotated_b_seq)
    print("Coarse Step - Adjusted Lev Distance: {}".format(adj_dist))

    # Save the rotation of Sequence B in FASTA format.
    rot_b_file = "coarse_step_rot_{}.fasta".format(b_name)
    SeqIO.write(rotated_b_seq, rot_b_file, "fasta")

    # Compute Levenshtein edit distance of sequence A and refined, rotated sequence B.
    adj_dist_refined = textdistance.levenshtein.distance(seq_a.seq, refined_rot_b_seq)
    print("Refined Step - Adjusted Lev Distance: {}".format(adj_dist_refined))

    # Save the rotation of Sequence B in FASTA format.
    refined_rot_b_file = "refined_rot_{}_for_{}.fasta".format(b_name, a_name)
    SeqIO.write(refined_rot_b_seq, refined_rot_b_file, "fasta")

if __name__ == "__main__":
    # Setting parameters for kgram length & MinHash window size
    k = int(input("Input a kgram length: "))
    w = int(input("Input a window size: "))

    # Uploading biological sequences to be aligned
    seq_a_file = input("Sequence A filepath: ")
    seq_b_file = input("Sequence B filepath: ")

    caLSH(k, w, seq_a_file, seq_b_file)

