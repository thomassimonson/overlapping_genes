import argparse
import re

import seq_tools as st
import alignement_weight as aw
from overlapping import arrangement, parse_arguments, read_gencode

# Description #################################################################
# This script perform the exact optimization of more than 2 protein sequences
# that are overlapping.
# Description #################################################################


class Quadon_t:
    """This class represents the quadruplet of nucleotides that codes for more than 3 sequence
    """

    def __init__(self, res_list, position, ext_5, ext_3, ent_list):
        "generalization of the gene overlapping code"
        self.ext_5 = ext_5
        self.ext_3 = ext_3
        self.res_list = res_list
        self.position = position
        self.score, self.quadon = self.optimize_midle(res_list, ent_list)
        self.before = None
        self.after = None
        # TODO: get residues
        self.protein = [
            tuple(translation(self.quadon, frame) for frame in FRAME_LIST)
        ]
        self.sequence = self.quadon

    def optimize_midle(self, res_list, ent_list):
        """Get the optimal quadruplet of nucleotides that optimize blosum score.
        """
        if self.ext_5 is None:
            results = [
                blosum(midle + self.ext_3, res_list, ent_list)
                for midle in arrangement(3)
            ]
        elif self.ext_3 is None:
            results = [
                blosum(self.ext_5 + midle, res_list, ent_list)
                for midle in arrangement(3)
            ]
        else:
            results = [
                blosum(self.ext_5 + midle + self.ext_3, res_list, ent_list)
                for midle in arrangement(2)
            ]
        return max(results, key=lambda (score, quadon): score)

    def set_cur_sequence(self, pre_quadon):
        """append to the beginning of the current sequence.
        """
        self.protein = pre_quadon.protein + self.protein
        self.sequence = pre_quadon.sequence[:-1] + self.sequence


def blosum(quadon, res_list, ent_list):
    """Return the quadon + similarity score.
    """
    residues = [translation(quadon, frame) for frame in FRAME_LIST]
    if None not in residues:
        score = lambda res_p, res_c, res_ent: SCORE_MATRIX[res_p][res_c] * res_ent
        return sum([
            score(res_prime, res, ent)
            for res_prime, res, ent in zip(residues, res_list, ent_list)
        ]), quadon
    else:
        return -INF, None


def cure_cur_position(cur_position):
    """Remove from the current list all the none optimal quadon.
    """
    results = {}
    for nuc in NUCLEOTIDE:
        results[nuc] = max(
            [position for position in cur_position if position.ext_3 == nuc],
            key=lambda pos: pos.score)
    return results


def translation(quadon, frame):
    """Depending on the frame, it will translate the quadon into two
    overlapping sequences.
    """
    if frame == 1:
        x_dna = quadon[:-1]
    elif frame == 2:
        x_dna = quadon[1:]
    elif frame == 3:
        x_dna = quadon[:-1][::-1]
    else:
        x_dna = quadon[1:][::-1]
    try:
        complement = st.complement if frame in [3, 4] else lambda nuc: nuc
        x_protein = "".join(
            [GENCODE[complement(codon)] for codon in st.split(x_dna)])
        return x_protein
    except KeyError:
        return None


def graph_optimization(pair_list, entropy):
    """Depending on the sign of the frame, create an overlapped pair of sequence x
    and y. Thus, it will create a graph of quadruplets of nucleotides
    (quadons). Both sequences are zipped into one list of tuples.
    """
    # Graph initialization for both side (first and last position)
    prev_position = ext_graph(pair_list.pop(0), entropy.pop(0))
    last_position = ext_graph(pair_list.pop(-1), entropy.pop(-1), False)

    # Bind corresponding 5' to 3'
    for position, res in enumerate(pair_list):
        cur_position = [
            connexion(
                Quadon_t(res, position, ext[0], ext[1], entropy[position]),
                prev_position) for ext in arrangement(2)
        ]
        prev_position = cure_cur_position(cur_position)

    # binding the last position
    for nuc in NUCLEOTIDE:
        last_position[nuc].score += prev_position[nuc].score
        last_position[nuc].set_cur_sequence(prev_position[nuc])

    best_quadon = max([last_position[nuc] for nuc in NUCLEOTIDE],
                      key=lambda quadon: quadon.score)
    return best_quadon


def connexion(cur_quadon, pre_quadons):
    """Connect all the previous quadons to the current quadon.
    """
    score = lambda c_quadon, p_quadon: c_quadon.score + p_quadon.score
    cur_quadon.score = score(cur_quadon, pre_quadons[cur_quadon.ext_5])
    cur_quadon.set_cur_sequence(pre_quadons[cur_quadon.ext_5])
    return cur_quadon


def ext_graph(res, entropy, init=True):
    """Graph initialization if init = True. We will use a list to represent this kind of graph.
    Otherwise it is the last position of pairs.
    """
    if init:
        graph = {
            nuc: Quadon_t(res, 0, None, nuc, entropy)
            for nuc in NUCLEOTIDE
        }
    else:
        graph = {
            nuc: Quadon_t(res, -1, nuc, None, entropy)
            for nuc in NUCLEOTIDE
        }
    return graph


def overlapping_g(sequence_list, ent_list=None):
    """This function perform the optimization of 2 sequences.
    """
    if ent_list is None:
        ent_list = [[1 for _ in range(len(seq))] for seq in sequence_list]

    for index, frame in enumerate(FRAME_LIST):
        if frame in [3, 4]:
            sequence_list[index] = sequence_list[index][::-1]
            ent_list[index] = ent_list[index][::-1]
    best = graph_optimization(zip(*sequence_list), zip(*ent_list))
    res = []
    for index, frame in enumerate(FRAME_LIST):
        if frame in [3, 4]:
            res.append(
                "".join("".join([el[index] for el in best.protein])[::-1]))
        else:
            res.append("".join([el[index] for el in best.protein]))
    return best


def main():
    """This module perform the triple gene overlapping.
    """
    # Initialize variables ####################################################
    global GENCODE, FRAME_LIST, SCORE_MATRIX, NUCLEOTIDE, STOPS, INF, ENTROPY, PFAM_DIR, GROUPS
    args = parse_arguments()
    PFAM_DIR = args.pfam_dir
    GROUPS = aw.get_group(args.group)
    GENCODE = read_gencode(args.gencode)
    SCORE_MATRIX = st.read_score(args.score_matrix)
    NUCLEOTIDE = ['A', 'C', 'T', 'G']
    STOPS = ["TAA", "TAG", "TGA"]
    INF = float("inf")
    FRAME_LIST = args.frame_list
    ENTROPY = True if args.entropy is not None else False
    fasta_seq = st.read_fasta_file(args.infile)
    test = lambda x_s, y_s: sum(SCORE_MATRIX[x][y] for x, y in zip(x_s, y_s))
    out = open(args.out_file, "w") if args.out_file is not None else False
    # Initialize variables ####################################################

    # Getting sequences #######################################################
    for name_x, seq_x_init in st.read_fasta_file(args.infile):
        for name_y, seq_y_init in st.read_fasta_file(args.infile):
            for name_z, seq_z_init in st.read_fasta_file(args.infile):
                # Getting weights from alignments + optimization ##########################
                if ENTROPY:
                    x_entropy = aw.read_alignement(name_x, PFAM_DIR, GROUPS)
                    y_entropy = aw.read_alignement(name_y, PFAM_DIR, GROUPS)
                    z_entropy = aw.read_alignement(name_z, PFAM_DIR, GROUPS)
                    best = overlapping_g([seq_x_init, seq_y_init, seq_z_init], [x_entropy, y_entropy, z_entropy])
                else:
                    best = overlapping_g([seq_x_init, seq_y_init, seq_z_init])

                # Sorting sequence (5' -> 3') #############################################
                res = []
                for index, frame in enumerate(FRAME_LIST):
                    if frame in [3, 4]:
                        res.append("".join("".join(
                            [el[index] for el in best.protein])[::-1]))
                    else:
                        res.append("".join([el[index] for el in best.protein]))

                # Pretty print ############################################################
                if args.verbose is not None or not out:
                    print "=============================================================="
                    print "SCORE :", best.score, "/ FRAME :", ",".join(map(str, FRAME_LIST))
                    print "{} => similarity X, X' = {}, frame position = {}".format(name_x, test(seq_x_init, res[0]), FRAME_LIST[0])
                    print "X : {}\nX': {}".format(seq_x_init, res[0])
                    print "{} => similarity Y, Y' = {}, frame position = {}".format(name_y, test(seq_y_init, res[1]), FRAME_LIST[1])
                    print "Y : {}\nY': {}".format(seq_y_init, res[1])
                    print "{} => similarity Z, Z' = {}, frame position = {}".format(name_z, test(seq_z_init, res[2]), FRAME_LIST[2])
                    print "Z : {}\nZ': {}\n".format(seq_z_init, res[2])
                    st.print_sequence_3(seq_x_init, seq_y_init, seq_z_init, res[0], res[1], res[2], best.sequence, FRAME_LIST)
                    print "\n=> Press enter to pursue"
                    raw_input()
                if out:
                    seq_name = ">{}_Vs_{}_Vs_{}".format(name_x, name_y, name_z)
                    out.write(seq_name+"_X\n{}\n".format(st.write_fasta_format(res[0])))
                    out.write(seq_name+"_Y\n{}\n".format(st.write_fasta_format(res[1])))
                    out.write(seq_name+"_Z\n{}\n".format(st.write_fasta_format(res[2])))
    if out:
        out.close()
    return 0


if __name__ == '__main__':
    main()
