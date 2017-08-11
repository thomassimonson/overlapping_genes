import argparse
import re
import numpy as np

import seq_tools as st
import alignement_weight as aw
import overlapping_g as og

# Description #################################################################
# This script perform the exact optimization of 2 protein sequences that are
# overlapping.
# Description #################################################################

class Quadon:
    """This class represents the quadruplet of nucleotides that are needed for a
    pair of overlapping amino-acids.
    """

    def __init__(self, x_res, y_res, position, ext_5, ext_3, x_ent=1, y_ent=1, frame=None):
        """
        Keyword Arguments:
        x_res    -- current residue from X
        y_res    -- current residue from Y
        position -- current position
        ext_5    -- A, T, C or G => 5' start
        ext_3    -- A, T, C or G => 3' end
        x_ent    -- entropy weight for x
        y_ent    -- entropy weight for y
        """
        self.ext_5 = ext_5
        self.ext_3 = ext_3
        self.x_res = x_res
        self.y_res = y_res
        self.position = position
        self.score, self.quadon = self.optimize_midle(x_res, y_res, x_ent, y_ent)
        self.before = None
        self.after = None
        self.protein = [translation(self.quadon)]
        self.sequence = self.quadon

    def optimize_midle(self, x_res, y_res, x_ent, y_ent):
        """Get the optimal quadruplet of nucleotides that optimize blosum score.
        """
        if self.ext_5 is None and self.ext_3 is None:
            results = [blosum(codon, x_res, y_res, x_ent, y_ent) for codon in arrangement(3)]
        elif self.ext_5 is None:
            results = [blosum(midle+self.ext_3, x_res, y_res, x_ent, y_ent) for midle in arrangement(3)]
        elif self.ext_3 is None:
            results = [blosum(self.ext_5+midle, x_res, y_res, x_ent, y_ent) for midle in arrangement(3)]
        else:
            results = [blosum(self.ext_5+midle+self.ext_3, x_res, y_res, x_ent, y_ent) for midle in arrangement(2)]
        return max(results, key=lambda (score, quadon): score)

    def set_cur_sequence(self, pre_quadon):
        """append to the beginning of the current sequence.
        """
        self.protein = pre_quadon.protein + self.protein
        self.sequence = pre_quadon.sequence[:-1] + self.sequence

    def __repr__(self):
        return """
         quadon : {}\nscore : {}\nposition : {}\nX:{}\nY:{}
               """.format(self.quadon, self.score, self.position,
                          self.protein[-1][0], self.protein[-1][1])


def arrangement(n, results=[], object_list=['A', 'C', 'T', 'G']):
    """Function for arrangement of a list with replacement.
    """
    if n > 0:
        if len(results) == 0:
            results = [elem for elem in object_list]
        else:
            results = ["".join(res+elem) for res in results for elem in object_list]
        return arrangement(n-1, results)
    else:
        return results

def blosum(quadon, x_res, y_res, x_ent, y_ent):
    """Return the quadon + similarity score.
    """
    residues = translation(quadon)
    if residues is not None:
        x_res_prime, y_res_prime = residues
        score = lambda X_p, X_c, Y_p, Y_c: SCORE_MATRIX[X_p][X_c] * x_ent + SCORE_MATRIX[Y_p][Y_c] * y_ent
        return score(x_res_prime, x_res, y_res_prime, y_res), quadon
    else:
        return -INF, None

def read_gencode(infile, excludes=['_']):
    """Read the genetic code from a file.
    """
    regexp = re.compile(r"\S+")
    gencode = {}
    with open(infile) as gencode_file:
        for line in gencode_file:
            val = regexp.findall(line)
            if val[1] not in excludes:
                dna, prot = val[0], val[1]
                gencode[dna] = prot
    return gencode

def graph_optimization(pair_list, x_entropy=None, y_entropy=None):
    """Depending on the sign of the frame, create an overlapped pair of sequence x
    and y. Thus, it will create a graph of quadruplets of nucleotides
    (quadons). Both sequences are zipped into one list of tuples.
    """
    # Graph initialization for both side (first and last position)
    prev_position = ext_graph(pair_list.pop(0), x_entropy.pop(0), y_entropy.pop(0))
    last_position = ext_graph(pair_list.pop(-1), False, x_entropy.pop(-1), y_entropy.pop(-1))

    # Bind corresponding 5' to 3'
    for position, (x_res, y_res) in enumerate(pair_list):
        cur_position = [connexion(Quadon(x_res, y_res, position, ext[0],
                                         ext[1]), prev_position, x_entropy[position],
                                  y_entropy[position]) for ext in arrangement(2)]
        prev_position = cure_cur_position(cur_position)

    # binding the last position
    for nuc in NUCLEOTIDE:
        last_position[nuc].score += prev_position[nuc].score
        last_position[nuc].set_cur_sequence(prev_position[nuc])

    best_quadon = max([last_position[nuc] for nuc in NUCLEOTIDE], key=lambda quadon: quadon.score)
    return best_quadon

def cure_cur_position(cur_position):
    """Remove from the current list all the none optimal quadon.
    """
    results = {}
    for nuc in NUCLEOTIDE:
        results[nuc] = max([position for position in cur_position if
                            position.ext_3 == nuc], key=lambda pos: pos.score)
    return results

def ext_graph((x_init_res, y_init_res), init=True, x_ent=1, y_ent=1):
    """Graph initialization if init = True. We will use a list to represent this kind of graph.
    Otherwise it is the last position of pairs.
    """
    if init:
        graph = {nuc: Quadon(x_init_res, y_init_res, 0, None, nuc, x_ent, x_ent) for nuc in NUCLEOTIDE}
    else:
        graph = {nuc: Quadon(x_init_res, y_init_res, -1, nuc, None, x_ent, y_ent) for nuc in NUCLEOTIDE}
    return graph

def connexion(cur_quadon, pre_quadons, x_ent, y_ent):
    """Connect all the previous quadons to the current quadon.
    """
    score = lambda c_quadon, p_quadon: c_quadon.score + p_quadon.score
    cur_quadon.score = score(cur_quadon, pre_quadons[cur_quadon.ext_5])
    cur_quadon.set_cur_sequence(pre_quadons[cur_quadon.ext_5])
    return cur_quadon

def translation(dna_sequence):
    """Depending on the frame, it will translate the dna_sequence into two
    overlapping sequences.
    """
    if FRAME == -2:
        x_dna = dna_sequence[:-1]
        y_dna = dna_sequence[1:][::-1]
    elif FRAME == -1:
        x_dna = dna_sequence[1:]
        y_dna = dna_sequence[:-1][::-1]
    elif FRAME == 1:
        x_dna = dna_sequence[:-1]
        y_dna = dna_sequence[1:]
    elif FRAME == 2:
        x_dna = dna_sequence[1:]
        y_dna = dna_sequence[:-1]
    else:
        x_dna = dna_sequence
        y_dna = dna_sequence[::-1]
    try:
        y_complement = st.complement if FRAME <= 0 else lambda nuc: nuc
        x_protein = "".join([GENCODE[codon] for codon in st.split(x_dna)])
        y_protein = "".join([GENCODE[y_complement(codon)] for codon in st.split(y_dna)])
        return x_protein, y_protein
    except KeyError:
        # Means, it is a stop codon
        return None

def optimize_zero(pair_list_x_y, x_entropy, y_entropy):
    """This function perform the optimization for the special case, frame = 0.
    """
    results = []
    for position, (x_res, y_res) in enumerate(pair_list_x_y):
        results.append(Quadon(x_res, y_res, position, None, None,
                              x_entropy[position], y_entropy[position]))
    x_protein = "".join([quadon.protein[0][0] for quadon in results])
    y_protein = "".join([quadon.protein[0][1] for quadon in results])[::-1]
    dna = "".join([quadon.quadon for quadon in results])
    return x_protein, y_protein, dna

def overlapping(x_sequence, y_sequence, frame, x_entropy=None, y_entropy=None):
    """This function perform the optimization of 2 sequences.
    """
    if x_entropy is None or y_entropy is None:
        x_entropy = [1 for _ in x_sequence]
        y_entropy = [1 for _ in y_sequence]

    if FRAME <= 0:
        # y is the anti-sens strand
        pair_list = zip(x_sequence, y_sequence[::-1])
        y_entropy = y_entropy[::-1]
    else:
        # sens => same direction for each strand
        pair_list = zip(x_sequence, y_sequence)

    if frame == 0:
        return optimize_zero(pair_list, x_entropy, y_entropy)
    else:
        best_path = graph_optimization(pair_list, x_entropy, y_entropy)
        x_protein = "".join([res[0] for res in best_path.protein])
        if FRAME < 0:
            y_protein = "".join([res[1] for res in best_path.protein[::-1]])
        else:
            y_protein = "".join([res[1] for res in best_path.protein])
        return x_protein, y_protein, best_path.sequence

def parse_arguments():
    """Parsing command line
    """
    parser = argparse.ArgumentParser(description="""This program is an implementation of an exact optimization for overlapping
    genes. From a fasta file, each sequence will overlap the others.""")
    parser.add_argument('-g', '--gencode', dest='gencode',
                        default='./data/gencode.dat')
    parser.add_argument('-i', '--infile',
                        dest='infile', help='The fasta file')
    parser.add_argument('-m', '--matrix', dest='score_matrix',
                        default='./data/BLOSUM62', help='the scoring matrix')
    parser.add_argument('-f', '--frame',
                        dest='frame', default=-2, type=int)
    parser.add_argument('-d', '--pfam_dir',
                        dest='pfam_dir', default='./data/pfam_family', help='directory that contains all the alignments files')
    parser.add_argument('-gr', '--group',
                        dest='group', default='./data/dict_aa', help='file that make correspond amino acids to groups')
    parser.add_argument('-e', '--entropy', dest='entropy', help='any value will make the program considering the entropy of the alignment')
    parser.add_argument('-v', '--verbose', dest='verbose', help='pretty print')
    parser.add_argument('-o', '--out_file', dest='out_file')
    parser.add_argument('-fl', '--frame_list', dest='frame_list', nargs="*", type=int)
    return parser.parse_args()

def main():
    """Perform the optimal algorithm for a fasta file
    """
    # Initialize Variables ####################################################
    global GENCODE, FRAME, SCORE_MATRIX, NUCLEOTIDE, STOPS, INF, PFAM_DIR, GROUPS
    global ENTROPY
    args = parse_arguments()
    FRAME = args.frame
    GENCODE = read_gencode(args.gencode)
    SCORE_MATRIX = st.read_score(args.score_matrix)
    NUCLEOTIDE = ['A', 'C', 'T', 'G']
    STOPS = ["TAA", "TAG", "TGA"]
    INF = float("inf")
    PFAM_DIR = args.pfam_dir
    GROUPS = aw.get_group(args.group)
    ENTROPY = True if args.entropy is not None else False
    out = open(args.out_file, "w") if args.out_file is not None else False
    test = lambda x_s, y_s: np.sum(SCORE_MATRIX[x][y] for x, y in zip(x_s, y_s))
    # Initialize Variables ####################################################

    # Perform triple genes ####################################################
    if args.frame_list is not None:
        return og.main()

    # Perform pair genes ######################################################
    for x_fasta in st.read_fasta_file(args.infile):
        for y_fasta in st.read_fasta_file(args.infile):
            name_x, seq_init_x = x_fasta
            name_y, seq_init_y = y_fasta
            # Getting weights from alignments + optimization ##################
            if ENTROPY:
                x_entropy = aw.read_alignement(name_x, PFAM_DIR, GROUPS)
                y_entropy = aw.read_alignement(name_y, PFAM_DIR, GROUPS)
                x_sequence, y_sequence, dna = overlapping(seq_init_x,
                                                          seq_init_y, FRAME, x_entropy, y_entropy)
            else:
                x_sequence, y_sequence, dna = overlapping(seq_init_x, seq_init_y, FRAME)

            # Compute similarity ##############################################
            similarity_x = test(x_sequence, seq_init_x)
            similarity_y = test(y_sequence[::-1], seq_init_y[::-1]) if FRAME <= 0 else test(y_sequence, seq_init_y)

            # Pretty print or write out put file ##############################
            if args.verbose is not None or not out:
                print "=============================================================="
                print "SCORE : {}".format(similarity_x + similarity_y), "/ FRAME : {}".format(FRAME)
                print "{} => similarity X, X' = {}".format(name_x, similarity_x)
                print "X : {}\nX': {}".format(seq_init_x, x_sequence)
                print "{} => similarity Y, Y' = {}".format(name_y, similarity_y)
                print "Y : {}\nY': {}\n".format(seq_init_y, y_sequence)
                st.print_sequence(seq_init_x, seq_init_y, x_sequence, y_sequence, dna, FRAME)
                print "\n=> Press enter to pursue"
                raw_input()
            if out:
                seq_name = ">{}_Vs_{}".format(name_x, name_y)
                out.write(seq_name+"_X\n{}\n".format(st.write_fasta_format(x_sequence)))
                out.write(seq_name+"_Y\n{}\n".format(st.write_fasta_format(y_sequence)))
    if out:
        out.close()
    return 0

if __name__ == '__main__':
    main()
