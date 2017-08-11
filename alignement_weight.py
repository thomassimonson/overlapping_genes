import argparse
import math
import glob
import re
import seq_tools as st

# Description #################################################################
# This script contains functions for retrieve amino acids frequency at each
# position of an alignment. These frequency are used as weights for nucleotides
# modifications.
# Description #################################################################

residues = [
    "A", "C", "E", "D", "F", "I", "H", "K", "M", "L", "N", "Q", "S", "R", "T",
    "W", "V", "Y", "P", "G", "_"
]

def get_pfam_id(query):
    """Get pfam ids from header like this 'PF0001.11_FHLKKJHLH/10-1000' -> PF0001
    """
    result = re.match(r"^(\w+)\.?", query)
    return result.group(1)

def read_alignement(family, family_dir, groups):
    """read the alignement file for each family and return a list of
    dictionaries that contain {aa : frequency}. The groups variable contains a
    dictionary {aa : group}
    """
    flag = True
    nb_seq = 0
    align = find_alig_file(get_pfam_id(family), family_dir)
    for name, seq in st.read_fasta_file(align):
        nb_seq += 1
        if flag:
            # the first sequence is the one that have been used to represents
            # the family.
            fst = (name, seq)
            flag = False
            results = [{aa: 0 for aa in set(groups.values())} for i in range(len(seq))]
        if not flag:
            for pos, res in enumerate(seq):
                try:
                    results[pos][groups[res]] += 1
                except:
                    "not in"
    name, seq = fst
    final_results = [entropy(count, nb_seq) for ind, count in enumerate(results) if
                     seq[ind] in residues]
    return final_results

def find_alig_file(family, family_dir):
    """find the corresponding alignment file in FASTA format.
    """
    return glob.glob("{}/{}*.fa".format(family_dir, family))[0]

def get_group(infile="./data/dict_aa"):
    """get a dictionary of groups for amino-acids.
    """
    reg_exp = re.compile("\S+")
    groups = {}
    with open(infile) as group_file:
        for line in group_file:
            if not line.startswith("#"):
                val = reg_exp.findall(line)
                groups[val[0]] = val[1]
    return groups

def get_freq(count, nb_seq):
    "Compute frequencies"
    return float(count) / nb_seq

def entropy(count, nb_seq):
    "compute entropy"
    if get_freq(count['-'], nb_seq) > 0.3:
        return 1. / 7
    else:
        return 1 / math.exp(- sum([get_freq(cnt, nb_seq) * math.log(get_freq(cnt, nb_seq)) for cnt in count.values() if cnt > 0]))
