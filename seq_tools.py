#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import argparse
from random import randint

# Description #################################################################
# This script contains functions for sequences analysis or creation
# Description #################################################################

# Global variables ############################################################
global GENCODE, COMPLEMENT
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

# Functions ###################################################################
def read_fasta_file(fasta_file):
    """read an alignement file and return a generator of tuple that represent the
    name and the sequence of a protein for the same family.
    """
    reg_exp = re.compile(r">(\S+)")
    seq_list = []
    name = ""
    seq = ""
    with open(fasta_file) as fasta:
        for line in fasta:
            if line.startswith(">"):
                if len(seq) > 0:
                    taille = len(seq)
                    yield (name, seq)
                    seq = ""
                name = reg_exp.match(line).group(1)
            else:
                seq += line.strip("\n\r")
        if len(seq) > 0:
            yield (name, seq)

def remove_invalid_aa(seq, inv_gencode):
    """this function return a sequence free of gaps and not normal amino acids.
    """
    return [aa for aa in seq if aa in inv_gencode.keys()]

def read_gencode(infile="./data/gencode.dat"):
    """Read the genetic code from a file.
    """
    regexp = re.compile(r"\S+")
    gencode = []
    inv_gencode = []
    with open(infile) as gencode_file:
        for line in gencode_file:
            val = regexp.findall(line)
            dna, prot = val[0], val[1]
            gencode.append((dna, prot))
            inv_gencode.append((prot, dna))
    return gencode

def write_fasta_format(seq, seq_len=60):
    """return a 60 char stk_line string
    """
    if len(seq) < seq_len:
        return seq

    temp = ""
    i = 0
    while i < len(seq):
        temp += seq[i]
        i += 1
        if i % seq_len is 0:
            temp += "\n"
    return temp

def read_score(score_mat):
    """get score from matrix
    """
    header = True
    head = []
    scores = {}
    regexp = re.compile(r"\S+")
    with open(score_mat) as mat:
        for line in mat:
            if not line.startswith("#"):
                if header:
                    head = regexp.findall(line)
                    header = False
                else:
                    val = regexp.findall(line)
                    scores[val[0]] = {}
                    for i, aa_name in enumerate(head):
                        scores[val[0]][aa_name] = float(val[i+1])
    return scores

def complement(codon):
    """Return the complementary nucleotide for each elem of codon
    """
    return "".join([COMPLEMENT[elem] for elem in codon])

def split(sequence, n=3):
    """Split sequence into bloc of k elem
    """
    size, rest = divmod(len(sequence), n)
    return [sequence[k*3: (k+1)*3] for k in range(size)]

def print_sequence(seq_x, seq_y, seq_x_p, seq_y_p, dna, FRAME):
    """Pretty print sequences
    """
    def split(seq, n):
        "split"
        results = []
        i = 0
        while i + n < len(seq):
            j = i+n
            results.append(seq[i:j])
            i = j
        results.append(seq[i:])
        return results

    n = 20
    min_size = len(seq_x) if len(seq_x) < len(seq_y) else len(seq_y)
    if FRAME in [-2, 1]:
        format_x = lambda x: " "+"  ".join(el for el in x)
        format_y = lambda y: "  "+"  ".join(el for el in y)
        format_x_i = lambda i_x: "   "+" ".join(map(lambda i: "%-2d"%i, i_x))
        format_y_i = lambda i_y: "   "+" ".join(map(lambda i: "%2d"%i, i_y))
    elif FRAME in [-1, 2]:
        format_x = lambda x: "  "+"  ".join(el for el in x)
        format_y = lambda y: " "+"  ".join(el for el in y)
        format_x_i = lambda i_x: "   "+" ".join(map(lambda i: "%2d"%i, i_x))
        format_y_i = lambda i_y: "   "+" ".join(map(lambda i: "%-2d"%i, i_y))
    else:
        format_x = lambda x: " "+"  ".join(el for el in x)
        format_y = lambda y: " "+"  ".join(el for el in y)
        format_x_i = lambda i_x: "   "+" ".join(map(lambda i: "%-2d"%i, i_x))
        format_y_i = lambda i_y: "  "+" ".join(map(lambda i: "%2d"%i, i_y))
    ind_x = range(len(seq_x))
    if FRAME <= 0:
        ind_y = range(len(seq_y))[::-1]
    else:
        ind_y = range(len(seq_y))

    ind_x = [format_x_i(el) for el in split(ind_x, n)]
    seq_x = [format_x(el) for el in split(seq_x, n)]

    ind_y = [format_y_i(el) for el in split(ind_y, n)]
    if FRAME <= 0:
        seq_y = [format_y(el) for el in split(seq_y[::-1], n)]
    else:
        seq_y = [format_y(el) for el in split(seq_y, n)]

    seq_x_p = [format_x(el) for el in split(seq_x_p, n)]
    if FRAME <= 0:
        seq_y_p = [format_y(el) for el in split(seq_y_p[::-1], n)]
    else:
        seq_y_p = [format_y(el) for el in split(seq_y_p, n)]

    dna_x = split(dna, n*3)
    dna_y = split(complement(dna), n*3)

    for s_x, s_x_p, d_x, d_y, s_y_p, s_y, i_x, i_y in zip(seq_x, seq_x_p, dna_x, dna_y, seq_y_p, seq_y, ind_x, ind_y):
        print "  " + "-" * 60
        print "\n".join([i_y, "Y "+s_y, "Y'"+s_y_p, "  "+d_y, "  "+d_x, "X'"+s_x_p, "X "+s_x, i_x])

def print_sequence_3(seq_x, seq_y, seq_z, seq_x_p, seq_y_p, seq_z_p, dna, FRAME_LIST):
    """Pretty print sequences
    """

    def split(seq, n):
        "split"
        results = []
        i = 0
        while i + n < len(seq):
            j = i + n
            results.append(seq[i:j])
            i = j
        results.append(seq[i:])
        return results

    n = 20
    if FRAME_LIST[0] not in [2, 4]:
        format_x = lambda x: " " + "  ".join(el for el in x)
    else:
        format_x = lambda x: "  " + "  ".join(el for el in x)

    if FRAME_LIST[1] not in [2, 4]:
        format_y = lambda y: " " + "  ".join(el for el in y)
    else:
        format_y = lambda y: "  " + "  ".join(el for el in y)

    if FRAME_LIST[2] not in [2, 4]:
        format_z = lambda z: " " + "  ".join(el for el in z)
    else:
        format_z = lambda z: "  " + "  ".join(el for el in z)

    if FRAME_LIST[0] not in [3, 4]:
        seq_x = [format_x(el) for el in split(seq_x, n)]
        seq_x_p = [format_x(el) for el in split(seq_x_p, n)]
    else:
        seq_x = [format_x(el) for el in split(seq_x[::-1], n)]
        seq_x_p = [format_x(el) for el in split(seq_x_p[::-1], n)]

    if FRAME_LIST[1] not in [3, 4]:
        seq_y = [format_y(el) for el in split(seq_y, n)]
        seq_y_p = [format_y(el) for el in split(seq_y_p, n)]
    else:
        seq_y = [format_y(el) for el in split(seq_y[::-1], n)]
        seq_y_p = [format_y(el) for el in split(seq_y_p[::-1], n)]


    if FRAME_LIST[2] not in [3, 4]:
        seq_z = [format_z(el) for el in split(seq_z, n)]
        seq_z_p = [format_z(el) for el in split(seq_z_p, n)]
    else:
        seq_z = [format_z(el) for el in split(seq_z[::-1], n)]
        seq_z_p = [format_z(el) for el in split(seq_z_p[::-1], n)]

    dna_x = split(dna, n * 3)
    dna_y = split(complement(dna), n * 3)

    for s_x, s_x_p, s_y, s_y_p, d_x, d_y, s_z_p, s_z in zip(
            seq_x, seq_x_p, seq_y, seq_y_p, dna_x, dna_y, seq_z_p, seq_z):

        print "  " + "-" * 60
        print "\n".join([
            "Z " + s_z, "Z'" + s_z_p, "  " + d_y, "  " + d_x, "Y'" + s_y_p,
            "Y " + s_y, "X'" + s_x_p, "X " + s_x
        ])

