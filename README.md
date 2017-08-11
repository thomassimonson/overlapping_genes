# overlapping_genes
Software to design overlappng gene sequences

            COMPUTATIONAL DESIGN OF FULLY OVERLAPPING CODING
                 SCHEMES FOR PROTEIN PAIRS AND TRIPLETS
           ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━


About
═════

  This repository contains an implementation of the algorithm presented
  in "Computational design of fully overlapping coding schemes for
  protein pairs and triplets" (Vaitea Opuu, Martin Silvert, Thomas
  Simonson).

  For a given set of sequences D (in fasta format), this program will
  perform D × D optimization for pairs or D × D × D for triplets. For
  all pair of protein sequences (x, y), the program will produce a pair
  (x', y') where x' and y' are overlapping.


Usage
═════

Generate overlapping pairs
──────────────────────────

  ┌────
  │ python overlapping.py -i test_data/test.fa -f -2 -v t
  └────
  • `-f' : The reading frame
  • `-v' : Pretty output
  • `-i' : Input fasta file
  • `-o' : Write sequences in an output fasta file
  • `-e' : set it to "y" or anything else in order to consider the
    alignment otherwise don't use this option.
  • `-h' : short help

  This is an example of output:
  ┌────
  │ ==============================================================
  │ SCORE : 639 / FRAME : -2
  │ PF00636.25_RNC_UREPA/43-140 => similarity X, X' = 309.0
  │ X : QRLEFLGDACVEWVISNFIFNYKIKDNEKMRSLDEGEMTRARSNMVRSEILSYAAKDL...
  │ X': RRLEFLGDQCMQWIYQNIIFPYELKDPEKVRSLDEGNLTRPRSSMHRTELLTYASKTL...
  │ PF00636.25_RNC_UREPA/43-140 => similarity Y, Y' = 330.0
  │ Y : QRLEFLGDACVEWVISNFIFNYKIKDNEKMRSLDEGEMTRARSNMVRSEILSYAAKDL...
  │ Y': RRLEFLGDLCIEWILSNIIFPYELTDPEKVRSLDEGSLSRARSSMRRSEVLSYASKTL...
  │ 
  │   ------------------------------------------------------------
  │    97 96 95 94 93 92 91 90 89 88 87 86 85 84 83 82 81 80 79 78
  │ Y   G  Q  D  Q  A  V  A  G  I  F  A  E  F  I  D  E  Y  I  K  E
  │ Y'  G  D  Q  T  G  L  L  G  I  Y  A  T  S  I  D  S  Y  L  K  E
  │   GCGGCAGAACTCAAGGATCCTCTGGTTACATACGTCACCTATATAGTCTTATATTAAAAA
  │   CGCCGTCTTGAGTTCCTAGGAGACCAATGTATGCAGTGGATATATCAGAATATAATTTTT
  │ X' R  R  L  E  F  L  G  D  Q  C  M  Q  W  I  Y  Q  N  I  I  F
  │ X  Q  R  L  E  F  L  G  D  A  C  V  E  W  V  I  S  N  F  I  F
  │    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
  │ ...
  └────


Generate overlapping triplets
─────────────────────────────

  ┌────
  │ python overlapping.py -i test_data/test.fa -fl 1 2 4 -v t
  └────

  This is an example of output
  ┌────
  │ ==============================================================
  │ SCORE : 664.0 / CADRE : 1,2,4
  │ PF00636.25_RNC_UREPA/43-140 => similarity X, X' = 276.0, frame position = 1
  │ X : QRLEFLGDACVEWVISNFIFNYKIKDNEKMRSLDEGEMTRARSNMVRSEILSYAAKDL...
  │ X': QRLEVVSDSCVHYISSNYHFSYVLSDKSSLLSLDEGERSRARSNNLQEQILRYASKDL...
  │ PF04297.13_Y142_UREPA/4-99 => similarity Y, Y' = 202.0, frame position = 2
  │ Y : KNKRWYLIALYDIYQGLLTTKQCEYFNLHYFKDLSFSEIAELKEISKSAISDCLNKVC...
  │ Y': NVSKWFLIAVFIIYQVIITSLTCSLINLHYFLSTKESEAERVQTISKSRFSDTLQKIC...
  │ PF04297.13_Y142_UREPA/4-99 => similarity Z, Z' = 186.0, frame position = 4
  │ Z : KNKRWYLIALYDIYQGLLTTKQCEYFNLHYFKDLSFSEIAELKEISKSAISDCLNKVC...
  │ Z': KISQWYQQTLRRYYRILLTSSQIEFFQYHYLRVRSQSVADLLKRIGESALGDCLNALC...
  │ 
  │   ------------------------------------------------------------
  │ Z   I  D  K  L  K  K  V  L  E  S  D  N  I  L  T  Y  L  D  N  R
  │ Z'  V  D  R  L  P  K  Q  Y  S  H  E  N  Y  I  L  Y  N  D  S  R
  │   GTTGCAGAGCTTCACCAAAGACTATCGACACAAGTAATATATAGTTCATTAATAGTGAAG
  │   CAACGTCTCGAAGTGGTTTCTGATAGCTGTGTTCATTATATATCAAGTAATTATCACTTC
  │ Y'  N  V  S  K  W  F  L  I  A  V  F  I  I  Y  Q  V  I  I  T  S
  │ Y   K  N  K  R  W  Y  L  I  A  L  Y  D  I  Y  Q  G  L  L  T  T
  │ X' Q  R  L  E  V  V  S  D  S  C  V  H  Y  I  S  S  N  Y  H  F
  │ X  Q  R  L  E  F  L  G  D  A  C  V  E  W  V  I  S  N  F  I  F
  │ ...
  └────


Documentations
══════════════

Requirement
───────────

  Python >= 2.7


Files and directories
─────────────────────

  • `Alignment_weight.py' : This module contains functions that retrieve
    the composition of amino-acids in order to apply a weight during the
    optimization . For each sequence in the input fasta file, you need
    to have the alignment file in the family directory
    (`data/pfam_family/' by default).

  • `overlapping.py' : This module performs the pairs optimization.

  • `overlapping_g.py' : This module is the generalization of the pairs
    optimization for triple genes.

  • `seq_tools.py' : This module contains a collection of functions used
    in the other modules.

  • `data/' : This directory contains all the data that is needed:
    • `BLOSUM62' : is the substitution matrix
    • `dict_aa' : is the file that make correspond amino-acids to groups
      (D, N, Q, E are in the same group for instance)
    • `gencode.dat' : contains a dictionary of codon: amino-acids
    • `pfam_family/' : contains the alignment files from the *PFAM*
      `seed' data base.


Details
═══════

Reading frame for pair optimization
───────────────────────────────────

  For the pair optimization, the reading frame (`-f' option) could take
  any value in [-2, -1, 0, 1, 2]. The output sequences x' and y' can be
  appended to a file (`-o' option).

  • x' fasta name in the output file is >x_Vs_y_X
  • y' fasta name in the output file is >x_Vs_y_Y


Reading frame for triplet optimization
──────────────────────────────────────

  For the triplets optimization we need to set an absolute position for
  each codon. For a given quartet, we can encode 4 different codons. The
  position of a codon regarding the others determines the reading frame.
  For instance, here is a quartet (and its complement):

  ┌────
  │ TAGC
  │ ATCG
  └────

  We can decompose this quartet in 4 different codons:
  • 1 = ATC / sens
  • 2 = TCG / sens
  • 3 = GAT / anti-sens
  • 4 = CGA / anti-sens

  Thus, if you choose the `frame_list' [1, 2, 4] (`-fl 1 2 4'), the
  sequence X (first element of the list) will be in the reading frame 1
  regarding the Y (second sequence) and the X will be in the reading
  frame -2 regarding Z (third sequence).

  • x' fasta name in the output file is >x_Vs_y_Vs_z_X
  • y' fasta name in the output file is >x_Vs_y_Vs_z_Y
  • z' fasta name in the output file is >x_Vs_y_Vs_z_Z


References
══════════

  • Computational design of fully overlapping coding schemes for protein
    pairs and triplets (Vaitea Opuu, Martin Silvert, Thomas Simonson)
