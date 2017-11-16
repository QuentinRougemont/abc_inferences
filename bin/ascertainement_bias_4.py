#!/usr/bin/env python
"""Filter loci to create an ascertainement bias

./ascertainement_bias.py num_wanted_loci block_size
"""

# Modules
from __future__ import print_function
from signal import signal, SIGPIPE, SIG_DFL
from random import random
import sys

# Fix broken pipe error
signal(SIGPIPE, SIG_DFL)

# Funcions
def filter_locus(infos):
    """Return True if pass, False if rejected
    """

    # Compute maf
    genotypes = infos[3:]
    num_genotypes = len(genotypes)

    allele_1 = len([x for x in genotypes if x == "0"])
    allele_2 = len([x for x in genotypes if x == "1"])
    minor_allele = min([allele_1, allele_2])
    maf = float(minor_allele) / float(num_genotypes)

    # Compute bias
    if maf == 0:
        bias = 0

    else:
        min_bias = 0.1
        bias = (((100. + 1 - 1 / (maf * 2)) / 100.) ** 8 + min_bias) / (1 + min_bias)

    passed = random() <+ bias
    return passed

def print_locus(infos):
    print()
    print("\n".join(infos))

# Parse user input
try:
    num_wanted_loci = int(sys.argv[1])
    block_size = int(sys.argv[2])
except:
    print(__doc__)
    sys.exit(1)

# Global variables
num_retained_loci = 0
num_treated_loci = 0

# Read input file
while True:
    line = sys.stdin.readline()
    if not line:
        break

    infos = []
    line1 = line.strip()

    #First line
    if not line.startswith("//"):
        continue 

    # Remember first line
    infos.append(line1)

    # Remember second line
    line2 = sys.stdin.readline().strip()
    assert line2.startswith("segsites")
    infos.append(line2)

    # Remember third line
    line3 = sys.stdin.readline().strip()
    assert line3.startswith("positions")
    infos.append(line3)
    
    # Retrieve genotype length
    length = int(line1.split()[1])
    genotype = []

    # Parse genotype
    for i in range(length):
        gene = sys.stdin.readline().strip()
        genotype.append(gene.strip())

    if num_retained_loci < num_wanted_loci and filter_locus(genotype):
        num_retained_loci += 1
        print_locus(infos + genotype)

    num_treated_loci += 1 

    #print(num_treated_loci)
    if num_treated_loci % block_size == 0:
        if num_retained_loci < num_wanted_loci:
            sys.stderr.write("Not enough input markers!!!\n")
            sys.exit(1)

        num_retained_loci = 0
        
sys.stderr.write("Ascertainement done\n")
