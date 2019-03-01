#!/usr/bin/env python
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os


# =====================================================
#  
#  Matt Gitzendanner
#  University of Florida
#  03/01/19
#
# Version 1.0: Initial version 
# =====================================================

import argparse
parser = argparse.ArgumentParser()
paerser.add_argument("file", help="Filename of GenBank formatte file to process")
parser.parse_args()


for seq_record in SeqIO.parse(args.file, "genbank"):
    print(seq_record.id, len(seq_record.seq) )
