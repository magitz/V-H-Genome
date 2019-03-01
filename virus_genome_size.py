#!/usr/bin/env python
import argparse
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

# Function based on https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
def index_genbank_features(gb_record, feature_type, qualifier) :
    feature_count= dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print ("WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, feature_count[value], index))
                    else :
                        feature_count[value] = index
    return feature_count

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("file", help="Filename of GenBank formatte file to process")
args = parser.parse_args()


for seq_record in SeqIO.parse(args.file, "genbank"):
    CDS_count = index_genbank_features(seq_record, "CDS", "gene")
    print(seq_record.id, len(seq_record.seq), CDS_count )

