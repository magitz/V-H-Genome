#!/usr/bin/env python
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os, re, sys
import numpy as np


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
parser.add_argument("tsvFile", help="Filename of virushostdb.tsv file")
parser.add_argument("GenBankFile", help="Filename of GenBank formatted file to process")
parser.add_argument("outputFile", help="Filename for output")

args = parser.parse_args()


# Function based on https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/
def index_genbank_features(gb_record, feature_type, qualifier) :
    feature_dict= dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in feature_dict :
                        #print ("WARNING - Duplicate key %s for %s features %i and %i" \
                        #   % (value, feature_type, feature_dict[value], index))
                        continue
                    else :
                        feature_dict[value] = index
    return feature_dict


def parse_GenBank_file(file):

    Record_summary={}

    for seq_record in SeqIO.parse(file, "genbank"):
        CDS_count = len( index_genbank_features(seq_record, "CDS", "gene") )
        Record_summary[seq_record.id]= (len(seq_record.seq), CDS_count )

    return Record_summary


try:
	IN=open(args.tsvFile, 'r')
except:
    print( "Can't open tsv file for reading: %s" %(args.tsvFile) )
    exit()

try:
	OUT=open(args.outputFile, 'w')
except:
    print("Can't open output file for writing: %s" %(args.outputFile))
    exit()



print("Parsing %s for genome information. This may take a while. You will be notified when done.\n\n" %(args.GenBankFile))
Virus_data = parse_GenBank_file(args.GenBankFile)
print("Done. \nAdding genome size and gene count to data file.\n\n")



count = 0

for Line in IN:
    if count == 0:
        Line_bits=Line.strip.split('\t')
        Header = (Line_bits [0:3],'genome_size','Gene_count',Line_bits[4:])
        OUT.write(Header)
        count+=1
    
    else:
        Line_bits=Line.strip.split('\t')
        refseq_id = Line_bits[3]
        if len(refseq_id) == 1:
            # If there is only one refseq genome, use its data
            Outline = (Line_bits [0:3],Virus_data[Line_bits[3]][0],Virus_data[Line_bits[3]][1],Line_bits[3:])
            OUT.write(Outline)

        else:
            # If there are more than one refseq genomes, need to average the data.
            Genome_array= np.array()
            for genome in Line_bits[3]:
                Genome_array.np.append(Virus_data[Line_bits[3]])
            Outline = (Line_bits [0:3],np.mean(Genome_array, axis=1), Line_bits[3:])

