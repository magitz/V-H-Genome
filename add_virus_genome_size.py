#!/usr/bin/env python
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os, re, sys
import time
import csv
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
        #CDS_count = len( index_genbank_features(seq_record, "CDS", "gene") )
        seq_name=re.split('\.', seq_record.id)[0]
        Record_summary[seq_name]= (len(seq_record.seq))

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


localtime = time.asctime( time.localtime(time.time()) )
print (localtime)

print("Parsing %s for genome information.\n  This may take a while.\n   You will be notified when done." %(args.GenBankFile))
Virus_data = parse_GenBank_file(args.GenBankFile)
print("   Done.\n")

localtime = time.asctime( time.localtime(time.time()) )
print (localtime)

print ("\nAdding genome size and gene count to data file.\n\n")



count = 0

for Line in IN:
    if count == 0:
        Line = Line.strip()
        Line_bits=re.split('\t', Line)

        for col in range(len(Line_bits)):
            if col == 3:
                Header=str(Line_bits[3]) + "\tGenome_size\tGene_count\t"
                OUT.write(Header)
            else:
                Header=str(Line_bits[col]) + "\t"
                OUT.write(Header)
        OUT.write("\n")
        count+=1

    else:
        Line = Line.strip()
        Line_bits=re.split('\t', Line)
        refseq_id = re.split(', ',Line_bits[3])

        if len(refseq_id) == 1:
            # If there is only one refseq genome, use its data
           for col in range(len(Line_bits)):
               if col == 3:
                   Header=str(Line_bits[3]) + "\t" + str(Virus_data[Line_bits[3]][0]) + "\t"
                        #str(Virus_data[Line_bits[3]][1]) + "\t"
                   OUT.write(Header)
               else:
                   Header=str(Line_bits[col]) + "\t"
                   OUT.write(Header)
           OUT.write("\n")

        else:
            # If there are more than one refseq genomes, need to average the data.
            Genome_size=[]
            #Gene_count=[]
            for genome in refseq_id:
                Genome_size.append(Virus_data[genome])
                #Gene_count.append(Virus_data[genome][1])

            for col in range(len(Line_bits)):
                if col == 3:
                    Header=str(Line_bits[col]) + "\t" + str(max(Genome_size)) + "\t"
                    OUT.write(Header)
                else:
                    Header=str(Line_bits[col]) + "\t"
                    OUT.write(Header)
            OUT.write("\n")
