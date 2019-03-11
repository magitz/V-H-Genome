#!/usr/bin/env python

from Bio import Entrez
from Bio import SeqIO
import os, sys
import argparse
import time

# =====================================================
#  Takes a list of species names and queries GenBank for 
#  that species. If any data are in GenBank, a file is 
#  written that has the GenBank IDs for that species.
#
#  Matt Gitzendanner
#  University of Florida
#  3/07/16
#
# =====================================================


#####################
# Options
#
# -i input file with list of accessions
# -e email address used for Entrez
# -o Output file for GB Taxon IDs
#####################

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with accessions")
parser.add_argument("-e", help="email address")
parser.add_argument("-o", help="Output file for GB Taxon IDs")

args = parser.parse_args()

infile = args.i
Entrez.email = args.e #sets the email for Entrez.
OutFile= args.o


    
try:
	IN=open(infile, 'r')
except IOError:
    print ("Can't open file %s" %(infile))
    sys.exit()

try:
    OUT=open(OutFile, 'w')
except:
    print ("Can't open file: %s" %(OutFile))
    sys.exit()

count=0

AccessionDict={} # Setup dictionary to hold Accession: TaxonID info

for Line in IN:
    if count == 0:
        #skip header row
        count+=1
        continue
    else:
        count+=1
        Accession=Line.strip('\n')
        
        if Accession in AccessionDict:
            # We've already looked up this Accession, just write stored value
            OUT.write(AccessionDict[Accession] + "\n")
       
        else:
            for i in range(3, 0, -1):
                try:
                    GBSeq=Entrez.efetch(db="nucleotide", id=Accession, rettype="gb", retmode="text")
                except:
                    if i > 0:
                        print('Failed to connect. Retrying')
                        time.sleep(5) #Wait 5 seconds and try again.
                    else:
                        print("Can't download data for %s" %(Accession))
                        break

            for Sequence in SeqIO.parse(GBSeq, "gb"):
                print(Sequence.annotations['source'])
                for feature in Sequence.features:
                    if feature.type=='source':
                        try: 
                            OUT.write(feature.qualifiers.get("db_xref")[0])
                            OUT.write( "\n")
                            AccessionDict[Accession]=feature.qualifiers.get("db_xref")[0]

                        except:
                            for item in feature.qualifiers.get("db_xref"):
                                OUT.write("%s\t" %(item))
                            AccessionDict[Accession]=feature.qualifiers.get("db_xref")

