#!/usr/bin/env python

from Bio import Entrez
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

for Line in IN:
    if count == 0:
        #skip header row
        continue
    else:
        Accession=Line.strip('\n')

        for i in range(3, 0, -1):
            try:
                GBSeq = Entrez.esearch(db="genome", term=Accession ) #Get the sequence
            except:
                if i > 0:
                    print('Failed to connect. Retrying')
                    time.sleep(5) #Wait 5 seconds and try again.
                else:
                    print("Can't download data for %s" %(Organism))
                    break

            Record= Entrez.read(GBSeq)
        
            if int(Record["Count"]) > 0:
                print ("%s had %d records in GenBank" %(Accession, int(Record["Count"])))
            for id in Record[IdList]:
                OUT.write(Accession + "\t" + Record.taxonomy + "\n")
