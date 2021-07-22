#! /usr/bin/python2.7

from collections import Counter
import csv
import os
import sys
import re
import numpy as np
from itertools import *
from sys import stdout
from numpy import genfromtxt
from Bio import AlignIO
import time
import argparse
import tempfile
from Bio import SeqIO
import math
import itertools
import jellyfish 


#from scipy import DataFrame
from multiprocessing import Lock, Process, Queue, current_process
#from progress_bar import ProgressBar
    #############################
 
def gen_sub(s, len_chunk):
   	for start in range(0, len(s)-len_chunk+1):
       		yield s[start:start+len_chunk]
   #############################
def get_distance(wc_seq1, wc_seq2):
	dd=0
	for kmer in wc_seq1:
		dd=dd+(wc_seq1[kmer]-wc_seq2[kmer])*(wc_seq1[kmer]-wc_seq2[kmer])
	return(dd)

#############################
#############################
parser = argparse.ArgumentParser(description='K-mers count calculation.')
parser.add_argument('-s', dest="seq_file", required=True, help='input fasta file')
parser.add_argument('-o', dest="output_file", required=True, help='output tabulated file')
parser.add_argument('-k', dest="k",   type=int, default=2, help='length of the words (default: 2)')
parser.add_argument('-f', dest="freq",   action='store_true', default=False, help='to print word frequencies instead of word counts (default: false)')
parser.add_argument('-p', dest="print_freq",   action='store_true', default=False, help='to print frequencies of word appareance among all sequences (default: false)')
parser.add_argument('-v', dest="verbose",   type=int, default=1, help='to be verbose (default: 1)')
parser.add_argument('-F', dest="minimal_freq",   type=float, default=0.0, help='the minimal frequency to select words for output. If set to 0.1, for example, only words observed in at least 10% of the sequences will be printed into the output file (default: 0.0).')



N=20

args = parser.parse_args()

if os.path.exists(args.output_file):
    os.remove(args.output_file)

handle1=open(args.seq_file, 'r')
seq_dict = SeqIO.to_dict(SeqIO.parse(handle1, "fasta"))
handle1.close()

if args.verbose > 0:
	print("sequences: ", len(seq_dict))




#kmers={}
#for ii in range(args.k1,args.k2+1): 
#	kmers[ii]=list(product('ACGT', repeat=ii))
#	print str(ii)+"-mers: ", len(kmers[ii])
#	
#	for jj in range(0, len(kmers[ii])):
#		kmers[ii][jj]=''.join(kmers[ii][jj])
	
all_words=[]
kmers_count={}
nb_kmers_count={}
for ss1 in seq_dict:
	kmers_count[ss1]=Counter([sub for sub in gen_sub(str(seq_dict[ss1].seq.upper()),args.k )])
	all_words= all_words+ kmers_count[ss1].keys()
	all_words=list(set(all_words))
	if args.freq:
		for ii in kmers_count[ss1]:
			kmers_count[ss1][ii]=float(kmers_count[ss1][ii])/(len(seq_dict[ss1].seq)-args.k+1.0)

#		print ss1, ii, kmers_count[ss1][ii]
		
#		print "#######################"
#		print ss1
#		print kmers_count[ss1]
#		print "#######################"
#	print all_words
		
#ss1='lcl|chr12_38292356_38301939_C:7635-7814'
#ss2='lcl|chr12_38292356_38301939_C:6265-6435'
#print kmers_count[ss1]-kmers_count[ss2]	
#xx=kmers_count[ss1]
#xx.subtract(kmers_count[ss2])
#print xx
#print "#######################"
all_words=sorted(all_words)

handle=open(args.output_file,'w')
handle.write("id")
nb_kmers_count= dict(zip(all_words , [0.0] * len(all_words)))
if args.minimal_freq > 0.0:
	if args.verbose > 0:
		print("initiale number of words: "+str(len(	all_words))	)

if args.minimal_freq > 0.0:
	for ss1 in seq_dict:
		for kkk in all_words:
			if kmers_count[ss1][kkk] > 0:
				nb_kmers_count[kkk]=nb_kmers_count[kkk]+1
	for kkk in all_words:
			nb_kmers_count[kkk]=nb_kmers_count[kkk]/len(seq_dict)
			if nb_kmers_count[kkk] < args.minimal_freq:
				if args.verbose > 0:
						print("removing: "+kkk+" with frequency: "+str(nb_kmers_count[kkk]))
				all_words=filter(lambda a: a != kkk, all_words)	
					
if args.minimal_freq > 0.0:
	if args.verbose > 0:
		print("final number of words: "+str(len(all_words))	)
for kkk in all_words:
	handle.write (" "+kkk)
handle.write ("\n")

for ss1 in seq_dict:
	handle.write(ss1)
	for kkk in all_words:
		if kkk in kmers_count[ss1]:
			handle.write (" "+str(kmers_count[ss1][kkk]))
		else:
			if args.freq:
				handle.write (" 0.0")
			else:
				handle.write (" 0")
	handle.write ("\n")
if args.print_freq:
		handle.write ("frequency")
		for kkk in all_words:
			handle.write (" "+str(nb_kmers_count[kkk]))
		handle.write ("\n")
handle.close()
if args.verbose > 0:
	print("Bye bye !!!")

exit()
