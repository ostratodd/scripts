#!/usr/bin/env python

#Created by Andrew Swafford. Developed and tested on Windows 7 in Python 2.7.*
#Contact Swafford.andrew@lifesci.ucsb.edu with questions/comments
import dendropy
import csv
from Bio import *
import argparse
import os, errno
from Bio import SeqIO

def inputs(seq_infile, tree_file):
    print("Organizing data..")
    records = SeqIO.index(seq_infile,'fasta')
    tree_tips = []
    treelist = dendropy.TreeList.get(path=tree_file, schema = 'newick',preserve_underscores=True)
    for tree in treelist:
        print('pruning tree..')
        for tip in tree.leaf_node_iter():
            tree_tips.append(tip.taxon.label)
    return (records,tree_tips)

def cut(records, tree_tips):
    print("Cutting alignment..")
    with open(aln_outfile,'w') as aln:
        for entry in records:
            print(records[entry].description)
            if records[entry].description in tree_tips:
                aln.write(">{0}\n{1}\n".format(records[entry].description,records[entry].seq))
    print("Done!")

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Cut alignment to match a tree")
    parser.add_argument('-p', metavar = 'prefix', type = str, required = True,
                        help = "Prefix to append onto output files")
    parser.add_argument('-o', metavar = 'output', type = str,
                        help = "Path to output dir.")
    parser.add_argument('-t', metavar = 'treeFile', type = str, required = True,
                        help = "Path to tree file (NEWICK).")
    parser.add_argument('-a', metavar = 'alnFile', type = str, required = True,
                        help = "Path to alignment file (FASTA).")
    args = parser.parse_args()
    if args.o == None:
        args.o = os.path.realpath(args.a)
        args.o = os.path.abspath(os.path.join(args.o,".."))
    args.o = os.path.normpath(args.o)
    print("Setting up output directory")
    mkdir_p(args.o)
    aln_outfile = os.path.join(args.o,args.p)
    if args.o == os.path.join(args.o,args.p):
        print(output_path)
    records,tips = inputs(args.a,args.t)
    cut(records,tips)

