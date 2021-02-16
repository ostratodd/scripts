from Bio import SeqIO
import argparse

# Construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-p", "--path", required=False, default='.', 
        help="path to fasta file default is ./")
ap.add_argument("-f", "--file", required=True, type=str,
	help="file name for fasta file")
ap.add_argument("-m", "--min", required=False, default=100, type=int,
	help="delete sequences below this minimum length")
args = vars(ap.parse_args())
path = args["path"]
file = args["file"]
min = args["min"]

short_sequences = []  # Setup an empty list
long_sequences = []  # Setup an empty list
for record in SeqIO.parse(file, "fasta"):
    if len(record.seq) < min:
        # Add this record to our list
        short_sequences.append(record)
    else:
        long_sequences.append(record)

print("Found %i short sequences" % len(short_sequences))
shortfile = 'short' + str(min) + '_' + str(file)
SeqIO.write(short_sequences, shortfile, "fasta")

print("Found %i long sequences" % len(long_sequences))
longfile = 'long' + str(min) + '_' + str(file)
SeqIO.write(long_sequences, longfile, "fasta")
