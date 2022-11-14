from Bio import SeqIO
import khmer
from khmer import khmer_args
from oxli.functions import build_graph
import argparse

def get_options():
    description = 'Generates Khmer count graph in gfa format.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python generate_graph.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    required=True,
                    help='List of files paths to fasta files to generate gfa from, one per line. ')
    IO.add_argument('--out',
                    default="reference",
                    type=str,
                    help='Output prefix for gfa file. ')
    IO.add_argument('--kmer',
                    default=15,
                    type=int,
                    help='Kmer size (default=15). ')
    IO.add_argument('--maxtable',
                    default=1e5,
                    type=int,
                    help='Maximum khmer table size (default=1e6). ')
    IO.add_argument('--numtable',
                    default=4,
                    type=int,
                    help='Number of khmer tables (default=4). ')
    IO.add_argument('--threads',
                    default=1,
                    type=int,
                    help='Number of threads (default=1). ')

    return parser.parse_args()

def main():
    options = get_options()

    infile = options.infile
    threads = options.threads
    outpref = options.out
    kmer = options.kmer
    maxtable = options.maxtable
    numtable = options.numtable

    dataset = []

    fasta_sequences = SeqIO.parse(open(infile),'fasta')
    for fasta in fasta_sequences:
        dataset.append(str(fasta.seq))

    countgraph = khmer.Countgraph(kmer, maxtable, numtable)

    build_graph(dataset, countgraph, num_threads=threads)

    countgraph.save(outpref + ".khmer")

    countgraph.save(outpref + ".khmer")

