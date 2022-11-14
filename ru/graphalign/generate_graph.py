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
                    default=None,
                    type=int,
                    help='Kmer size (default=15). ')
    IO.add_argument('--maxtable',
                    default=None,
                    type=int,
                    help='Maximum khmer table size (default=None). ')
    IO.add_argument('--numtable',
                    default=None,
                    type=int,
                    help='Number of khmer tables (default=None). ')
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

    with open(infile, "r") as f:
        for line in f:
            dataset.append(line.strip())

    parser = khmer_args.build_counting_args("Aligns reads to count graph")

    args = parser.parse_args()

    if kmer is None:
        kmer = args.ksize
    if maxtable is None:
        maxtable = args.max_tablesize
    if numtable:
        numtable = args.n_tables

    countgraph = khmer.Countgraph(kmer, maxtable, numtable)

    build_graph(dataset, countgraph, num_threads=threads)

    countgraph.save(outpref + ".khmer")

