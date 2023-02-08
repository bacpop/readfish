import argparse
import query_cpp
import sys

def get_options():
    description = 'Generates Bifrost graph in gfa format.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='ru_generate_graph')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--refs',
                    default=None,
                    help='List of files paths to assembly fasta files to generate gfa from, one per line. ')
    IO.add_argument('--reads',
                    default=None,
                    help='List of files paths to reads fasta files to generate gfa from, one per line. ')
    IO.add_argument('--out',
                    default="reference",
                    type=str,
                    help='Output prefix for gfa file. ')
    IO.add_argument('--kmer',
                    default=11,
                    type=int,
                    help='Kmer size (default=11). ')
    IO.add_argument('--gap',
                    default=1,
                    type=int,
                    help='Gap size for split-kmers (default=1). ')
    IO.add_argument('--threads',
                    default=1,
                    type=int,
                    help='Number of threads (default=1). ')

    return parser.parse_args()

def main():
    options = get_options()

    refs = options.refs
    reads = options.refs
    threads = options.threads
    out = options.out
    kmer = options.kmer
    gap = options.gap

    graph = query_cpp.Graph()

    # if refs file specified for building
    if (refs is not None) and (reads is None):
        graph.build(refs, kmer, threads, True, "NA", out)
    # if reads file specified for building
    elif (refs is None) and (reads is not None):
        graph.build(reads, kmer, threads, False, "NA", out)
    # if both reads and refs file specified for building
    elif (refs is not None) and (reads is not None):
        graph.build(refs, kmer, gap, threads, False, reads, out)
    else:
        print("Error: incorrect number of input files specified. Please only specify the below combinations:\n"
              "- List of assembly files to '--refs'.\n"
              "- List of read files to '--reads'.\n"
              "- A list of reference files and a list of read files to '--refs' and '--reads' respectively.")
        sys.exit(1)

    sys.exit(0)

if __name__ == "__main__":
    main()