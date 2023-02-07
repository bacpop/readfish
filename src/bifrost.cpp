#include "bifrost.h"

CompactedDBG<> buildGraph (const std::string& infile_1,
                          const std::string& infile_2,
                          const bool is_ref,
                          const int kmer,
                          const int threads,
                          const std::string& output_prefix)
{
    std::ifstream infile1(infile_1);
    std::ifstream infile2(infile_2);
    CDBG_Build_opt opt;


    opt.k = kmer;
    opt.nb_threads = threads;
    opt.verbose = false;
    opt.prefixFilenameOut = output_prefix;

    std::string filename;
    if (is_ref && (infile_2 == "NA")) {
        while (std::getline(infile1, filename))
        {
            opt.filename_ref_in.push_back(filename);
        }
    } else if (!is_ref && (infile_2 == "NA"))
    {
        while (std::getline(infile1, filename))
        {
            opt.filename_seq_in.push_back(filename);
        }
    } else {
        while (std::getline(infile1, filename))
        {
            opt.filename_ref_in.push_back(filename);
        }
        while (std::getline(infile2, filename))
        {
            opt.filename_seq_in.push_back(filename);
        }
    }

    CompactedDBG<> cdbg(opt.k);
    cdbg.build(opt);
    cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);

    cdbg.write(opt.prefixFilenameOut, opt.nb_threads);

    return cdbg;
}

void Graph::build (const std::string& infile1,
                   const int kmer,
                   size_t num_threads,
                   bool is_ref,
                   const std::string& infile2,
                   const std::string& outpref) {
    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // read in compact DBG
    cout << "Building compacted DBG..." << endl;

    // initialise persistent variables
    int overlap = kmer - 1;

    if (infile2 != "NA") {
        is_ref = 0;
    }

    // generate graph
    _cdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, outpref);
    _kmer = kmer;

    cout << "Graph written to " << outpref << ".gfa" << endl;
}

// read existing graph and index
void Graph::read (const std::string& graphfile) {
    // read in graph
    _cdbg.read(graphfile);
    _kmer = _cdbg.getK();
}

double Graph::query (const std::string& query, const int gap) {

    // hold number of kmers
    const size_t num_kmers = query.length() - _kmer + 1;
    const size_t num_split_kmers = num_kmers - (_kmer + gap);
    int total_matches = 0;

    // convert query to string for search in graph
    const char *query_str = query.c_str();

    std::vector<bool> match_vec;

    for (KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km)
    {
        // count all matches
        auto um = _cdbg.find(it_km->first);
        match_vec.push_back(!um.isEmpty);
    }

    // iterate over matches, linking split-kmers
    int ind1 = 0;
    int ind2 = ind1 + _kmer + gap;

    // iterate over all split-kmers
    for (; ind2 < num_kmers; ind2++)
    {
        if (match_vec[ind1] && match_vec[ind2])
        {
            total_matches++;
        }
        ind1++;
    }

    double prop_match = (double)total_matches / (double)num_split_kmers;

    return prop_match;
}

