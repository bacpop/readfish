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

PathVec iter_nodes(CompactedDBG<>& cdbg,
                   const Kmer& head_kmer,
                   const bool strand,
                   const size_t length,
                   const size_t overlap)
{
    PathVec complete_paths;

    std::stack<PathTuple> node_stack;
    std::vector<std::pair<Kmer, bool>> node_vector;

    // create stack
    node_stack.push({0, head_kmer, strand, 0});

    while(!node_stack.empty())
    {
        // pop node in stack
        auto node_tuple = node_stack.top();
        node_stack.pop();

        // unpack tuple
        const size_t & pos_idx = std::get<0>(node_tuple);
        const Kmer & node_head = std::get<1>(node_tuple);
        const bool & node_strand = std::get<2>(node_tuple);
        const size_t & path_length = std::get<3>(node_tuple);

        // slice path, unless at first node
        if (pos_idx != 0)
        {
            node_vector = std::vector<std::pair<Kmer, bool>> (node_vector.begin(), node_vector.begin() + pos_idx);
        }

        // add node to path
        node_vector.push_back({node_head, node_strand});

        // get length of vector for new pos_idx
        const size_t new_pos_idx = node_vector.size();

        auto um = cdbg.find(node_head, true);

        um.strand = node_strand;

        // iterate over neighbours until length satistified
        for (auto& neighbour_um : um.getSuccessors())
        {
            Kmer neighbour_head = neighbour_um.getUnitigHead();
            bool neighbour_strand = neighbour_um.strand;

            const size_t updated_path_length = path_length + (neighbour_um.size - overlap);


            // if at full length, return
            if (updated_path_length >= length)
            {
                std::vector<std::pair<Kmer, bool>> return_path = node_vector;
                return_path.push_back({neighbour_head, neighbour_strand});

                // add to final list and move onto next node
                complete_paths.push_back(std::move(return_path));
                continue;
            }

            node_stack.push({new_pos_idx, neighbour_head, neighbour_strand, updated_path_length});
        }
    }

    return complete_paths;
}

// generate strings and split k-mers from complete paths
KmerMap generate_split_kmers(CompactedDBG<>& cdbg,
                              const PathVec& complete_paths,
                              const size_t length,
                              const size_t overlap,
                              const bool reverse,
                              const bool forward_empty)
{
    KmerMap kmer_map;

    std::vector<uint64_t> kmer_vec;

    // generate initial string for start node
    {
        std::vector<uint64_t> start_kmer_vec;

        const auto& start_node = complete_paths[0][0];
        auto um = cdbg.find(start_node.first, true);
        std::string start_str;

        if (start_node.second)
        {
            start_str = um.referenceUnitigToString();
        } else
        {
            start_str = reverse_complement(um.referenceUnitigToString());
        }

        // determine number of kmers
        const size_t num_kmers = start_str.size() - overlap;

        size_t kmer_count = 0;

        const char *um_cstr = start_str.c_str();
        for (KmerIterator it_km(um_cstr), it_km_end; it_km != it_km_end; ++it_km)
        {
            // generate vector of all kmer strings
            auto kmer_str = it_km->first.toString();
            start_kmer_vec.push_back(to_binary(kmer_str));

            // add to kmer_vec if within range of end of node
            int diff = num_kmers - kmer_count;
            if (diff <= length)
            {
                kmer_vec.push_back(to_binary(kmer_str));
            }

            kmer_count++;
        }

        // iterate and pair kmers that are less than kmer + gap away from end
        // only do in forward direction to avoid duplicate effort, or if forward is empty
        if (!reverse || forward_empty)
        {
            int ind1 = 0;
            int ind2 = ind1 + length;

            // iterate over all split-kmers
            for (; ind2 < num_kmers; ind2++)
            {
                // add forward match...
                kmer_map[start_kmer_vec[ind1]].insert(start_kmer_vec[ind2]);

                // ...and reverse complement
                kmer_map[ReverseComp64(start_kmer_vec[ind2], overlap + 1)].insert(ReverseComp64(start_kmer_vec[ind1], overlap + 1));

                ind1++;
            }
        }
    }

    // now iterate through extended paths, only add in forward direction as will be covered by other unitig traversal
    const int stop_idx = kmer_vec.size() - 1;
    {
        for (const auto& path : complete_paths)
        {
            bool complete = false;

            int ind1 = 0;

            size_t kmer_count = stop_idx;

            for (int nidx = 1; nidx < path.size(); nidx++)
            {
                const auto& start_node = path[nidx];
                auto um = cdbg.find(start_node.first, true);

                // determine position of kmers in path
                const size_t pos = (um.size - overlap) + kmer_count;

                // if unitig too short for gap, ignore and move on
                if (pos < length)
                {
                    kmer_count += (um.size - overlap);
                    continue;
                }

                std::string seq;
                if (start_node.second)
                {
                    seq = um.referenceUnitigToString();
                } else
                {
                    seq = reverse_complement(um.referenceUnitigToString());
                }

                // get the sequence and iterate through kmers until satisified
                const char *seq_cstr = seq.c_str();

                int ind2 = (length + ind1) - (kmer_count + 1);

                // start at matched kmer
                KmerIterator it_km(seq_cstr),it_km_end;
                it_km += ind2;
                for (; it_km != it_km_end; ++it_km)
                {
                    // generate vector of all kmer strings
                    auto kmer_str = it_km->first.toString();

                    // add kmers until all kmers in original unitig can be paired
                    kmer_map[kmer_vec[ind1]].insert(to_binary(kmer_str));

                    ind1++;
                    // all kmers paired
                    if (ind1 > stop_idx)
                    {
                        complete = true;
                        break;
                    }
                }

                if (complete){
                    break;
                }

                kmer_count += (um.size - overlap);
            }
        }
    }

    return kmer_map;
}

KmerMap traverse_unitig(CompactedDBG<>& cdbg,
                        const Kmer& head_kmer,
                        const size_t length,
                        const size_t overlap)
{
    KmerMap kmer_map;
    bool forward_empty = false;

    // iterate forward
    {
        auto complete_paths = iter_nodes(cdbg, head_kmer, true, length, overlap);
        if (!complete_paths.empty())
        {
            kmer_map = generate_split_kmers(cdbg, complete_paths, length, overlap, false, forward_empty);
        } else
        {
            forward_empty = true;
        }
    }

    // iterate reverse
    {
        auto complete_paths = iter_nodes(cdbg, head_kmer, false, length, overlap);
        if (!complete_paths.empty())
        {
            auto kmer_map_temp = generate_split_kmers(cdbg, complete_paths, length, overlap, false, forward_empty);
            kmer_map.insert(make_move_iterator(kmer_map_temp.begin()), make_move_iterator(kmer_map_temp.end()));
        }
    }

    return kmer_map;
}

void Graph::index_split_kmer()
{
    // initialise return values
    KmerVec head_kmer_arr;

    // determine length
    const size_t length = _gap + _kmer;
    const size_t overlap = _kmer - 1;

    // iterate over each unitig to get head kmer array.
    for (auto& um : _cdbg)
    {
        head_kmer_arr.push_back(um.getUnitigHead());
    }

    // set up progress bar
    progressbar bar(head_kmer_arr.size());
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("|");
    bar.set_closing_bracket_char("|");

    // reiterate, use openMP
    #pragma omp parallel
    {
        KmerMap kmer_map_temp;
        #pragma omp for nowait
        for (auto it = head_kmer_arr.begin(); it < head_kmer_arr.end(); it++)
        {
            auto kmer_map_temp_temp = traverse_unitig(_cdbg, *it, length, overlap);
            kmer_map_temp.insert(make_move_iterator(kmer_map_temp_temp.begin()), make_move_iterator(kmer_map_temp_temp.end()));
            // update progress bar
            #pragma omp critical
            {
                bar.update();
                _kmermap.insert(make_move_iterator(kmer_map_temp.begin()), make_move_iterator(kmer_map_temp.end()));
            }
        }
    }

    // add new line to account for progress bar
    cout << endl;
}

void Graph::build (const std::string& infile1,
                   const size_t kmer,
                   const size_t gap,
                   size_t num_threads,
                   bool is_ref,
                   const std::string& infile2,
                   const std::string& outpref,
                   bool splitk) {
    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // set OMP number of threads
    omp_set_num_threads(num_threads);

    // read in compact DBG
    cout << "Building compacted DBG..." << endl;

    if (infile2 != "NA") {
        is_ref = 0;
    }

    // generate graph
    _cdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, outpref);
    _kmer = kmer;
    _gap = gap;

    cout << "Graph written to " << outpref << ".gfa" << endl;

    // generate kmermap
    if (splitk)
    {
        cout << "Generating split-kmer index..." << endl;
        index_split_kmer();

        // write sk_index

        // generate map that can be serialised
        StdKmerMap kmermapstd;
        for (const auto& i1 : _kmermap)
        {
            for (const auto& i2 : i1.second){
                kmermapstd[i1.first].insert(i2);
            }
        }

        std::tuple<int, int, StdKmerMap> for_writing = {_gap, _kmer, kmermapstd};

        std::ofstream os(outpref + ".sk", std::ios::binary);
        cereal::BinaryOutputArchive oarchive(os);
        // write class instance to archive
        oarchive(for_writing);

        cout << "Split-kmer index written to " << outpref << ".sk" << endl;
    }
}

void Graph::read_skindex(const std::string& sk_index)
{
    std::tuple<int, int, StdKmerMap> for_writing;

    std::ifstream is(sk_index, std::ios::binary);
    cereal::BinaryInputArchive iarchive(is);
    iarchive(for_writing);

    // regenenerate _kmermap
    for (const auto& i1 : std::get<2>(for_writing))
    {
        for (const auto& i2 : i1.second){
            _kmermap[i1.first].insert(i2);
        }
    }


    _gap = std::get<0>(for_writing);
    _gapd = (double)_gap;
    _kmer = std::get<1>(for_writing);
    _kmerd = (double)_kmer;

    _param1 = (-1/((2 * _kmerd) + _gapd));

    _cdbg = CompactedDBG<> (_kmer);
}

// read existing graph and index
void Graph::read (const std::string& infile)
{
    if (infile.find(".sk") != string::npos)
    {
        read_skindex(infile);
        _splitk = true;
    } else {
        _cdbg.read(infile);
        _kmer = _cdbg.getK();
        _kmerd = (double) _kmer;
        _param1 = (-1/_kmerd);
        _splitk = false;
    }
}

std::pair<double, double> Graph::query (const std::string& query) {

    // hold number of kmers
    int total_matches = 0;
    int total_mismatches = 0;

    const size_t num_kmers = query.size() - _kmer + 1;

    // convert query to string for search in graph
    const char *query_str = query.c_str();

    std::vector<uint64_t> kmer_vec;

    for (KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km)
    {
        // just count kmers
        if (!_splitk)
        {
            auto um = _cdbg.find(it_km->first);

            if (!um.isEmpty)
            {
                total_matches++;
            } else
            {
                total_mismatches += 2;
            }
        } else
        {
            // count all matches
            auto kmer_str = it_km->first.toString();
            kmer_vec.push_back(to_binary(kmer_str));
        }
    }

    // count split-kmers
    if (_splitk)
    {
        // iterate over matches, linking split-kmers
        int ind1 = 0;
        int ind2 = ind1 + _kmer + _gap;

        // iterate over all split-kmers
        for (; ind2 < num_kmers; ind2++)
        {
            // find kmer1 in _kmermap
            auto kmer1_found = _kmermap.find(kmer_vec[ind1]);
            if (kmer1_found != _kmermap.end())
            {
                const auto& kmer_set = kmer1_found->second;
                if (kmer_set.find(kmer_vec[ind2]) != kmer_set.end())
                {
                    total_matches++;
                } else
                {
                    total_mismatches += 2;
                }
            } else
            {
                total_mismatches += 2;
            }

            ind1++;
        }
    }

    // calcaulte jaccard index
    double jaccard = (double)total_matches / (double)(total_matches + total_mismatches);

    // calculate mash distance, https://mash.readthedocs.io/en/latest/distances.html and https://www.biorxiv.org/content/10.1101/453142v1.full.pdf
    double mash_sim = 0;
    if (jaccard > 0)
    {
        double param2 = log((2 * jaccard) / (1 + jaccard));
        mash_sim = 1 - (_param1 * param2);
    }

    double prop_kmers = (double)total_matches / (double)num_kmers;

    std::pair<double, double> out_pair = {mash_sim, prop_kmers};

    return out_pair;
}

// SKA2 functions
uint64_t to_binary(std::string& current_kmer) {
    // convert k-mer to bitvector
    uint64_t packed_int = 0;
    for (auto it = current_kmer.cbegin(); it != current_kmer.cend(); ++it) {
        packed_int = packed_int << 2;
        packed_int += look_up_table[*(it)];
    }
    std::bitset<64> x(packed_int);
    return packed_int;
}

uint64_t ReverseComp64(const uint64_t mer, uint8_t kmerSize)
{
    uint64_t res = ~mer;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);

    return(res >> (2ULL * (32 - kmerSize)));
}

