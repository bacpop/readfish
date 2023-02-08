#ifndef BIFROST_H
#define BIFROST_H

#include <string>
#include <iostream>
#include <sstream>
#include <bifrost/CompactedDBG.hpp>
#include <vector>
#include <math.h>
#include <unordered_set>
#include <unordered_map>
#include <stack>

// openmp
#include <omp.h>

// boost serialising headers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/dynamic_bitset/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include "serialize_tuple.h"

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// global variable declaration
namespace py = pybind11;

typedef std::vector<Kmer> KmerVec;
typedef std::unordered_map<std::string, std::unordered_set<std::string>> KmerMap;

// vector of paths, which contain the head kmer and the strand
typedef std::vector<std::vector<std::pair<Kmer, bool>>> PathVec;

// tuple containing information for path (1st Path index, 2nd head Kmer, 3rd strand, 4th path length)
typedef std::tuple<int, Kmer, bool, size_t> PathTuple;

class Graph {
    public:
//    // constructors for debugging
//    // default constructor
//    Graph()
//    {
//        cout << "Graph default constructor" << endl;
//    };
//
//    // copy
//    Graph(const Graph & rhs)
//    {
//          std::cout << "Graph copy constructor" << std::endl;
//          _ccdbg = rhs._ccdbg;
//    };
//    // copy operator
//    Graph & operator = (const Graph & rhs)
//    {
//          std::cout << "Graph copy operator" << std::endl;
//          if(this != &rhs)
//          {
//              _ccdbg = rhs._ccdbg;
//          }
//          return *this;
//    };
//
//    // move
//    Graph(const Graph && rhs) noexcept
//    {
//          std::cout << "Graph move constructor" << std::endl;
//          _ccdbg = std::move(rhs._ccdbg);
//    };
//    // move operator
//    Graph & operator = (const Graph && rhs) noexcept
//    {
//          std::cout << "Graph move operator" << std::endl;
//          if(this != &rhs)
//          {
//              _ccdbg = std::move(rhs._ccdbg);
//          }
//          return *this;
//    };
//
//    // destructor
//    ~Graph()
//    {
//        std::cout << "Graph destructor" << std::endl;
//    };

    // build new bifrost graph and index
    void build(const std::string& infile1,
               const size_t kmer,
               const size_t gap,
               size_t num_threads,
               bool is_ref,
               const std::string& infile2,
               const std::string& outpref);

    // read existing graph and index
    void read(const std::string& sk_index);

    double query(const std::string& query);

    private:
    // stored bifrost DBG
    CompactedDBG<> _cdbg;
    size_t _kmer;
    size_t _gap;
    double _kmerd;

    KmerMap _kmermap;

    void index_split_kmer();
    void read_skindex(const std::string& sk_index);
};

CompactedDBG<> buildGraph (const std::string& infile_1,
                          const std::string& infile_2,
                          const bool is_ref,
                          const int kmer,
                          const int threads,
                          const std::string& output_prefix);

#endif //BIFROST_H