#ifndef BIFROST_H
#define BIFROST_H

#include <vector>
#include <string>
#include <bitset>
#include <robin_hood.h>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <progressbar.h>
#include <bifrost/CompactedDBG.hpp>
#include <math.h>
#include <stack>

// openmp
#include <omp.h>

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/unordered_set.hpp>
#include <cereal/types/unordered_map.hpp>

// global variable declaration
namespace py = pybind11;

typedef std::vector<Kmer> KmerVec;
typedef std::unordered_map<uint64_t, std::unordered_set<uint64_t>> StdKmerMap;
typedef robin_hood::unordered_map<uint64_t, robin_hood::unordered_set<uint64_t>> KmerMap;

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
    void build(const std::string &infile1,
               const size_t kmer,
               const size_t gap,
               size_t num_threads,
               bool is_ref,
               const std::string &infile2,
               const std::string &outpref,
               const bool splitk);

    // read existing graph and index
    void read(const std::string &infile);

    std::pair<double, double> query(const std::string &query);

private:
    // stored bifrost DBG
    CompactedDBG<> _cdbg;
    size_t _kmer;
    size_t _gap;
    double _gapd;
    double _kmerd;
    bool _splitk;
    double _param1;

    KmerMap _kmermap;

    void index_split_kmer();

    void read_skindex(const std::string &sk_index);

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(_gap, _kmer, _kmermap);
    };
};

CompactedDBG<> buildGraph (const std::string& infile_1,
                          const std::string& infile_2,
                          const bool is_ref,
                          const int kmer,
                          const int threads,
                          const std::string& output_prefix);

// SKA2 functions

//// 64-bit random seeds corresponding to bases and their complements
static const uint64_t seedA = 0b00;
static const uint64_t seedC = 0b01;
static const uint64_t seedG = 0b10;
static const uint64_t seedT = 0b11;
//static const uint64_t seedN = 0b100; // this won't be used

static const uint64_t seedM = 0b0100;
static const uint64_t seedR = 0b0101;
static const uint64_t seedW = 0b0110;
static const uint64_t seedS = 0b0111;
static const uint64_t seedY = 0b1000;
static const uint64_t seedK = 0b1001;
static const uint64_t seedV = 0b1010;
static const uint64_t seedH = 0b1011;
static const uint64_t seedD = 0b1100;
static const uint64_t seedB = 0b1101;
static const uint64_t seedN = 0b1110;
static const uint64_t seedX = 0b1111;

static const uint64_t look_up_table[256] = {
        seedX, seedT, seedX, seedG, seedA, seedA, seedX, seedC, // 0..7
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 8..15
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 16..23
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 24..31
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 32..39
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 40..47
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 48..55
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 56..63
        seedX, seedA, seedB, seedC, seedD, seedX, seedX, seedG, // 64..71
        seedH, seedX, seedX, seedK, seedX, seedM, seedN, seedX, // 72..79
        seedX, seedX, seedR, seedS, seedT, seedT, seedV, seedW, // 80..87
        seedX, seedY, seedX, seedX, seedX, seedX, seedX, seedX, // 88..95
        seedX, seedA, seedB, seedC, seedD, seedX, seedX, seedG, // 96..103
        seedH, seedX, seedX, seedK, seedX, seedM, seedN, seedX, // 104..111
        seedX, seedX, seedR, seedS, seedT, seedT, seedV, seedW, // 112..119
        seedX, seedY, seedX, seedX, seedX, seedX, seedX, seedX, // 120..127
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 128..135
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 136..143
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 144..151
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 152..159
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 160..167
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 168..175
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 176..183
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 184..191
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 192..199
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 200..207
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 208..215
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 216..223
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 224..231
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 232..239
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX, // 240..247
        seedX, seedX, seedX, seedX, seedX, seedX, seedX, seedX  // 248..255
};

uint64_t to_binary(std::string& current_kmer);

uint64_t ReverseComp64(const uint64_t mer, uint8_t kmerSize);
#endif //BIFROST_H