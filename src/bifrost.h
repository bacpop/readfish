#include <bifrost/CompactedDBG.hpp>

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
               const int kmer,
               size_t num_threads,
               bool is_ref,
               const std::string& infile2,
               const std::string& outpref);

    // read existing graph and index
    void read(const std::string& graphfile);

    double query(const std::string& query);

    // clear graph object
    void clear() {_cdbg.clear();};

    private:
    // stored bifrost DBG
    CompactedDBG<> _cdbg;
    int _kmer;
};

CompactedDBG<> buildGraph (const std::string& infile_1,
                          const std::string& infile_2,
                          const bool is_ref,
                          const int kmer,
                          const int threads,
                          const std::string& output_prefix);

