#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include "hypergraph.h"
typedef  std::unordered_map<size_t, size_t> intIntMap;
typedef  std::map<std::string, size_t> strIntMap;
typedef  std::map<std::string, std::vector<size_t>> strvIntMap;
typedef  std::map<std::string, std::set<size_t>> strsIntMap;
typedef  std::map<size_t, std::set<std::string>> intsStrMap;
typedef  std::map<std::string, std::vector<std::string>> strvStrMap;
typedef  std::map<size_t, std::vector<std::string>> intvStrMap;
typedef  std::vector<std::string> strvec;
typedef  std::set<std::string> strset;
typedef  std::vector<size_t> intvec;
typedef std::vector<std::pair<std::string, std::string>> strstrprvec;
typedef std::map<std::string, std::string> strstrMap;
typedef std::unordered_map <size_t,bool> intboolMap;
typedef std::unordered_set<size_t> uintSet;
typedef std::vector<uintSet > uintsetvec;
typedef std::unordered_map<size_t, uintSet> intuSetintMap;
typedef std::map<std::string, std::string> strstrMap;
typedef std::vector< intvec > intintvec;
typedef std::pair<size_t,size_t> intpair;
typedef std::tuple<size_t,size_t,size_t> inttriplet;
typedef std::vector<inttriplet> vinttriplet;
class Algorithm{
    Hypergraph hg;
    public:
    intIntMap core;
    vinttriplet kdcores;
    double exec_time = 0;
    double core_exec_time = 0;
    double correction_time = 0;
    size_t nu_cu = 0;
    size_t num_nbr_queries = 0;
    strstrMap output;
    std::vector< strstrMap > hnlog;
    strstrMap timelogs;
    Algorithm( Hypergraph &H);
    ~Algorithm();
    void printcore();
    void writecore(std::string folder="../../output/");
    void writekdcore_distinct(std::string folder="../../output/", bool verbose=false);
    void write_time(std::string folder="../../output/");
};
void print_bucket( intuSetintMap&, intvec&);
void local_core_OPTIV( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log);
void kdCorehybrid( std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log);
#endif