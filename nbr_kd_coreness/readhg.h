#ifndef READHG_H
#define READHG_H

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include "hypergraph.h"

template <typename Out>
inline void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

inline void getHg(const std::string& dataname, Hypergraph &hg){
    const std::string filepath = "../../datasets/" + dataname + "/network.hyp";
    // std::cout << "[DEBUG] filepath = " << filepath << std::endl;
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        std::cerr << "[getHg] Failed to open: " << filepath << "\n";
        return;
    }

    std::string line;
    size_t i = 0;

    while (std::getline(infile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.find_first_not_of(" \t") == std::string::npos) continue;

        std::vector<std::string> tokens;
        if (line.find(',') != std::string::npos) {
            tokens = split(line, ',');
        } else if (line.find('\t') != std::string::npos) {
            tokens = split(line, '\t');
        } else {
            std::istringstream iss(line);
            std::string token;
            while (iss >> token) tokens.push_back(token);
        }

        tokens.erase(std::remove_if(tokens.begin(), tokens.end(),
                      [](const std::string& s){ return s.empty(); }),
                      tokens.end());

        if (!tokens.empty()) {
            hg.addEdge(i, tokens);
            ++i;
        }
    }
}

#endif // READHG_H