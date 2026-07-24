// #include <bits/stdc++.h>
#include <iostream>
#include <sstream>
#include <ctime>
#include <set>
#include <tuple>
#include "hypergraph.h"
#include "readhg.h"
#include "algorithms.h"

// typedef std::map<std::string, std::string> strstrMap;
typedef std::map<std::string, std::string> strstrMap;
typedef std::tuple<size_t,size_t,size_t> inttriplet;
typedef std::vector<inttriplet> vinttriplet;

template <typename T>
bool isSubsetOrEqual(std::set<T> const& a, std::set<T> const& b) {
   for(auto const& av:a){
      if(std::find(b.begin(),b.end(),av)==b.end())
          return false;
   }
   return true;
}
int check_conditiondeg(Hypergraph& h, intIntMap& core){
    /*
    Checks that the sub-hypergraph induced by all nodes v with c(v)>=k has at least k incident hyperedges
    in that sub-hypergraph \forall k \in [min_v c(v), max_v c(v) ]. (Coreness condition)
    */
    std::set<size_t> core_values;
    int incorrect = 0;
    for (auto node : h.init_nodes){
        core_values.insert(core[node]);
    }
    bool condition_true = true;
    for(auto p: core_values){
        std::set<std::string> subnodes;
        for(auto node: h.init_nodes){
            auto node_str = std::to_string(node);
            if (core[node]>= p){
                subnodes.insert(node_str);
            }
        }
        Hypergraph subH;
        size_t count = 0;
        for(auto y: h.hyperedges){
            std::vector<std::string> strvecE(y.size());
            for(auto u:y){
                strvecE.push_back(std::to_string(u));
            }
            if (isSubsetOrEqual(std::set<std::string>(strvecE.begin(),strvecE.end()),subnodes)){
                subH.addEdge(count++,strvecE);
            }
        }
        intIntMap node_deg; 
        for(auto y: subH.hyperedges){
            for(auto u:y){
                if (node_deg.find(u) == node_deg.end())
                    node_deg[u] = 1;
                else 
                    node_deg[u] += 1;
            }
        }
        condition_true = true;
        std::string violoating_node;
        for(auto node: subH.init_nodes){
            if (node_deg[node]< p){
                condition_true = false;
                violoating_node = node;
            }
        }
        if (!condition_true){
            std::cout<< "violating core: ("<<p<<")\n";
            std::cout << "violating node: "<<violoating_node<<"\n";
            incorrect+=1;
        }
    }
    if (condition_true){
        std::cout<< "No violation\n";
    }
    return incorrect;
}
int check_conditionnbr(Hypergraph& h, intIntMap& core){
    /*
    Checks that the sub-hypergraph induced by all nodes v with c(v)>=k has at least k neighhbors
    in that sub-hypergraph \forall k \in [min_v c(v), max_v c(v) ]. (Coreness condition)
    */
    std::set<size_t> core_values;
    int incorrect = 0;
    for (auto node : h.init_nodes){
        core_values.insert(core[node]);
    }
    bool condition_true = true;
    for(auto p: core_values){
        // auto p = x.second;
        // std::cout<<p<<" - "<<s<<"\n";
        std::set<std::string> subnodes;
        for(auto node: h.init_nodes){
            // if (core[node]>= p){
            if (core[node]>= p){
                auto node_str = std::to_string(node);
                subnodes.insert(node_str);
            }
        }
        // for(auto u: subnodes)   std::cout<<u<<" ";
        // std::cout<<"\n";
        Hypergraph subH;
        size_t count = 0;
        for(auto y: h.hyperedges){
            std::vector<std::string> strvecE(y.size());
            for(auto u:y){
                strvecE.push_back(std::to_string(u));
            }
            if (isSubsetOrEqual(std::set<std::string>(strvecE.begin(),strvecE.end()),subnodes)){
                subH.addEdge(count++,strvecE);
            }
        }
        // subH.initialise();
        // count neighbors
        intuSetintMap init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
        // intintvec edges( e_id_to_edge.size() ,intvec{});
        for(auto elem:subH.hyperedges){
            // auto elem = e_id_to_edge[eid];
            auto edge_sz = elem.size();
            for(auto v_id: elem){
                if(init_nbr.find(v_id) == init_nbr.end()){
                    init_nbr[v_id] = std::unordered_set<size_t>();
                }
                else{
                    auto _tmp = &init_nbr[v_id];
                    for (auto u: elem){
                        if (u!=v_id){
                            _tmp->insert(u);
                        }
                    }
                }
            }
        }
        condition_true = true;
        std::string violoating_node;
        for(auto node: subH.init_nodes){
            // if (subH.init_nbrsize[node]< p){
            if (init_nbr[node].size() < p){
                condition_true = false;
                violoating_node = node;
            }
        }
        if (!condition_true){
            std::cout<< "violating core: ("<<p<<")\n";
            std::cout << "violating node: "<<violoating_node<<"\n";
            incorrect+=1;
            // std::cout<<"sub-hyp.\n";
            // subH.printHypergraph();
            // std::cout<<"violating node nbr: ";
            // for (auto u: subH.init_nbr[violoating_node]) std::cout<< u<<" ";
            // std::cout<<"\n";
        }
    }
    if (condition_true){
        std::cout<< "No violation\n";
    }
    return incorrect;
}

int check_conditionkd(Hypergraph& h, vinttriplet& kdcores){
    /*
    Checks that the sub-hypergraph induced by all nodes v with c(v)>=k has at least k neighhbors
    in that sub-hypergraph \forall k \in [min_v c(v), max_v c(v) ]. (Coreness condition)
    */
    std::set<std::pair<size_t,size_t>> core_pairs;
    int incorrect = 0;
    for (auto [node,f,s] : kdcores){
        std::pair<size_t,size_t> x = std::make_pair(f,s);
        core_pairs.insert(x);
    }
    bool condition_true = true;
    for(auto x: core_pairs){
        auto p = x.first;
        auto s = x.second;
        std::set<std::string> subnodes;
        for(auto [node,f,sec]: kdcores){
            // if (core[node]>= p){
            if (f>= p && sec>= s){
                auto node_str = std::to_string(node);
                subnodes.insert(node_str);
            }
        }
        Hypergraph subH;
        size_t count = 0;
        for(auto y: h.hyperedges){
            std::vector<std::string> strvecE(y.size());
            for(auto u:y){
                strvecE.push_back(std::to_string(u));
            }
            if (isSubsetOrEqual(std::set<std::string>(strvecE.begin(),strvecE.end()),subnodes)){
                subH.addEdge(count++,strvecE);
            }
        }
        intuSetintMap init_nbr;  //# key => node id, value => List of Neighbours. (use hashtable instead of dictionary => Faster on large |V| datasets. )
        std::map<size_t,size_t> node_deg;
        for(auto elem:subH.hyperedges){
            // auto elem = e_id_to_edge[eid];
            auto edge_sz = elem.size();
            for(auto v_id: elem){
                if(init_nbr.find(v_id) == init_nbr.end()){
                    init_nbr[v_id] = std::unordered_set<size_t>();
                }
                else{
                    auto _tmp = &init_nbr[v_id];
                    for (auto u: elem){
                        if (u!=v_id){
                            _tmp->insert(u);
                        }
                    }
                }
                if (node_deg.find(v_id) == node_deg.end())
                    node_deg[v_id] = 1;
                else 
                    node_deg[v_id] += 1;
            }
        }
        size_t violoating_node;
        for(auto node: subH.init_nodes){
            // if (subH.init_nbrsize[node]< p){
            if (init_nbr[node].size() < p || node_deg[node]< s){
                condition_true = false;
                violoating_node = node;
            }
        }
        if (!condition_true){
            std::cout<< "violating core: ("<<p<<","<<s<<")\n";
            std::cout << "violating node: "<<violoating_node<<"\n";
            incorrect+=1;
        }
    }

    if (condition_true){
        std::cout<< "No violation\n";
    }
    return incorrect;
}
int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <dataset> <algorithm>\n";
        return 1;
    }

    std::string dataset = argv[1];
    std::string algorithm = argv[2];

    Hypergraph h;
    getHg(dataset, h);
    h.dataset = dataset;
    h.initialise();
    Algorithm a(h);
    a.output["algo"] = algorithm;

    if (algorithm == "nbr_k_coreness"){  
        local_core_OPTIV(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, false);
        a.writecore();
        a.write_time();
    }
    if (algorithm=="kd_coreness"){
        kdCorehybrid(h.dataset, h.hyperedges, h.init_nodes, h.node_index, a, false);
        a.writekdcore_distinct();
        a.write_time();
    }
}
