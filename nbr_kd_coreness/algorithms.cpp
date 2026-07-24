#include <filesystem>
#include <vector>
#include <set>
#include <map>
#include <string> 
#include <iostream>
#include <set>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include "hypergraph.h"
#include "algorithms.h"
#include "utils.h"

void Algorithm::printcore(){
    //std::cout << "core: \n";
    for(const auto& elem : core)
    {
    std::cout << elem.first << "->"<<elem.second<<"\n";
    }
}

using KDTriplet = std::tuple<size_t,size_t,size_t>;         // (v,k,d)
using KDList    = std::vector<std::pair<size_t,size_t>>;    // [(k,d)]
using KDMap     = std::unordered_map<size_t, KDList>;       // v -> [(k,d)]

std::string getDatasetName(const std::string& path) {
    std::filesystem::path p(path);
    return p.filename().string();
}

static KDMap build_per_node_kd_unique(const std::vector<KDTriplet>& kdcores){
    struct H{
        size_t operator()(const KDTriplet&t)const{
            auto [v,k,d]=t;
            // 간단 해시
            return std::hash<unsigned long long>{}(
                ((unsigned long long)v<<42) ^ ((unsigned long long)k<<21) ^ (unsigned long long)d
            );
        }
    };
    std::unordered_set<KDTriplet,H> uniq; uniq.reserve(kdcores.size()*2);
    KDMap per; per.reserve(kdcores.size()/2+1);
    for(const auto& t: kdcores){
        if(uniq.insert(t).second){
            size_t v,k,d; std::tie(v,k,d)=t;
            per[v].emplace_back(k,d);
        }
    }
    return per;
}

static KDList pareto_skyline_one(const KDList& in, bool verbose=false, size_t node_id=(size_t)-1){
    if(in.empty()) return {};
    KDList v=in;

    // k desc, d desc
    std::sort(v.begin(), v.end(),
        [](const auto& A, const auto& B){
            return (A.first!=B.first)? A.first>B.first : A.second>B.second;
        });

    KDList keep; keep.reserve(v.size());
    long long best_d=-1; size_t best_k=0;

    for(const auto& kd: v){
        size_t k=kd.first, d=kd.second;
        if((long long)d > best_d){
            keep.emplace_back(k,d);
            best_d=(long long)d; best_k=k;
        }else if(verbose){
            std::cout<<"[DROP] node "<<node_id<<" ("<<k<<","<<d<<") <= dominated by ("
                     <<best_k<<","<<best_d<<")\n";
        }
    }

    // 보기 좋게 k asc, d asc
    std::sort(keep.begin(), keep.end(),
        [](const auto& A, const auto& B){
            return (A.first!=B.first)? A.first<B.first : A.second<B.second;
        });
    return keep;
}

static KDMap pareto_filter_kd_map(const KDMap& per, bool verbose=false){
    KDMap out; out.reserve(per.size());
    for(const auto& kv: per){
        out[kv.first]=pareto_skyline_one(kv.second, verbose, kv.first);
    }
    return out;
}

static void write_kd_map_csv(const std::string& folder, const KDMap& mp, const std::string& fname){
    namespace fs=std::filesystem;
    fs::create_directories(folder);
    std::ofstream f(fs::path(folder)/fname, std::ios::trunc);
    if(!f) return;

    f<<"node_id,coreness\n";
    for(const auto& kv: mp){
        f<<kv.first<<",\"{";
        const auto& pairs=kv.second;
        for(size_t i=0;i<pairs.size();++i){
            auto [k,d]=pairs[i];
            f<<"("<<k<<", "<<d<<")"<<(i+1<pairs.size()? ", ":"");
        }
        f<<"}\"\n";
    }
}

void Algorithm::writekdcore_distinct(std::string folder, bool verbose){
    auto per_node = build_per_node_kd_unique(this->kdcores);
    auto skyline  = pareto_filter_kd_map(per_node, verbose);
    const std::string fname = getDatasetName(hg.dataset) + "_" + output["algo"] + ".csv";
    write_kd_map_csv(folder, skyline, fname);
    std::cout<<"File saved to: "<<folder<<fname<<"\n";
}

void Algorithm::writecore(std::string folder){

    std::string file = folder + getDatasetName(hg.dataset)+"_"+output["algo"]+".csv";
    std::cout<<"writing to: "<<file<<"\n";

    std::stringstream ss;
    ss << "node_id,coreness\n";
    for (auto elem : core) {
        ss << std::to_string(elem.first) << "," << std::to_string(elem.second) << "\n";
    }
    std::ofstream out(file.c_str());
    if(out.fail())
    {
        out.close();
    }
    out << ss.str();
    out.close();
    std::cout<<"File saved to: "<<file<<"\n";
}

void Algorithm::write_time(std::string folder) {
    namespace fs = std::filesystem;
    fs::create_directories(folder);

    const std::string fname = "time.csv";
    const fs::path fpath = fs::path(folder) / fname;

    const std::string dataset = getDatasetName(hg.dataset);
    const std::string algo = output["algo"];
    const std::string time = output["execution time"];

    bool file_exists = fs::exists(fpath);

    std::ofstream out(fpath, std::ios::app);  // append 모드
    if (!out) return;

    if (!file_exists) {
        out << "algorithm,dataset,time\n";  // 헤더 추가
    }

    out << algo << "," << dataset << "," << time << "\n";
    out.close();
    std::cout<<"File saved to: "<<fpath.string()<<"\n";
}

Algorithm::Algorithm(Hypergraph &H){
    hg = H;
    output["dataset"] = H.dataset;
}
Algorithm::~Algorithm(){}



void iterate_nbrs(size_t v, intvec & nbrs, uintsetvec & inc_dict, intintvec& e_id_to_edge){
    /* Returns the set of neighbours of v.
        implements a traversal from vertex v to each of its neighbours in contrast to set in neighbors(). 
        It also returns an iterator. So it avoids creating the neighborhood list explicitely.
        Overall complexity: O(d(v) * |e_max|), where e_max = largest hyperedge 
    */
    auto incident_edges = inc_dict[v];  //# {O(1)}
    if (incident_edges.size()){
        intboolMap visited_dict;
        for (auto e_id : incident_edges){  //# { O(d(v)) }
            for (auto u : e_id_to_edge[e_id]){  //# { O(|e|)}
                if (u != v){
                    if (visited_dict.find(u) == visited_dict.end()){ // u does not exists
                        visited_dict[u] = true;
                        nbrs.push_back(u);
                    }
                    if (!visited_dict[u]){
                        visited_dict[u] = true;
                        nbrs.push_back(u);
                    }
                }
            }
        }
    }
}

void removeV_transform(size_t v, uintsetvec & inc_dict, intintvec &e_id_to_edge){
    // For every edge e_id incident on v, remove e_id from every vertex u in edge[e_id] distinct from v
    for (auto e_id : inc_dict[v]){
        for(auto u: e_id_to_edge[e_id]){
            if (u==v)   continue;
            inc_dict[u].erase(e_id);
        }
    }
    inc_dict[v] = uintSet();
}

size_t get_number_of_nbrs(size_t v, uintsetvec & inc_dict, intintvec &e_id_to_edge){
    std::unordered_set<int> nbrs;
    auto incident_edges = inc_dict[v];  //# {O(1)}
    if (incident_edges.size()){
        intboolMap visited_dict;
        for (auto e_id : incident_edges){  //# { O(d(v)) }
            for (auto u : e_id_to_edge[e_id]){  //# { O(|e|)}
                if (u != v){
                    nbrs.insert(u);
                }
            }
        }
    }
    return nbrs.size();
}

bool LCCSAT_check_OPTIV(size_t inc_edges_F[], size_t inc_edges_N[],intvec& min_hindices, std::vector<intvec>&edges, size_t u_id, size_t core_u){

    std::set <size_t> Nplus;
    for (size_t i = inc_edges_N[u_id]; i<inc_edges_N[u_id+1]; i++){
        size_t e_id = inc_edges_F[i];
        if (min_hindices[e_id] >= core_u){
            for(auto v_id: edges[e_id]){
                Nplus.insert(v_id);
            }
            if(Nplus.size()-1>=core_u)
            return true;
        }
    }
    auto sz_Nplus = Nplus.size();
    if (sz_Nplus>0){
        if (sz_Nplus-1 >= core_u)
            return true;
        else
            return false;
    }
    else
        return false;
}

bool LCCSAT_OPTIV(std::vector<intvec> & min_hindex_to_edge, std::vector<intvec>&edges, size_t core_u, std::set<size_t> & Nplus){

    for(auto eid : min_hindex_to_edge[core_u]){
        for(auto v : edges[eid]){
            Nplus.insert(v);
        }
    }
    auto sz_Nplus = Nplus.size();
    if (sz_Nplus>0){
        if (sz_Nplus-1 >= core_u)
            return true;
        else
            return false;
    }
    else
        return false;
}

size_t core_correct_OPTIV(size_t inc_edges_F[], size_t inc_edges_N[],intvec& min_hindices, std::vector<intvec>&edges, size_t u_id, size_t core_u){   // This function is used in optimised local_correct_optIII
        // """ Finds the correct \hat{h} by traversing in descending order from core_u, core_u-1,...,until correct."""

        // time_t start,end,start1, end1;
        core_u = core_u - 1;

        // start1 = clock();
        std::vector<intvec> min_hindex_to_edge(core_u+1);
        for (size_t i = inc_edges_N[u_id]; i<inc_edges_N[u_id+1]; i++){
            size_t eid = inc_edges_F[i];
            if(min_hindices[eid]>=core_u){
                min_hindex_to_edge[core_u].push_back(eid);
            }else 
                {
                min_hindex_to_edge[min_hindices[eid]].push_back(eid);
            }
        }
        // end1 = clock();
        std::set<size_t> Nplus;
        
        while (LCCSAT_OPTIV(min_hindex_to_edge, edges, core_u, Nplus) == false) {
            core_u = core_u - 1;
        }

        return core_u;
}

void local_core_OPTIV( std::string dataset, intintvec &e_id_to_edge, intvec& init_nodes, intIntMap& node_index, Algorithm& a, bool log){
    clock_t start, end;
    clock_t start1,end1,start2,end2,start3,end3;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    size_t N = init_nodes.size();
    size_t M = e_id_to_edge.size();
    intvec pcore(N); //
    // strIntMap node_index; //key = node id (string), value = array index of node (integer)
    intvec llb(N,0); // key => node id (v), value => max(|em|-1) for all edge em incident on v 
    size_t glb = std::numeric_limits<size_t>::max();
    intintvec edges( M ,intvec{}); // i = edge_id, value = vector of vertices in e[edge_id]
    intvec min_e_hindex(M);
    intintvec inc_edges(N, intvec{}); // i=node_id, value = vector of edge ids incident on node_id
    uintsetvec nbrs(N, std::unordered_set<size_t>{});

    // compute initial neighbors and number of neighbors
    start3 = clock();
    for(size_t eid= 0; eid<M; eid++){
        auto elem = e_id_to_edge[eid];
        sz_inc_edge += elem.size();
        for(auto v_id: elem){
            auto j = node_index[v_id];
            llb[j] = std::max(elem.size()-1,llb[j]);
            inc_edges[j].push_back(eid);
            edges[eid].push_back(j);
            auto _tmp = &nbrs[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
        // std::cout<<"\n";
    }
    for(size_t _i = 0; _i< N; _i ++){
        sz_init_nbrs += nbrs[_i].size();
    }

    // std::cout<<"Init nbrs "<<sz_init_nbrs<<" Inc_edges "<<sz_inc_edge<<"\n"; 
    end3 = clock();
    // std::cout<<"Time for init_nbr calculation "<<double(end3 - start3) / double(CLOCKS_PER_SEC)<<"\n";
    time_t start4,end4;
    start4 = clock();

    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    size_t *inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));
    inc_edges_N[0] = 0;
    nbrs_N[0] = 0;
    for (int _i = 1; _i<= N; _i ++){
        auto nbr_i = nbrs[_i-1].size();
        nbrs_N[_i] = nbrs_N[_i-1] + nbr_i;
        glb = std::min(glb, nbr_i);
		inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[_i-1].size();
    }
    
    // Calculate csr representation for incident edges
    for (int _i = 1; _i<= N; _i ++){
		auto _index = nbrs_N[_i-1];
		for(auto u: nbrs[_i-1]){
			nbrs_F[_index++] = u;
		}
		_index = inc_edges_N[_i-1];
		for(auto eid : inc_edges[_i-1])
			inc_edges_F[_index++] = eid;
	}
    a.output["init_time"] = std::to_string(double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();

    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrs[i].size(); // initialize pcore
        llb[i] = std::max(llb[i],glb);
    }
    
    if (log){
        strstrMap h0;
        for(int i=0;i<N;i++){
            h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
        }
        h0["Time"] = "0";
        a.hnlog.push_back(h0);
    }

    // std::vector<size_t> hn(N);
    size_t iterations = 0;
    size_t correction_number=0, check=0;
    time_t start_main, end_main,start_h,end_h, start_minh, end_minh;
    double hindext = 0, minht = 0;
    start_main = clock();
    while (1){
        iterations+=1;
        bool flag = true;
        // compute h-index and update core
        start_h = clock();
        for(size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            size_t H_value = hIndex_csr(nbrs_N[i],nbrs_N[i+1],nbrs_F,pcore);
            if (H_value < pcore[i]) 
                pcore[i] = H_value;     //pcore[i] is same as hvn here
        }
        end_h = clock();
        hindext += double(end_h - start_h) / double(CLOCKS_PER_SEC);
        // for every edge update is minimum of hindex(constituent vertices)
        start_minh = clock(); 
        for(size_t i = 0; i< M; i++){
            size_t _min = N+1;      //Why M+1 and why calculate 
            for (auto u_id: edges[i])   _min = std::min(_min,pcore[u_id]);
            min_e_hindex[i] = _min;
        }
        end_minh = clock();
        minht += double(end_minh - start_minh) / double(CLOCKS_PER_SEC);
        start1 = clock();
        for (size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            bool lccsat = LCCSAT_check_OPTIV(inc_edges_F,inc_edges_N,min_e_hindex,edges,i,pcore[i]);
            if (lccsat == false){ 
                start2 = clock();  
                correction_number++;
                flag = flag && false;   
                auto hhatn = core_correct_OPTIV(inc_edges_F,inc_edges_N,min_e_hindex,edges,i,pcore[i]);
                pcore[i]  = hhatn;
                for (size_t j = inc_edges_N[i]; j<inc_edges_N[i+1]; j++){
                    size_t e_id = inc_edges_F[j];
                    if (min_e_hindex[e_id] >= pcore[i]){
                        min_e_hindex[e_id] = pcore[i];
                    }
                }
                end2 = clock();
                a.correction_time+= double(end2-start2)/ double(CLOCKS_PER_SEC);
            }
        }
        end1 = clock();
        a.core_exec_time += double(end1 - start1) / double(CLOCKS_PER_SEC);
        end1 = clock();
        if (log){
            strstrMap h0;
            for(size_t i=0;i<N;i++){
                h0[std::to_string(init_nodes[i])] = std::to_string(pcore[i]);
            }
            h0["Time"] = std::to_string(double(end1 - start_main) / double(CLOCKS_PER_SEC));
            a.hnlog.push_back(h0);
        }
        if (flag)
            break;

    }
    end_main = clock();
    end = clock();
    for(size_t i=0; i<N; i++)
    {
        auto node = init_nodes[i];
        a.core[node] = pcore[i];
        a.nu_cu += nbrs[i].size() - pcore[i];
    }
    a.exec_time = double(end - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);

}

void print_bucket(intuSetintMap& degbucket, intvec& init_nodes){
    for(auto pr: degbucket) {
        std::cout<<pr.first<<": ";
        for(auto u: pr.second) std::cout<< init_nodes[u]<<", "; 
        std::cout<<"\n";
    }
}
void kdCorehybrid(std::string dataset, intintvec e_id_to_edge, intvec init_nodes, intIntMap& node_index, Algorithm& a, bool log){
    clock_t start, end;
    double init_tm=0, actual_tm=0;
    /* Recording the starting clock tick.*/
    start = clock();
    size_t sz_init_nbrs = 0;    // stores the number of initial neighbours for all vertices
    size_t sz_inc_edge = 0;     // stores the number of incident edges for all vertices
    size_t N = init_nodes.size();
    size_t M = e_id_to_edge.size();
    intvec pcore(N); //
    intvec llb(N,0); // key => node id (v), value => max(|em|-1) for all edge em incident on v 
    size_t glb = std::numeric_limits<size_t>::max();
    intintvec edges( M ,intvec{}); // i = edge_id, value = vector of vertices in e[edge_id]
    intvec min_e_hindex(M);
    uintsetvec inc_edges(N, uintSet{});
    uintsetvec nbrs(N, std::unordered_set<size_t>{});
    // std::vector<std::vector<intpair>> score_mult(N,std::vector<intpair>{});
    // compute initial neighbors and number of neighbors
    for(size_t eid= 0; eid<M; eid++){
        auto elem = e_id_to_edge[eid];
        sz_inc_edge += elem.size();
        for(auto v_id: elem){
            auto j = node_index[v_id];
            inc_edges[j].insert(eid);
            edges[eid].push_back(j);
            auto _tmp = &nbrs[j];
            for (auto u: elem){
                if (u!=v_id){
                    _tmp->insert(node_index[u]);
                }
            }
        }
    }
    for(size_t _i = 0; _i< N; _i ++){
        sz_init_nbrs += nbrs[_i].size();
    }

    size_t* inc_edges_F = (size_t*)malloc(sz_inc_edge*sizeof(size_t));
    size_t *inc_edges_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_N = (size_t*)malloc((N+1)*sizeof(size_t));
    size_t* nbrs_F = (size_t*)malloc(sz_init_nbrs*sizeof(size_t));
    inc_edges_N[0] = 0;
    nbrs_N[0] = 0;
    for (int _i = 1; _i<= N; _i ++){
        auto nbr_i = nbrs[_i-1].size();
        nbrs_N[_i] = nbrs_N[_i-1] + nbr_i;
        glb = std::min(glb, nbr_i);
		inc_edges_N[_i] = inc_edges_N[_i-1] + inc_edges[_i-1].size();
    }
    
    // Calculate csr representation for incident edges
    for (int _i = 1; _i<= N; _i ++){
		auto _index = nbrs_N[_i-1];
		for(auto u: nbrs[_i-1]){
			nbrs_F[_index++] = u;
		}
		_index = inc_edges_N[_i-1];
		for(auto eid : inc_edges[_i-1])
			inc_edges_F[_index++] = eid;
	}
    init_tm = (double(clock() - start) / double(CLOCKS_PER_SEC));
    start = clock();

    // initialise core to a upper bound
    for (size_t i = 0; i < N; i++){
        pcore[i] = nbrs[i].size(); // initialize pcore
        llb[i] = std::max(llb[i],glb);
    }

    // std::vector<size_t> hn(N);
    size_t iterations = 0;
    while (1){
        iterations+=1;
        bool flag = true;
        // compute h-index and update core
        for(size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            size_t H_value = hIndex_csr(nbrs_N[i],nbrs_N[i+1],nbrs_F,pcore);
            if (H_value < pcore[i]) 
                pcore[i] = H_value;     //pcore[i] is same as hvn here
        }
        for(size_t i = 0; i< M; i++){
            size_t _min = N+1;      //Why M+1 and why calculate 
            for (auto u_id: edges[i])   _min = std::min(_min,pcore[u_id]);
            min_e_hindex[i] = _min;
        }
        for (size_t i = 0; i<N; i++){
            if (pcore[i] == llb[i]) continue;
            bool lccsat = LCCSAT_check_OPTIV(inc_edges_F,inc_edges_N,min_e_hindex,edges,i,pcore[i]);
            if (lccsat == false){ 
                flag = flag && false;   
                auto hhatn = core_correct_OPTIV(inc_edges_F,inc_edges_N,min_e_hindex,edges,i,pcore[i]);
                pcore[i]  = hhatn;
                for (size_t j = inc_edges_N[i]; j<inc_edges_N[i+1]; j++){
                    size_t e_id = inc_edges_F[j];
                    if (min_e_hindex[e_id] >= pcore[i]){
                        min_e_hindex[e_id] = pcore[i];
                    }
                }
            }
        }
        if (flag)
            break;

    }
    // Peeling iteration to find secondary core.
    // intuSetintMap nbrbucket;
	// initialize every nodes initial bucket to the primary core-number.
	size_t min_cv = N+1;
    size_t max_cv = 0;
    for (size_t i = 0; i<N; i++){
        auto cv = pcore[i];
        // auto v = i;
		// if (nbrbucket.find(cv) == nbrbucket.end())
		// 	nbrbucket[cv] = uintSet({v});
		// else
		// 	nbrbucket[cv].insert(v);
        max_cv = std::max(max_cv, cv);
        min_cv = std::min(min_cv, cv);
    }
    // std::cout<<min_cv<<":"<<max_cv<<"\n";
	for (size_t pk = 1; pk<= max_cv; pk++){
        // intvec score(N,-1); //
        if(log) std::cout<<"pcore="<<pk<<"\n";
		// deg bucket init 
        intvec inverse_bucket(N);
		intuSetintMap degbucket;
		size_t max_deg = 0;
        // std::vector<bool>stop(N,false);
        inc_edges = uintsetvec(N, uintSet{});
        edges = intintvec(M, intvec{});
        nbrs = uintsetvec(N, std::unordered_set<size_t>{});
    //         intintvec edges( M ,intvec{}); // i = edge_id, value = vector of vertices in e[edge_id]
    // intvec min_e_hindex(M);
    // uintsetvec inc_edges(N, uintSet{});
        for(size_t eid= 0; eid<M; eid++){
            auto elem = e_id_to_edge[eid];
            bool flag = false;
            for(auto v_id: elem){
                auto j = node_index[v_id];
                if (pcore[j]<pk){
                    flag = true;
                    break;
                }
            }
            if (!flag){
                for(auto v_id: elem){
                    auto j = node_index[v_id];
                    inc_edges[j].insert(eid);
                    edges[eid].push_back(j);
                    auto _tmp = &nbrs[j];
                    for (auto u: elem){
                        if (u!=v_id){
                            _tmp->insert(node_index[u]);
                        }
                    }
                }
            }
        }
		for (size_t u = 0; u<N; u++){
            if(pcore[u]>=pk){
                auto d = inc_edges[u].size();
                if (degbucket.find(d) == degbucket.end()) degbucket[d] = uintSet();
                degbucket[d].insert(u);
                inverse_bucket[u] = d;
                max_deg = std::max(d,max_deg);
            }
		}
        if(log) {
            std::cout<<"max_deg: "<<max_deg<<"\n";
            std::cout<<"init degbucket: \n"; 
            print_bucket(degbucket,init_nodes);
            std::cout<<"inc_edges: \n";
            for(size_t i =0; i<N; i++){
                std::cout<<init_nodes[i]<<": ";
                for(auto u: inc_edges[i])  std::cout<<init_nodes[u]<<", ";
                std::cout<<"\n";
            }
        }
		// bool stop = false;
		// size_t maximal_dk = 1;
		for(size_t dk = 1; dk<= max_deg; dk++){
            if(log) std::cout<<"dk = "<<dk<<"\n";
            if(degbucket.find(dk)==degbucket.end()) continue;
			while (degbucket[dk].size()!=0){
				// Pop v from degbucket[dk];
                auto set_it = degbucket[dk].begin();  //# get first element in the bucket
                auto v = *set_it;
                degbucket[dk].erase(set_it);
                if (log){
                    std::cout<<"pop: "<<v<<"/"<<init_nodes[v]<<"\n";
                }
				// score[v] = dk; // assign secondary core-num to v
                a.kdcores.push_back(std::make_tuple(init_nodes[v],pk,dk));
                // stop[v] = true;
                intvec nbrs_v;
                iterate_nbrs(v, nbrs_v, inc_edges, edges);
                if(log){
                    std::cout<<"iterate_nbrs: \n";
                    for(auto u: nbrs_v) std::cout<<u<<"/"<<init_nodes[u]<<",";
                    std::cout<<"\n";
                }
                removeV_transform(v,inc_edges, edges);
                if(log){
                    std::cout<<"inc_edge after removal: \n";
                    for(size_t i =0; i<N; i++){
                        std::cout<<init_nodes[i]<<": ";
                        for(auto u: inc_edges[i])  std::cout<<init_nodes[u]<<", ";
                        std::cout<<"\n";
                    }
                    std::cout<<"done\n";
                }
				for(auto u: nbrs_v){
                    if(log) std::cout<<" -- "<<u<<"/"<<init_nodes[u]<<"\n"; 
                    // if(stop[u]) continue;
					// if |N(u)| in residual hyp < primary core , stop 
					if (get_number_of_nbrs(u, inc_edges, edges)< pk){
						// stop  peeling v caused nbr u's |N(u)| in the current subhyp. < pk
                        if(log) std::cout<<"stop\n"; 
                        if(log){
                            std::cout<<"current deg: \n";
                            for(size_t i = 0; i<N; i++){
                                std::cout<<init_nodes[i]<<": "<<inc_edges[i].size()<<"\n";
                            }
                        }
                        degbucket[inverse_bucket[u]].erase(u); // erase u from previous bucket index
                        // if (degbucket.find(dk) == degbucket.end()) degbucket[dk] = uintSet();
                        degbucket[dk].insert(u);
                        inverse_bucket[u] = dk;
                        if (log) {std::cout<<"bucket: \n"; print_bucket(degbucket,init_nodes);}
                        if(log) std::cout<<"batch delete\n";
                        /* We peel remaining nodes with pcore[u] == pk one by one without 
                        doing expensive nbr traversal for efficiency. 
                        We could have taken induced subhypergrpah {u: pcore[u]>=pk} but that would 
                        require constructing sub-hyp. from scratch which is again more expensive than just
                        peeling the remainder nodes.
                        */
                        // for(size_t ddk = dk; ddk<=max_deg; ddk++){
                        //     if (log) std::cout<<"ddk: "<<ddk<<"\n";
                        //     if(degbucket.find(ddk)!=degbucket.end()){
                        //         while(degbucket[ddk].size()){
                        //             auto set_it = degbucket[ddk].begin();  //# get first element in the bucket
                        //             auto u = *set_it;
                        //             if (log) std::cout<<init_nodes[u]<<",";
                        //             degbucket[ddk].erase(set_it);
                        //             if(pcore[u]==pk)
                        //                 removeV_transform(u,inc_edges, edges);
                        //             score[u] = dk;
                        //         }
                        //     }
                        //     if(log) std::cout<<"\n";
                        // }
                        // break;
					}
					else{ // else, update index in degree bucket for u \in N(v) 
                        // only update bucket position of nodes in nbr pk-core.
                        // pk+1, and higher core-nodes will be processed in later time.
                        // if (nbrbucket[pk].find(u) != nbrbucket[pk].end()){
                            auto d = inc_edges[u].size();
                            d = std::max(d,dk);
                            degbucket[inverse_bucket[u]].erase(u); // erase u from previous bucket index
                            if (degbucket.find(d) == degbucket.end()) degbucket[d] = uintSet();
                            degbucket[d].insert(u);
                            if (log){std::cout<< "bucket update: \n";    print_bucket(degbucket,init_nodes);}
                            inverse_bucket[u] = d;
                        // }
					}
				}
                if (log) std::cout<<"done traversing nbrs\n";
			}
		}
    }
    a.exec_time = double(clock() - start) / double(CLOCKS_PER_SEC);
    a.output["execution time"]= std::to_string(a.exec_time);
    a.output["total iteration"] = std::to_string(iterations);
    a.output["init_time"] = std::to_string(init_tm);

}
