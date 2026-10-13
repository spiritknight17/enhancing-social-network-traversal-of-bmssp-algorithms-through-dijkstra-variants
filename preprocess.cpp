#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <random>
#include <cmath>
#include <chrono>
using namespace std;

struct RawEdge { uint32_t u, v; };

void process_dataset(const string& name, const string& txt_path, const string& base_dir) {
    cout << "\n=======================================================\n";
    cout << "Processing Topology: " << name << "\n";
    cout << "=======================================================\n";

    ifstream infile(txt_path);
    if (!infile) { cout << "ERROR: Could not find " << txt_path << "\n"; return; }

    string line;
    vector<RawEdge> edges;
    uint32_t max_node_id = 0;

    cout << "1. Parsing raw text file...\n";
    while (getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line); uint32_t u, v;
        if (ss >> u >> v) {
            if (u == v) continue;
            edges.push_back({u, v});
            max_node_id = max({max_node_id, u, v});
        }
    }
    infile.close();
    uint32_t num_nodes = max_node_id + 1;

    cout << "2. Sorting edges for CSR (structural preparation)...\n";
    sort(edges.begin(), edges.end(), [](const RawEdge& a, const RawEdge& b){
        return a.u != b.u ? a.u < b.u : a.v < b.v;
    });
    edges.erase(unique(edges.begin(), edges.end(), [](const RawEdge& a, const RawEdge& b){
        return a.u == b.u && a.v == b.v;
    }), edges.end());
    uint32_t num_edges = (uint32_t)edges.size();
    cout << "   -> " << num_nodes << " Nodes | " << num_edges << " Edges\n";

    vector<uint32_t> indptr(num_nodes+1,0), indices(num_edges),
                     w_diffusion(num_edges);
    mt19937 rng(42);
    uniform_real_distribution<double>  prob_dist(0.01,0.99);

    cout << "3. Assigning diffusion edge weights (Information Distance w=-log(p)*10000)...\n";
    for (uint32_t i = 0; i < num_edges; ++i) {
        indptr[edges[i].u + 1]++;
        indices[i] = edges[i].v;
        double p = prob_dist(rng);
        w_diffusion[i] = max(1u, (uint32_t)(-log(p)*10000.0));
    }
    for (uint32_t i = 0; i < num_nodes; ++i) indptr[i+1] += indptr[i];

    auto write_bin = [&](const string& type, const vector<uint32_t>& w){
        string path = base_dir+"/"+name+"_"+type+".bin";
        ofstream out(path, ios::binary);
        uint32_t hdr[2]={num_nodes,num_edges};
        out.write((char*)hdr,8);
        out.write((char*)indptr.data(),(num_nodes+1)*4);
        out.write((char*)indices.data(),num_edges*4);
        out.write((char*)w.data(),num_edges*4);
        cout << "   -> Saved: " << path << "\n";
    };
    cout << "4. Writing binary CSR file to Google Drive (diffusion only)...\n";
    write_bin("diffusion",  w_diffusion);
}

int main() {
    auto t0 = chrono::high_resolution_clock::now();
    string base_dir = getenv("GRAPH_DIR") ? getenv("GRAPH_DIR") : "/content/drive/MyDrive/GraphDatasets";
    process_dataset("soc_pokec",   base_dir+"/soc_pokec.txt",   base_dir);
    process_dataset("livejournal", base_dir+"/livejournal.txt", base_dir);
    process_dataset("orkut",       base_dir+"/orkut.txt",       base_dir);
    chrono::duration<double> dt = chrono::high_resolution_clock::now()-t0;
    cout << "\nPreprocessing complete in " << dt.count() << "s\n";
    cout << "Diffusion .bin files saved to Google Drive (structural CSR is internal preparation).\n";
    return 0;
}
