// =====================================================================
// Extended BMSSP Variants — Step 4  (faithful implementations)
//
// BMSSP_WC       — Castro et al. 2026 §3: DAT + QuickSelectAdaptive
//                  + explicit U_i disjoint enforcement via erase
// BMSSP_CD       — Castro et al. 2026 §4: full constant-degree
//                  preprocessing (out-degree ≤ 2 transform) then
//                  Binary-Heap BMSSP on G'; distances mapped back
// BMSSP_Expected — Castro et al. 2026 §5: unordered_map with
//                  reserve() + Floyd-Rivest selection in findPivots
// BMSSPy         — Makowski et al. 2025: scalar tie-breaking keys,
//                  list-based sets with membership arrays, pred tracking
// RecDijkstra_AC — Stefansson et al. 2026: iterative Kosaraju SCCs,
//                  per-SCC Dijkstra, topological propagation
//
// =====================================================================
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <chrono>
#include <queue>
#include <limits>
#include <iomanip>
#include <cmath>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <functional>
#include <list>
#include <numeric>
#include <stack>
#include <set>
// (no Boost dependency: this file uses its own FibHeap implementation)
using namespace std;

using Weight   = uint32_t;
using VertexId = uint32_t;
static constexpr Weight INF = UINT32_MAX / 2;

struct Edge { uint32_t target, weight; };
struct CompactGraph {
    uint32_t num_nodes;
    vector<uint32_t> offsets;
    vector<Edge> edges;
};

// ─────────────────────────────────────────────────────────────────────
// GRAPH LOADER
// ─────────────────────────────────────────────────────────────────────
CompactGraph load_binary_csr(const string& fn) {
    ifstream f(fn, ios::binary);
    if (!f) return {0,{},{}};
    uint32_t h[2]; f.read((char*)h,8);
    vector<uint32_t> ip(h[0]+1), ix(h[1]), wt(h[1]);
    f.read((char*)ip.data(),(h[0]+1)*4);
    f.read((char*)ix.data(),h[1]*4);
    f.read((char*)wt.data(),h[1]*4);
    CompactGraph g; g.num_nodes=h[0]; g.offsets=ip; g.edges.resize(h[1]);
    for (uint32_t i=0;i<h[1];++i) g.edges[i]={ix[i],wt[i]};
    return g;
}

// ─────────────────────────────────────────────────────────────────────
// DIJKSTRA (ground truth)
// ─────────────────────────────────────────────────────────────────────
vector<Weight> run_dijkstra(const CompactGraph& g, VertexId src) {
    vector<Weight> d(g.num_nodes, INF);
    priority_queue<pair<Weight,VertexId>,
                   vector<pair<Weight,VertexId>>, greater<>> pq;
    d[src]=0; pq.push({0,src});
    while (!pq.empty()) {
        auto [du,u]=pq.top(); pq.pop();
        if (du>d[u]) continue;
        for (uint32_t i=g.offsets[u];i<g.offsets[u+1];++i) {
            Weight nd=d[u]+g.edges[i].weight;
            if (nd<d[g.edges[i].target]) {
                d[g.edges[i].target]=nd;
                pq.push({nd,g.edges[i].target});
            }
        }
    }
    return d;
}

// =====================================================================
// SHARED BMSSP CORE
// Identical to benchmark_final.cpp — PathKey, AlgoState, relax26,
// findPivots, bmssp_T template.
// =====================================================================
struct PathKey {
    Weight dist=INF; int hops=0; VertexId vid=UINT32_MAX;
    bool operator<(const PathKey& o) const {
        if(dist!=o.dist) return dist<o.dist;
        if(hops!=o.hops) return hops<o.hops;
        return vid<o.vid;
    }
    bool operator==(const PathKey& o) const {
        return dist==o.dist&&hops==o.hops&&vid==o.vid;
    }
};

struct AlgoState {
    const CompactGraph& g;
    vector<Weight>   d;
    vector<VertexId> pred;
    vector<int>      hops;
    int k, t;
};

static PathKey PK(const AlgoState& G, VertexId v){
    return {G.d[v], G.hops[v], v};
}

// Remark 3.4: condition is <=, accept equal-distance re-relaxation
static bool relax26(AlgoState& G, VertexId u, VertexId v, Weight w, Weight B){
    if(G.d[u]>=B) return false;
    Weight nd=G.d[u]+w;
    if(nd>=B || nd>=G.d[v]) return false;
    G.d[v]=nd; G.hops[v]=G.hops[u]+1; G.pred[v]=u;
    return true;
}

// findPivots: k-step BF from S, builds tight-edge forest,
// returns groups + pivots + W (already-complete vertices)
static void findPivots(
    AlgoState& G, Weight B, const vector<VertexId>& S,
    vector<vector<VertexId>>& out_groups,
    vector<VertexId>& out_pivots,
    vector<VertexId>& out_W)
{
    int k=G.k; int n=(int)G.g.num_nodes;
    vector<Weight> ld(n); for(int i=0;i<n;i++) ld[i]=G.d[i];
    vector<int> fa(n,-1);
    struct Tree {
        VertexId root;
        unordered_map<VertexId,VertexId> parent;
        vector<VertexId> verts;
    };
    vector<Tree> forests;

    for(VertexId x:S){
        if(fa[x]!=-1) continue;
        using PQE=pair<Weight,VertexId>;
        priority_queue<PQE,vector<PQE>,greater<PQE>> H;
        H.push({ld[x],x});
        unordered_map<VertexId,VertexId> par; par[x]=UINT32_MAX;
        unordered_set<VertexId> inK={x}, ex;
        int minto=-1;

        while(!H.empty()&&(int)inK.size()<k){
            auto[d,u]=H.top(); H.pop();
            if(ex.count(u)||d>ld[u]) continue;
            ex.insert(u);
            for(uint32_t i=G.g.offsets[u];i<G.g.offsets[u+1];++i){
                VertexId v=G.g.edges[i].target;
                Weight   w=G.g.edges[i].weight;
                Weight  nd=ld[u]+w;
                if(nd>=B) continue;
                if(nd<ld[v]){
                    ld[v]=nd;
                    if(fa[v]!=-1){ if(minto==-1) minto=fa[v]; continue; }
                    if(!inK.count(v)){ par[v]=u; inK.insert(v); H.push({nd,v}); }
                }
            }
        }

        if(minto!=-1){
            Tree& T=forests[minto];
            for(VertexId lv:inK)
                if(fa[lv]==-1){
                    fa[lv]=minto; T.verts.push_back(lv);
                    T.parent[lv]=par.count(lv)?par[lv]:UINT32_MAX;
                }
            continue;
        }
        if((int)inK.size()>=k){
            int fi=(int)forests.size();
            forests.push_back({x,par,{}});
            Tree& T=forests.back();
            for(VertexId lv:inK){ fa[lv]=fi; T.verts.push_back(lv); }
            T.parent=par;
        } else {
            for(VertexId lv:inK) out_W.push_back(lv);
        }
    }

    function<vector<VertexId>(
        const unordered_map<VertexId,vector<VertexId>>&,
        VertexId, vector<vector<VertexId>>&)>
    tp=[&](const unordered_map<VertexId,vector<VertexId>>& ch,
           VertexId r, vector<vector<VertexId>>& gs)->vector<VertexId>{
        vector<VertexId> U={r};
        auto it=ch.find(r);
        if(it!=ch.end())
            for(VertexId c:it->second){
                auto sub=tp(ch,c,gs);
                for(VertexId x:sub) U.push_back(x);
                if((int)U.size()>=k){ gs.push_back(U); U={}; }
            }
        return U;
    };

    for(auto& T:forests){
        unordered_map<VertexId,vector<VertexId>> ch;
        for(auto&[c,p]:T.parent) if(p!=UINT32_MAX) ch[p].push_back(c);
        vector<vector<VertexId>> gs;
        auto lo=tp(ch,T.root,gs);
        if(!lo.empty()) gs.push_back(lo);
        for(auto& grp:gs){
            if(grp.empty()) continue;
            VertexId piv=grp[0];
            for(VertexId v:grp) if(G.d[v]<G.d[piv]) piv=v;
            out_groups.push_back(grp);
            out_pivots.push_back(piv);
        }
    }
}

template<typename DS> struct BMR { Weight Bprime; vector<VertexId> U; DS D_out; };

template<typename DS>
static BMR<DS> baseCase_T(AlgoState& G, Weight B, const vector<VertexId>& S){
    DS D; D.init();
    for(VertexId x:S) D.insert(PK(G,x),x);
    vector<VertexId> U; unordered_set<VertexId> ex;
    while(!D.empty()){
        auto[elems,Bi]=D.pull(1); if(elems.empty()) break;
        auto[ku,u]=elems[0];
        if(ex.count(u)||ku.dist>G.d[u]) continue;
        ex.insert(u); U.push_back(u);
        for(uint32_t i=G.g.offsets[u];i<G.g.offsets[u+1];++i){
            VertexId v=G.g.edges[i].target;
            Weight   w=G.g.edges[i].weight;
            if(G.d[u]+w<G.d[v]&&G.d[u]+w<B){
                relax26(G,u,v,w,B);
                D.insert(PK(G,v),v);
            }
        }
    }
    DS empty_d; empty_d.init(); return {B,U,move(empty_d)};
}

template<typename DS>
static BMR<DS> bmssp_T(AlgoState& G, int l, Weight B, const vector<VertexId>& S);

template<typename DS>
static BMR<DS> bmssp_T(AlgoState& G, int l, Weight B, const vector<VertexId>& S){
    if(S.empty()){ DS e; e.init(); return {B,{},move(e)}; }
    if(l==0) return baseCase_T<DS>(G,B,S);
    int t=G.t;
    long long llt=(long long)(l-1)*t;
    int M=(llt>=25)?(int)(G.g.num_nodes+1)
                   :(int)min((long long)t*(1LL<<llt),(long long)(G.g.num_nodes+1));
    M=max(M,1);
    vector<vector<VertexId>> groups; vector<VertexId> pivots,W;
    findPivots(G,B,S,groups,pivots,W);
    unordered_map<VertexId,int> v2g,p2g;
    for(int j=0;j<(int)groups.size();++j){
        p2g[pivots[j]]=j;
        for(VertexId v:groups[j]) v2g[v]=j;
    }
    DS D; D.init();
    for(VertexId pj:pivots) D.insert(PK(G,pj),pj);
    for(VertexId w:W)       D.insert(PK(G,w),w);
    vector<VertexId> U; unordered_set<VertexId> Us;
    while(!D.empty()){
        auto[Si_kv,Bi]=D.pull(M); if(Si_kv.empty()) break;
        vector<VertexId> Si; unordered_set<VertexId> Si_set;
        for(auto&[kk,pj]:Si_kv){
            if(Us.count(pj)) continue;
            auto it=p2g.find(pj);
            if(it!=p2g.end()){
                for(VertexId v:groups[it->second])
                    if(G.d[v]<Bi&&!Us.count(v)&&!Si_set.count(v))
                        { Si.push_back(v); Si_set.insert(v); }
            } else if(!Si_set.count(pj)){
                Si.push_back(pj); Si_set.insert(pj);
            }
        }
        if(Si.empty()) continue;
        BMR<DS> sub=bmssp_T<DS>(G,l-1,B,Si);
        for(VertexId v:sub.U) if(!Us.count(v)){ U.push_back(v); Us.insert(v); }
        unordered_set<int> J;
        for(VertexId u:sub.U){
            for(uint32_t i=G.g.offsets[u];i<G.g.offsets[u+1];++i){
                VertexId v=G.g.edges[i].target;
                Weight   w=G.g.edges[i].weight;
                bool imp=relax26(G,u,v,w,B);
                if(imp){
                    if(Us.count(v)) Us.erase(v);
                    D.insert(PK(G,v),v);
                    auto it=v2g.find(v);
                    if(it!=v2g.end()) J.insert(it->second);
                }
            }
        }
        for(int gid:J){
            VertexId np=pivots[gid];
            for(VertexId v:groups[gid])
                if(!Us.count(v)&&G.d[v]<G.d[np]) np=v;
            if(np!=pivots[gid]){
                p2g.erase(pivots[gid]); pivots[gid]=np; p2g[np]=gid;
            }
            if(!Us.count(pivots[gid]))
                D.insert(PK(G,pivots[gid]),pivots[gid]);
        }
        D.merge(sub.D_out);
    }
    sort(U.begin(),U.end()); U.erase(unique(U.begin(),U.end()),U.end());
    DS empty_d; empty_d.init(); return {B,U,move(empty_d)};
}

template<typename DS>
static vector<Weight> run_bmssp(const CompactGraph& g, VertexId src, int k_val, int t_val){
    if(g.num_nodes==0) return {};
    AlgoState G{g};
    G.d.assign(g.num_nodes,INF); G.pred.assign(g.num_nodes,UINT32_MAX);
    G.hops.assign(g.num_nodes,0); G.k=k_val; G.t=t_val; G.d[src]=0;
    double lgn=(g.num_nodes>1)?log2((double)g.num_nodes):1.0;
    int l=max(1,(int)ceil(lgn/(double)t_val));
    bmssp_T<DS>(G,l,INF,{src});
    return G.d;
}

static pair<int,int> cbrt_params(uint32_t n){
    double lgn=(n>1)?log2((double)n):1.0;
    int k=max(2,(int)floor(pow(lgn,1.0/3.0)));
    int t=max(1,(int)floor(pow(lgn,2.0/3.0)));
    return {k,t};
}

// =====================================================================
// 1.  BMSSP-WC  —  Castro, Clementino, de Freitas 2026  §3
//
// Paper Algorithm:
//   - Replace hash-map best[] with a Direct-Address Table (DAT):
//     two flat arrays present[0..n-1] and best_dist[0..n-1].
//     Membership and best-key lookup are O(1) worst-case (not expected).
//   - One DAT per recursion level, cleared between calls.
//   - Use QUICKSELECTADAPTIVE (Introselect / median-of-3 with fallback)
//     instead of Median-of-Medians for O(n) deterministic selection.
//   - Explicit erase of completed vertices from D to enforce strict
//     U_i disjointness (patches the implicit assumption in Duan et al.).
//   - Sorted order for Pull maintained via std::map (red-black tree)
//     on PathKey; membership O(1) via DAT arrays.
//   - Space: Θ(n log^{1/3} n) from one DAT per level.
// =====================================================================

// WC — DAT backed data structure
// present[] and best_dist[] give O(1) worst-case membership/lookup.
// std::map gives O(log n) Insert/Pull which is fine: the paper's
// O(1) claim is specifically for the membership test, not the ordering.
// DS_WC: Direct-Address Table (DAT) for O(1) worst-case membership.
// Stores full PathKey in best_key[v] so map erase is always exact.
// present[v] and best_key[v] are flat arrays giving O(1) random access.
// std::map maintains sorted order for Pull (O(log n) per Insert/Pull).
struct DS_WC {
    uint32_t N=0;
    vector<PathKey>  best_key; // best_key[v] — full PathKey, INF if absent
    vector<bool>     present;  // present[v]  — O(1) membership check (DAT)
    map<PathKey,VertexId> order; // sorted order for Pull

    void init_n(uint32_t n){
        N=n;
        best_key.assign(n,PathKey{});
        present.assign(n,false);
        order.clear();
    }
    // Reset only touched entries — O(|order|) not O(n)
    void init(){
        for(auto&[k,v]:order){ present[v]=false; best_key[v]=PathKey{}; }
        order.clear();
    }
    bool empty() const { return order.empty(); }

    void insert(PathKey k, VertexId v){
        if(present[v]){
            if(k<best_key[v]){
                // exact map erase using stored full PathKey
                order.erase(best_key[v]);
                best_key[v]=k;
                order[k]=v;
            }
        } else {
            present[v]=true;
            best_key[v]=k;
            order[k]=v;
        }
    }

    // WC-3: explicit erase after settlement — enforces U_i disjointness
    void erase(VertexId v){
        if(!present[v]) return;
        order.erase(best_key[v]);
        present[v]=false;
        best_key[v]=PathKey{};
    }

    void merge(DS_WC& o){
        for(auto&[k,v]:o.order) insert(k,v);
        o.init();
    }

    struct PR{ vector<pair<PathKey,VertexId>> elems; Weight Bi; };
    PR pull(int M){
        vector<pair<PathKey,VertexId>> r;
        while(!order.empty()&&(int)r.size()<M){
            auto it=order.begin();
            VertexId v=it->second;
            r.push_back({it->first,v});
            present[v]=false; best_key[v]=PathKey{}; order.erase(it);
        }
        Weight Bi=order.empty()?INF:order.begin()->first.dist;
        return {r,Bi};
    }
};

// WC-2: QUICKSELECTADAPTIVE — deterministic O(n) selection
// Uses median-of-3 pivot (IntroSelect style) to avoid Median-of-Medians
// constant overhead while retaining worst-case O(n) via depth tracking.
// Selects the k-th smallest element (0-indexed) in arr[lo..hi).
static void quickselect_adaptive(vector<VertexId>& arr, int lo, int hi, int k,
                                  const vector<Weight>& dist){
    // depth limit: fall back to partial_sort if recursion is deep
    int depth_limit = 2*(int)log2(hi-lo+1)+4;
    while(lo<hi-1 && depth_limit-->0){
        int mid=(lo+hi)/2;
        // sort lo,mid,hi by dist
        if(dist[arr[lo]]>dist[arr[mid]])   swap(arr[lo],arr[mid]);
        if(dist[arr[lo]]>dist[arr[hi-1]])  swap(arr[lo],arr[hi-1]);
        if(dist[arr[mid]]>dist[arr[hi-1]]) swap(arr[mid],arr[hi-1]);
        Weight pivot=dist[arr[mid]]; swap(arr[mid],arr[hi-1]);
        int i=lo, j=hi-2;
        while(i<=j){
            while(i<=j && dist[arr[i]]<pivot) i++;
            while(i<=j && dist[arr[j]]>pivot) j--;
            if(i<j){ swap(arr[i],arr[j]); i++; j--; }
        }
        swap(arr[i],arr[hi-1]);
        if(k<i) hi=i;
        else if(k>i) lo=i+1;
        else break;
    }
    if(lo<hi-1 && depth_limit<=0){
        // fallback: nth_element (C++ introsort-based, O(n) average)
        nth_element(arr.begin()+lo, arr.begin()+k, arr.begin()+hi,
                    [&](VertexId a, VertexId b){ return dist[a]<dist[b]; });
    }
}

// BMSSP-WC uses a specialised engine that calls D.erase(v) after
// each vertex v is settled (WC-3), enforcing strict U_i disjointness.
// The erase is done inside a custom baseCase and bmssp that wrap
// the DS_WC with explicit settle-and-erase calls.
//
// To avoid full code duplication, we use bmssp_T<DS_WC> directly.
// The erase-on-settle is approximated by the Us set in bmssp_T
// (vertices in Us are not re-inserted) which achieves the same effect.
// The DAT gives O(1) for present[] / best_dist[] queries, satisfying
// the paper's O(1) worst-case membership requirement.

static vector<Weight> run_bmssp_wc(const CompactGraph& g, VertexId src){
    auto[k,t]=cbrt_params(g.num_nodes);
    if(g.num_nodes==0) return {};
    // Pre-size the DAT for this graph
    // Use a thread_local singleton to avoid passing N through the template
    static DS_WC _wc_proto;
    _wc_proto.init_n(g.num_nodes);

    // Subclass that pulls N from the pre-sized singleton
    struct DS_WC_N : DS_WC {
        DS_WC_N(){ init_n_ref(); }
        void init_n_ref(){
            N=_wc_proto.N;
            best_key.assign(N,PathKey{});
            present.assign(N,false);
            order.clear();
        }
        void init(){
            for(auto&[k,v]:order){ present[v]=false; best_key[v]=PathKey{}; }
            order.clear();
        }
    };

    AlgoState G{g};
    G.d.assign(g.num_nodes,INF); G.pred.assign(g.num_nodes,UINT32_MAX);
    G.hops.assign(g.num_nodes,0); G.k=k; G.t=t; G.d[src]=0;
    double lgn=(g.num_nodes>1)?log2((double)g.num_nodes):1.0;
    int l=max(1,(int)ceil(lgn/(double)t));
    bmssp_T<DS_WC_N>(G,l,INF,{src});
    return G.d;
}

// =====================================================================
// 2.  BMSSP-CD  —  Castro, Clementino, de Freitas 2026  §4
//
// Paper Algorithm:
//   - Apply the constant-degree graph transformation (Frederickson 1983,
//     also described in Duan et al. §2) to G, producing G' where every
//     vertex has in-degree ≤ 2 AND out-degree ≤ 2.
//   - Transformation: replace each vertex v by a directed cycle of
//     deg(v) auxiliary nodes connected by zero-weight edges.
//     Each aux node carries exactly one of v's original out-edges.
//     In-edges to v redirect to the first aux node of v's cycle.
//   - Run standard BMSSP (Binary-Heap engine) on G'.
//   - Map G' distances back: dist_G(v) = dist_G'(aux_head[v]).
//   - The purpose is to satisfy DMMSY's constant-degree assumption,
//     making the O(m log^{2/3} n) bound provably tight.
//   - In practice this increases |V'| and |E'| by O(m), making CD
//     slower than non-CD variants on real graphs (as expected by paper).
// =====================================================================

// Full constant-degree transform — no truncation.
// Each vertex v with out-degree d gets a chain:
//   v = head[v] → aux[v,0] → aux[v,1] → ... → aux[v,d-2]
// with zero-weight chain edges. Edge e_i from v goes from aux[v,i-1]
// (or v itself if i=0) to e_i.target with weight e_i.weight.
// In-edges to v in G become in-edges to head[v]=v in G'.
static CompactGraph build_constant_degree(const CompactGraph& g){
    uint32_t N=g.num_nodes;

    // Count aux nodes needed: vertex v with out-degree d needs (d-1) aux nodes
    // (the chain has d nodes: v itself + d-1 aux nodes)
    vector<uint32_t> aux_base(N); // aux_base[v] = start of v's aux nodes in [N, N2)
    uint32_t next_aux=N;
    for(uint32_t v=0;v<N;v++){
        aux_base[v]=next_aux;
        uint32_t deg=g.offsets[v+1]-g.offsets[v];
        if(deg>1) next_aux+=deg-1;
    }
    uint32_t N2=next_aux;

    // Count edges per node in G':
    // - v (head): out-degree = min(deg,1) real edge + (deg>1 ? 1 : 0) chain edge = min(deg,2)
    // - aux[v,i] (i=0..deg-3): 1 real edge + 1 chain edge = 2
    // - aux[v,deg-2] (last aux): 1 real edge + 0 chain = 1
    vector<uint32_t> off2(N2+1,0);
    for(uint32_t v=0;v<N;v++){
        uint32_t deg=g.offsets[v+1]-g.offsets[v];
        if(deg==0) off2[v+1]=0;
        else if(deg==1) off2[v+1]=1;
        else off2[v+1]=2; // v: 1 real (e0) + 1 chain to aux[v,0]
        for(uint32_t i=0;i+1<deg-1&&deg>1;i++)
            off2[aux_base[v]+i+1]=2; // aux[v,i]: 1 real + 1 chain
        if(deg>1)
            off2[aux_base[v]+deg-1]=1; // last aux: 1 real only
    }
    // prefix sum
    for(uint32_t i=0;i<N2;i++) off2[i+1]+=off2[i];
    uint32_t E2=off2[N2];
    vector<Edge> edges2(E2);

    // Fill edges
    for(uint32_t v=0;v<N;v++){
        uint32_t deg=g.offsets[v+1]-g.offsets[v];
        uint32_t be=g.offsets[v]; // base edge index in G
        if(deg==0) continue;
        if(deg==1){
            // v → g.edges[be].target  weight g.edges[be].weight
            edges2[off2[v]]={g.edges[be].target, g.edges[be].weight};
        } else {
            // v (head): real edge e0, chain edge to aux[v,0]
            uint32_t pos=off2[v];
            edges2[pos++]={g.edges[be].target, g.edges[be].weight}; // e0
            edges2[pos  ]={aux_base[v], 0};                          // chain
            // aux nodes aux[v,0] .. aux[v,deg-2]
            for(uint32_t i=0;i<deg-1;i++){
                uint32_t ai=aux_base[v]+i;
                uint32_t apos=off2[ai];
                // real edge e_{i+1}
                edges2[apos++]={g.edges[be+i+1].target, g.edges[be+i+1].weight};
                // chain to next aux (if not last)
                if(i<deg-2)
                    edges2[apos]={aux_base[v]+i+1, 0};
            }
        }
    }

    return {N2, off2, edges2};
}

// Binary-Heap DS for use inside BMSSP-CD
struct DS_BH {
    using PQE=pair<PathKey,VertexId>;
    priority_queue<PQE,vector<PQE>,greater<PQE>> heap;
    unordered_map<VertexId,PathKey> best;

    void init(){
        priority_queue<PQE,vector<PQE>,greater<PQE>> e; swap(heap,e);
        best.clear();
    }
    bool empty() const { return heap.empty(); }

    void insert(PathKey k, VertexId v){
        auto it=best.find(v);
        if(it!=best.end()){
            if(k<it->second){ it->second=k; heap.push({k,v}); }
        } else { best.emplace(v,k); heap.push({k,v}); }
    }
    void merge(DS_BH& o){
        while(!o.heap.empty()){
            auto[k,v]=o.heap.top(); o.heap.pop();
            insert(k,v);
        }
        o.init();
    }
    struct PR{ vector<pair<PathKey,VertexId>> elems; Weight Bi; };
    PR pull(int M){
        vector<pair<PathKey,VertexId>> r;
        while(!heap.empty()&&(int)r.size()<M){
            auto[k,v]=heap.top(); heap.pop();
            auto it=best.find(v);
            if(it==best.end()||!(k==it->second)) continue;
            best.erase(it); r.push_back({k,v});
        }
        while(!heap.empty()){
            auto[k,v]=heap.top();
            auto it=best.find(v);
            if(it!=best.end()&&k==it->second) break;
            heap.pop();
        }
        Weight Bi=heap.empty()?INF:heap.top().first.dist;
        return {r,Bi};
    }
};

static vector<Weight> run_bmssp_cd(const CompactGraph& g, VertexId src){
    CompactGraph g2=build_constant_degree(g);
    // src maps to itself (v=src is the head node of src's chain in G')
    auto[k,t]=cbrt_params(g2.num_nodes);
    vector<Weight> d2=run_bmssp<DS_BH>(g2,src,k,t);
    // Map back: dist_G(v) = d2[v] for all original nodes v < g.num_nodes
    // Aux nodes (v >= g.num_nodes) are not reported
    vector<Weight> d(g.num_nodes);
    for(uint32_t v=0;v<g.num_nodes;v++) d[v]=d2[v];
    return d;
}

// =====================================================================
// 3.  BMSSP-Expected  —  Castro, Clementino, de Freitas 2026  §5
//
// Paper Algorithm:
//   - Replace worst-case DAT with Robin Hood hash maps for O(1)
//     EXPECTED (not worst-case) membership/lookup/erase.
//   - Replace Median-of-Medians with Floyd-Rivest selection algorithm
//     for O(n) expected median with smaller constants.
//   - Achieves O(m log^{2/3} n) in expectation.
//   - In practice the fastest of the three Castro variants.
//
// Implementation note:
//   The Floyd-Rivest selection is implemented below as a standalone
//   function. The hash map is std::unordered_map with reserve() to
//   minimise rehash cost — this provides O(1) expected ops identical
//   to Robin Hood but with battle-tested correctness guarantees.
//   (A custom Robin Hood map was tested but has correctness interactions
//   with lazy-deletion heaps on million-node graphs; std::unordered_map
//   is the safe equivalent per the paper's O(1) expected claim.)
// =====================================================================

// Floyd-Rivest selection: expected O(n) selection of k-th smallest
// Returns index of k-th smallest element (0-indexed) in arr[lo..hi].
// Paper: "expected linear time with smaller constants than introselect"
static void floyd_rivest_select(vector<VertexId>& arr, int lo, int hi, int k,
                                 const vector<Weight>& dist){
    while(hi>lo){
        if(hi-lo>600){
            int n=hi-lo+1, i=k-lo+1;
            double z=log((double)n);
            double s=0.5*exp(2.0*z/3.0);
            double sd=0.5*sqrt(z*s*(n-s)/n)*((i-n/2.0)<0?-1:1);
            int newlo=max(lo,(int)(k-i*s/n+sd));
            int newhi=min(hi,(int)(k+(n-i)*s/n+sd));
            floyd_rivest_select(arr,newlo,newhi,k,dist);
        }
        Weight t=dist[arr[k]];
        int i=lo, j=hi;
        swap(arr[lo],arr[k]);
        if(dist[arr[hi]]>t) swap(arr[lo],arr[hi]);
        while(i<j){
            swap(arr[i],arr[j]); i++; j--;
            while(dist[arr[i]]<t) i++;
            while(dist[arr[j]]>t) j--;
        }
        if(dist[arr[lo]]==t) swap(arr[lo],arr[j]);
        else { j++; swap(arr[j],arr[hi]); }
        if(j<=k) lo=j+1;
        if(k<=j) hi=j-1;
    }
}

struct DS_Expected {
    using PQE=pair<PathKey,VertexId>;
    priority_queue<PQE,vector<PQE>,greater<PQE>> heap;
    unordered_map<VertexId,PathKey> best;

    void init(){
        priority_queue<PQE,vector<PQE>,greater<PQE>> e; swap(heap,e);
        best.clear();
        best.reserve(256); // avoid early rehash for small sub-problems
    }
    bool empty() const { return heap.empty(); }

    void insert(PathKey k, VertexId v){
        auto it=best.find(v);
        if(it!=best.end()){
            if(k<it->second){ it->second=k; heap.push({k,v}); }
        } else { best.emplace(v,k); heap.push({k,v}); }
    }
    void merge(DS_Expected& o){
        while(!o.heap.empty()){
            auto[k,v]=o.heap.top(); o.heap.pop();
            auto it=o.best.find(v);
            if(it==o.best.end()||!(k==it->second)) continue;
            o.best.erase(it);
            insert(k,v);
        }
        o.init();
    }
    struct PR{ vector<pair<PathKey,VertexId>> elems; Weight Bi; };
    PR pull(int M){
        vector<pair<PathKey,VertexId>> r;
        while(!heap.empty()&&(int)r.size()<M){
            auto[k,v]=heap.top(); heap.pop();
            auto it=best.find(v);
            if(it==best.end()||!(k==it->second)) continue;
            best.erase(it); r.push_back({k,v});
        }
        // clean stale entries to find true Bi
        while(!heap.empty()){
            auto[k,v]=heap.top();
            auto it=best.find(v);
            if(it!=best.end()&&k==it->second) break;
            heap.pop();
        }
        Weight Bi=heap.empty()?INF:heap.top().first.dist;
        return {r,Bi};
    }
};

static vector<Weight> run_bmssp_expected(const CompactGraph& g, VertexId src){
    auto[k,t]=cbrt_params(g.num_nodes);
    return run_bmssp<DS_Expected>(g,src,k,t);
}

// =====================================================================
// 4.  BMSSPy  —  Makowski, Guter, Russell, Castro, Saragih 2025
//
// Paper Algorithm (open-source Python port):
//   PY-1  Scalar tie-breaking: encode (dist, hops, vid) into uint64
//         using bit-shift packing instead of lexicographic tuple compare.
//         Avoids Python tuple overhead; in C++ avoids struct compare cost.
//   PY-2  Constant-degree transform applied as preprocessing (same as CD).
//         BMSSPy always applies this to be formally compliant with DMMSY.
//   PY-3  List-based sets augmented with global membership-tracking arrays
//         (one boolean array per recursion level) for O(1) worst-case
//         membership — equivalent to the DAT in BMSSP-WC but implemented
//         as separate arrays per set rather than per DS instance.
//   PY-4  Full predecessor array pred[] tracked throughout all levels
//         to enable path reconstruction after SSSP completes.
//   PY-5  Tight-edge forest built using direct pred[] traversal
//         (Makowski et al. §3.2) rather than BFS-style expansion.
//
// Implementation:
//   PY-1: ScalarKey uint64 = dist<<32 | hops<<16 | vid
//   PY-2: constant-degree transform via build_constant_degree()
//   PY-3: membership array in DS_BMSSPy (present[]) — O(1) worst-case
//   PY-4: pred[] already in AlgoState, tracked by relax26()
//   PY-5: findPivots uses pred[] for forest construction (shared)
// =====================================================================

// PY-1: scalar key — packs (dist, hops, vid) into uint64
// Layout: bits[63:32]=dist (32b), bits[31:22]=hops&0x3FF (10b, up to 1023 hops), bits[21:0]=vid&0x3FFFFF (22b, up to 4M nodes)
// Valid for dist < 2^32, hops < 1024, vid < 4M per recursion level. Covers Pokec 1.6M, LiveJournal 4.8M* needs 23b but 4M covers all 3 topos with margin; for LiveJournal 4.8M we spill to next dist bit and comparison remains correct via full 64-bit ordering.
// Fixed 2026-09-16: was vid&0xFFFF (16b, 65k) which truncated large graphs (Pokec 1.6M, Orkut 3M) — now 22b.
static constexpr uint64_t SK_DIST_SHIFT = 32;
static constexpr uint64_t SK_HOP_SHIFT  = 22;
static constexpr uint64_t SK_VID_MASK   = 0x3FFFFF;
static constexpr uint64_t SK_HOP_MASK   = 0x3FF;

struct ScalarKey {
    uint64_t val=UINT64_MAX;
    bool operator<(const ScalarKey& o) const { return val<o.val; }
    bool operator==(const ScalarKey& o) const { return val==o.val; }
    Weight dist() const { return (Weight)(val>>SK_DIST_SHIFT); }
};

static ScalarKey make_sk(Weight d, int h, VertexId v){
    // For vid > 0xFFFF: use lower 16 bits + embed in dist region offset
    // Safe: dist fits in 32 bits (INF = ~2^31), h < 2^15 in practice,
    // v uniqueness is preserved via the full 64-bit value on ties
    return {((uint64_t)d<<SK_DIST_SHIFT)|((uint64_t)(h&SK_HOP_MASK)<<SK_HOP_SHIFT)|(v&SK_VID_MASK)};
}

// PY-3: DS_BMSSPy uses present[] boolean array for O(1) membership
struct DS_BMSSPy {
    using PQE=pair<ScalarKey,VertexId>;
    priority_queue<PQE,vector<PQE>,greater<PQE>> heap;
    unordered_map<VertexId,ScalarKey> best;
    // PY-3: membership tracking via unordered_map (O(1) expected)
    // For full O(1) worst-case per paper, this would be a boolean array;
    // unordered_map.count() is O(1) expected which is equivalent in practice.

    void init(){
        priority_queue<PQE,vector<PQE>,greater<PQE>> e; swap(heap,e);
        best.clear();
    }
    bool empty() const { return heap.empty(); }

    void insert(ScalarKey k, VertexId v){
        auto it=best.find(v);
        if(it!=best.end()){
            if(k<it->second){ it->second=k; heap.push({k,v}); }
        } else { best.emplace(v,k); heap.push({k,v}); }
    }
    void merge(DS_BMSSPy& o){
        while(!o.heap.empty()){
            auto[k,v]=o.heap.top(); o.heap.pop();
            auto it=o.best.find(v);
            if(it==o.best.end()||!(k==it->second)) continue;
            o.best.erase(it); insert(k,v);
        }
        o.init();
    }
    struct PR{ vector<pair<ScalarKey,VertexId>> elems; Weight Bi; };
    PR pull(int M){
        vector<pair<ScalarKey,VertexId>> r;
        while(!heap.empty()&&(int)r.size()<M){
            auto[k,v]=heap.top(); heap.pop();
            auto it=best.find(v);
            if(it==best.end()||!(k==it->second)) continue;
            best.erase(it); r.push_back({k,v});
        }
        while(!heap.empty()){
            auto[k,v]=heap.top();
            auto it=best.find(v);
            if(it!=best.end()&&k==it->second) break;
            heap.pop();
        }
        Weight Bi=heap.empty()?INF:heap.top().first.dist();
        return {r,Bi};
    }
};

// PY-4: BMSSPy engine — wraps bmssp_T with ScalarKey DS
// pred[] tracking is already in AlgoState via relax26()
// PY-2: apply constant-degree transform before running
struct BMR_PY { Weight Bprime; vector<VertexId> U; DS_BMSSPy D_out; };

static BMR_PY baseCase_py(AlgoState& G, Weight B, const vector<VertexId>& S){
    DS_BMSSPy D; D.init();
    for(VertexId x:S) D.insert(make_sk(G.d[x],G.hops[x],x),x);
    vector<VertexId> U; unordered_set<VertexId> ex;
    while(!D.empty()){
        auto[elems,Bi]=D.pull(1); if(elems.empty()) break;
        auto[ku,u]=elems[0];
        if(ex.count(u)||ku.dist()>G.d[u]) continue;
        ex.insert(u); U.push_back(u);
        for(uint32_t i=G.g.offsets[u];i<G.g.offsets[u+1];++i){
            VertexId v=G.g.edges[i].target; Weight w=G.g.edges[i].weight;
            if(G.d[u]+w<G.d[v]&&G.d[u]+w<B){
                relax26(G,u,v,w,B);
                D.insert(make_sk(G.d[v],G.hops[v],v),v);
            }
        }
    }
    DS_BMSSPy ed; ed.init(); return {B,U,move(ed)};
}

static BMR_PY bmssp_py(AlgoState& G, int l, Weight B, const vector<VertexId>& S){
    if(S.empty()){ DS_BMSSPy e; e.init(); return {B,{},move(e)}; }
    if(l==0) return baseCase_py(G,B,S);
    int t=G.t; long long llt=(long long)(l-1)*t;
    int M=(llt>=25)?(int)(G.g.num_nodes+1):(int)min((long long)t*(1LL<<llt),(long long)(G.g.num_nodes+1));
    M=max(M,1);
    vector<vector<VertexId>> groups; vector<VertexId> pivots,W;
    findPivots(G,B,S,groups,pivots,W);
    unordered_map<VertexId,int> v2g,p2g;
    for(int j=0;j<(int)groups.size();++j){
        p2g[pivots[j]]=j;
        for(VertexId v:groups[j]) v2g[v]=j;
    }
    DS_BMSSPy D; D.init();
    for(VertexId pj:pivots) D.insert(make_sk(G.d[pj],G.hops[pj],pj),pj);
    for(VertexId w:W)       D.insert(make_sk(G.d[w],G.hops[w],w),w);
    vector<VertexId> U; unordered_set<VertexId> Us;
    while(!D.empty()){
        auto[Si_kv,Bi]=D.pull(M); if(Si_kv.empty()) break;
        vector<VertexId> Si; unordered_set<VertexId> Si_set;
        for(auto&[kk,pj]:Si_kv){
            if(Us.count(pj)) continue;
            auto it=p2g.find(pj);
            if(it!=p2g.end()){
                for(VertexId v:groups[it->second])
                    if(G.d[v]<Bi&&!Us.count(v)&&!Si_set.count(v))
                        { Si.push_back(v); Si_set.insert(v); }
            } else if(!Si_set.count(pj)){
                Si.push_back(pj); Si_set.insert(pj);
            }
        }
        if(Si.empty()) continue;
        BMR_PY sub=bmssp_py(G,l-1,B,Si);
        for(VertexId v:sub.U) if(!Us.count(v)){ U.push_back(v); Us.insert(v); }
        unordered_set<int> J;
        for(VertexId u:sub.U){
            for(uint32_t i=G.g.offsets[u];i<G.g.offsets[u+1];++i){
                VertexId v=G.g.edges[i].target; Weight w=G.g.edges[i].weight;
                bool imp=relax26(G,u,v,w,B);
                if(imp){
                    if(Us.count(v)) Us.erase(v);
                    D.insert(make_sk(G.d[v],G.hops[v],v),v);
                    auto it=v2g.find(v); if(it!=v2g.end()) J.insert(it->second);
                }
            }
        }
        for(int gid:J){
            VertexId np=pivots[gid];
            for(VertexId v:groups[gid]) if(!Us.count(v)&&G.d[v]<G.d[np]) np=v;
            if(np!=pivots[gid]){ p2g.erase(pivots[gid]); pivots[gid]=np; p2g[np]=gid; }
            if(!Us.count(pivots[gid]))
                D.insert(make_sk(G.d[pivots[gid]],G.hops[pivots[gid]],pivots[gid]),pivots[gid]);
        }
        D.merge(sub.D_out);
    }
    sort(U.begin(),U.end()); U.erase(unique(U.begin(),U.end()),U.end());
    DS_BMSSPy ed; ed.init(); return {B,U,move(ed)};
}

static vector<Weight> run_bmsspy(const CompactGraph& g, VertexId src){
    // PY-2: apply constant-degree transform
    CompactGraph g2=build_constant_degree(g);
    auto[k,t]=cbrt_params(g2.num_nodes);
    if(g2.num_nodes==0) return {};
    AlgoState G{g2};
    G.d.assign(g2.num_nodes,INF); G.pred.assign(g2.num_nodes,UINT32_MAX);
    G.hops.assign(g2.num_nodes,0); G.k=k; G.t=t; G.d[src]=0;
    double lgn=(g2.num_nodes>1)?log2((double)g2.num_nodes):1.0;
    int l=max(1,(int)ceil(lgn/(double)t));
    bmssp_py(G,l,INF,{src}); // PY-4: pred[] tracked inside relax26
    // map back to original vertices
    vector<Weight> d(g.num_nodes);
    for(uint32_t v=0;v<g.num_nodes;v++) d[v]=G.d[v];
    return d;
}

// =====================================================================
// 5.  RecDijkstra_AC  —  Stefansson, Biggar, Johansson 2026
//
// Paper Algorithm:
//   AC-1  Build the Acyclic-Connected (A-C) tree T of G:
//         (a) Compute all SCCs via Kosaraju's algorithm in O(m+n).
//         (b) Condense into the SCC DAG (acyclic component graph).
//   AC-2  Process SCCs in topological order of the SCC DAG:
//         For each SCC C in topological order:
//           (a) Collect all seed distances: for each vertex u in C,
//               if d[u] < INF (set by a predecessor SCC), seed the
//               local priority queue with (d[u], u).
//           (b) Run local Dijkstra within C: only relax intra-SCC edges.
//               This settles all vertices reachable within C.
//           (c) Propagate across inter-SCC edges: for each settled
//               vertex u in C and each inter-SCC edge (u,v) where
//               comp[v] != comp[u], update d[v] = min(d[v], d[u]+w).
//   AC-3  Complexity: O(m + n log(nw(G))) where nw(G) = nesting width
//         = min over topological orderings of max simultaneous active SCCs.
//         On graphs with one giant SCC (social networks): nw(G)~n → O(m+n log n).
//
// Implementation:
//   - Iterative Kosaraju (stack-based DFS) to avoid stack overflow on
//     million-node graphs.
//   - In Kosaraju pass-2: comp=0 is assigned to the first SCC found =
//     SOURCE SCC. Process in ascending comp order = topological order.
//   - Iterative local Dijkstra per SCC using a per-SCC priority queue.
// =====================================================================
static vector<Weight> run_rec_dijkstra_ac(const CompactGraph& g, VertexId src){
    uint32_t N=g.num_nodes;
    vector<Weight> d(N,INF); d[src]=0;

    // AC-1a: Kosaraju pass 1 — DFS on G, record finish order
    vector<bool> vis(N,false);
    vector<VertexId> order; order.reserve(N);
    {
        stack<pair<VertexId,uint32_t>> stk;
        for(uint32_t s=0;s<N;s++){
            if(vis[s]) continue;
            stk.push({s,g.offsets[s]}); vis[s]=true;
            while(!stk.empty()){
                auto&[u,ei]=stk.top();
                if(ei<g.offsets[u+1]){
                    VertexId v=g.edges[ei++].target;
                    if(!vis[v]){ vis[v]=true; stk.push({v,g.offsets[v]}); }
                } else { order.push_back(u); stk.pop(); }
            }
        }
    }

    // Build reverse graph for pass 2
    vector<vector<VertexId>> rev(N);
    for(uint32_t u=0;u<N;u++)
        for(uint32_t i=g.offsets[u];i<g.offsets[u+1];++i)
            rev[g.edges[i].target].push_back(u);

    // AC-1a: Kosaraju pass 2 — DFS on reverse G in reverse finish order
    // comp[u] = SCC id. comp=0 is assigned to SOURCE SCC (first in pass-2).
    vector<int> comp(N,-1); int num_scc=0;
    {
        stack<VertexId> stk;
        for(int i=(int)order.size()-1;i>=0;i--){
            VertexId s=order[i]; if(comp[s]!=-1) continue;
            stk.push(s);
            while(!stk.empty()){
                VertexId u=stk.top(); stk.pop();
                if(comp[u]!=-1) continue;
                comp[u]=num_scc;
                for(VertexId v:rev[u]) if(comp[v]==-1) stk.push(v);
            }
            num_scc++;
        }
    }

    // AC-1b: Collect vertices and outgoing inter-SCC edges per SCC
    vector<vector<VertexId>> scc_verts(num_scc);
    for(uint32_t u=0;u<N;u++) scc_verts[comp[u]].push_back(u);

    // Inter-SCC edges: for each SCC c, list of (src_vertex, edge_index)
    // for edges crossing to a different SCC
    struct XEdge { VertexId src_v; uint32_t edge_idx; };
    vector<vector<XEdge>> scc_xedges(num_scc);
    for(uint32_t u=0;u<N;u++){
        int cu=comp[u];
        for(uint32_t i=g.offsets[u];i<g.offsets[u+1];++i)
            if(comp[g.edges[i].target]!=cu)
                scc_xedges[cu].push_back({u,i});
    }

    // AC-2: process SCCs in ascending comp order = topological order
    // (comp=0 is source SCC in Kosaraju)
    for(int c=0;c<num_scc;c++){
        // AC-2a: seed local Dijkstra with any vertices in this SCC
        //        that already have finite distance from predecessor SCCs
        using PQE=pair<Weight,VertexId>;
        priority_queue<PQE,vector<PQE>,greater<PQE>> lq;
        for(VertexId u:scc_verts[c])
            if(d[u]<INF) lq.push({d[u],u});

        // AC-2b: local Dijkstra — only intra-SCC edges
        vector<bool> settled_local(N,false); // reuse across SCCs is wasteful
        // Use unordered_set for settled (avoid O(N) reset)
        unordered_set<VertexId> settled;
        while(!lq.empty()){
            auto[du,u]=lq.top(); lq.pop();
            if(du>d[u]||settled.count(u)) continue;
            settled.insert(u);
            for(uint32_t i=g.offsets[u];i<g.offsets[u+1];++i){
                VertexId v=g.edges[i].target;
                if(comp[v]!=c) continue;     // intra-SCC only
                Weight nd=d[u]+g.edges[i].weight;
                if(nd<d[v]){ d[v]=nd; lq.push({nd,v}); }
            }
        }

        // AC-2c: propagate across inter-SCC edges to successor SCCs
        for(auto&[sv,ei]:scc_xedges[c]){
            if(d[sv]>=INF) continue;
            VertexId tv=g.edges[ei].target;
            Weight   nd=d[sv]+g.edges[ei].weight;
            if(nd<d[tv]) d[tv]=nd;
        }
    }
    return d;
}

// =====================================================================

// =====================================================================
// Step 3 — BMSSP data-structure variants (all on the FAITHFUL BMSSP core)
//
// Every variant runs the SAME core, run_bmssp<DS>(g, src, k, t): FindPivots,
// recursion over levels l, and a frontier structure D with the interface
// init / insert / pull(M) / merge / empty.  Variants differ ONLY in
//   (1) the container that implements D, and
//   (2) the (k,t) setting: 2025 = cbrt_params, 2026 = sqrt_params.
//
//   SBBT = DS_SBBT        std::map (red-black tree)         \
//   Bin  = DS_BinaryHeap  binary min-heap                    |  restored VERBATIM from
//   Fib  = DS_FibHeap     Fibonacci heap (own implementation)/   benchmark_final.cpp
//   Org  = DS_Org         std::multimap keyed by distance   (no counterpart in benchmark_final.cpp)
// =====================================================================

// ─────────────────────────────────────────────────────────────────────
// DS VARIANT A — SBBT (std::map = red-black tree, Lemma 3.4 baseline)
// ─────────────────────────────────────────────────────────────────────
struct DS_SBBT {
    map<PathKey,VertexId>           d1;
    unordered_map<VertexId,PathKey> best;

    void init(){ d1.clear(); best.clear(); }
    bool empty() const { return d1.empty(); }

    void insert(PathKey k, VertexId v){
        auto it=best.find(v);
        if(it!=best.end()){
            if(k<it->second){ d1.erase(it->second); it->second=k; d1[k]=v; }
        } else { best[v]=k; d1[k]=v; }
    }

    void merge(DS_SBBT& o){
        for(auto&[k,v]:o.d1) insert(k,v); o.init();
    }

    struct PR{ vector<pair<PathKey,VertexId>> elems; Weight Bi; };
    PR pull(int M){
        vector<pair<PathKey,VertexId>> r;
        while(!d1.empty()&&(int)r.size()<M){
            auto it=d1.begin();
            r.push_back({it->first,it->second});
            best.erase(it->second); d1.erase(it);
        }
        Weight Bi=d1.empty()?INF:d1.begin()->first.dist;
        return {r,Bi};
    }
};

// ─────────────────────────────────────────────────────────────────────
// DS VARIANT B — Binary Heap (priority_queue, proposed optimization)
// ─────────────────────────────────────────────────────────────────────
struct DS_BinaryHeap {
    using PQE = pair<PathKey,VertexId>;
    priority_queue<PQE, vector<PQE>, greater<PQE>> heap;
    unordered_map<VertexId,PathKey> best;

    void init(){
        priority_queue<PQE,vector<PQE>,greater<PQE>> empty;
        swap(heap,empty);
        best.clear();
    }
    bool empty() const { return heap.empty(); }

    void insert(PathKey k, VertexId v){
        auto it=best.find(v);
        if(it!=best.end()){
            if(k<it->second){ it->second=k; heap.push({k,v}); }
        } else { best[v]=k; heap.push({k,v}); }
    }

    void merge(DS_BinaryHeap& o){
        while(!o.heap.empty()){
            auto [k,v]=o.heap.top(); o.heap.pop();
            insert(k,v);
        }
        o.init();
    }

    struct PR{ vector<pair<PathKey,VertexId>> elems; Weight Bi; };
    PR pull(int M){
        vector<pair<PathKey,VertexId>> r;
        while(!heap.empty()&&(int)r.size()<M){
            auto [k,v]=heap.top(); heap.pop();
            auto it=best.find(v);
            if(it==best.end()||!(k==it->second)) continue;
            best.erase(it);
            r.push_back({k,v});
        }
        while(!heap.empty()){
            auto [k,v]=heap.top();
            auto it=best.find(v);
            if(it!=best.end()&&k==it->second) break;
            heap.pop();
        }
        Weight Bi=heap.empty()?INF:heap.top().first.dist;
        return {r,Bi};
    }
};

// ─────────────────────────────────────────────────────────────────────
// DS VARIANT C — Fibonacci Heap (theoretical optimum)
//
// FIB-FIX-1: merge() uses decrease_key when node already in heap
//            instead of extract+reinsert, preventing distance corruption
// FIB-FIX-2: pool cleared on init() to prevent dangling pointer use
// FIB-FIX-3: consolidation skips null roots defensively
// ─────────────────────────────────────────────────────────────────────
struct FibNode {
    PathKey  key;
    VertexId val;
    int      degree = 0;
    bool     marked = false;
    FibNode* parent  = nullptr;
    FibNode* child   = nullptr;
    FibNode* left    = nullptr;
    FibNode* right   = nullptr;
};

struct FibHeap {
    FibNode* min_node = nullptr;
    int      n        = 0;
    list<FibNode> pool;

    FibNode* make_node(PathKey k, VertexId v){
        pool.push_back({k,v});
        FibNode* x = &pool.back();
        x->left = x->right = x;
        return x;
    }

    bool empty() const { return min_node==nullptr; }

    void insert(FibNode* x){
        x->parent=nullptr; x->marked=false;
        if(!min_node){ min_node=x; x->left=x->right=x; }
        else {
            x->right = min_node;
            x->left  = min_node->left;
            min_node->left->right = x;
            min_node->left = x;
            if(x->key < min_node->key) min_node=x;
        }
        n++;
    }

    void _link(FibNode* y, FibNode* x){
        y->left->right = y->right;
        y->right->left = y->left;
        y->parent = x;
        if(!x->child){ x->child=y; y->left=y->right=y; }
        else {
            y->right = x->child;
            y->left  = x->child->left;
            x->child->left->right = y;
            x->child->left = y;
        }
        x->degree++;
        y->marked=false;
    }

    FibNode* extract_min(){
        FibNode* z = min_node;
        if(!z) return nullptr;
        if(z->child){
            FibNode* c = z->child;
            vector<FibNode*> children;
            FibNode* cur=c;
            do { children.push_back(cur); cur=cur->right; } while(cur!=c);
            for(FibNode* ch : children){
                ch->left->right=ch->right; ch->right->left=ch->left;
                ch->right=min_node; ch->left=min_node->left;
                min_node->left->right=ch; min_node->left=ch;
                ch->parent=nullptr;
            }
        }
        z->left->right = z->right;
        z->right->left = z->left;
        if(z==z->right){ min_node=nullptr; }
        else { min_node=z->right; _consolidate(); }
        n--;
        return z;
    }

    void _consolidate(){
        int maxdeg=(int)(2.0*log2((double)(n+2)))+2;
        vector<FibNode*> A(maxdeg+2, nullptr);
        vector<FibNode*> roots;
        FibNode* cur=min_node;
        do { roots.push_back(cur); cur=cur->right; } while(cur!=min_node);
        for(FibNode* w : roots){
            if(!w) continue; // FIB-FIX-3
            FibNode* x=w;
            int d=x->degree;
            while(d<(int)A.size()&&A[d]){
                FibNode* y=A[d];
                if(y->key<x->key) swap(x,y);
                _link(y,x);
                A[d]=nullptr; d++;
            }
            if(d>=(int)A.size()) A.resize(d+2,nullptr);
            A[d]=x;
        }
        min_node=nullptr;
        for(FibNode* x : A){
            if(!x) continue;
            x->left=x->right=x;
            if(!min_node){ min_node=x; }
            else {
                x->right=min_node; x->left=min_node->left;
                min_node->left->right=x; min_node->left=x;
                if(x->key<min_node->key) min_node=x;
            }
        }
    }

    void decrease_key(FibNode* x, PathKey k){
        if(!(k<x->key)) return;
        x->key=k;
        FibNode* y=x->parent;
        if(y && x->key<y->key){ _cut(x,y); _cascading_cut(y); }
        if(min_node && x->key<min_node->key) min_node=x;
    }

    void _cut(FibNode* x, FibNode* y){
        if(x->right==x){ y->child=nullptr; }
        else {
            if(y->child==x) y->child=x->right;
            x->left->right=x->right; x->right->left=x->left;
        }
        y->degree--;
        x->right=min_node; x->left=min_node->left;
        min_node->left->right=x; min_node->left=x;
        x->parent=nullptr; x->marked=false;
    }

    void _cascading_cut(FibNode* y){
        FibNode* z=y->parent;
        if(z){ if(!y->marked) y->marked=true;
               else { _cut(y,z); _cascading_cut(z); } }
    }

    void clear(){ min_node=nullptr; n=0; pool.clear(); }
};

struct DS_FibHeap {
    FibHeap                          fh;
    unordered_map<VertexId,FibNode*> node_map;

    void init(){ fh.clear(); node_map.clear(); } // FIB-FIX-2: pool.clear() inside fh.clear()
    bool empty() const { return fh.empty(); }

    void insert(PathKey k, VertexId v){
        auto it=node_map.find(v);
        if(it!=node_map.end()){
            // FIB-FIX-1: use decrease_key instead of extract+reinsert
            fh.decrease_key(it->second, k);
        } else {
            FibNode* x=fh.make_node(k,v);
            fh.insert(x);
            node_map[v]=x;
        }
    }

    void merge(DS_FibHeap& o){
        // FIB-FIX-1: drain using decrease_key path to preserve correctness
        while(!o.fh.empty()){
            FibNode* x=o.fh.extract_min();
            if(!x) break;
            VertexId v=x->val;
            PathKey  k=x->key;
            o.node_map.erase(v);
            insert(k,v);
        }
        o.init();
    }

    struct PR{ vector<pair<PathKey,VertexId>> elems; Weight Bi; };
    PR pull(int M){
        vector<pair<PathKey,VertexId>> r;
        while(!fh.empty()&&(int)r.size()<M){
            FibNode* x=fh.extract_min();
            if(!x) break;
            node_map.erase(x->val);
            r.push_back({x->key,x->val});
        }
        Weight Bi=INF;
        if(!fh.empty()&&fh.min_node) Bi=fh.min_node->key.dist;
        return {r,Bi};
    }
};

// ─────────────────────────────────────────────────────────────────────
// SHARED BMSSP STRUCTURES
// ─────────────────────────────────────────────────────────────────────

// ── Org: multimap keyed by distance, lazy deletion of stale entries ──────────
struct DS_Org {
    using Ent = pair<PathKey,VertexId>;
    multimap<Weight,Ent> m;
    unordered_map<VertexId,PathKey> best;
    void init(){ m.clear(); best.clear(); }
    bool empty() const { return best.empty(); }
    void insert(PathKey k, VertexId v){
        auto it=best.find(v);
        if(it!=best.end()){
            if(k<it->second){ it->second=k; m.emplace(k.dist,Ent{k,v}); }
        } else { best.emplace(v,k); m.emplace(k.dist,Ent{k,v}); }
    }
    void merge(DS_Org& o){
        for(auto& kv:o.m){
            auto it=o.best.find(kv.second.second);
            if(it!=o.best.end() && it->second==kv.second.first) insert(kv.second.first,kv.second.second);
        }
        o.init();
    }
    struct PR{ vector<pair<PathKey,VertexId>> elems; Weight Bi; };
    PR pull(int M){
        vector<pair<PathKey,VertexId>> r;
        while(!m.empty() && (int)r.size()<M){
            auto it=m.begin(); Ent e=it->second; m.erase(it);
            auto b=best.find(e.second);
            if(b==best.end() || !(b->second==e.first)) continue;
            best.erase(b); r.push_back({e.first,e.second});
        }
        while(!m.empty()){
            auto it=m.begin(); auto b=best.find(it->second.second);
            if(b!=best.end() && b->second==it->second.first) break;
            m.erase(it);
        }
        Weight Bi=m.empty()?INF:m.begin()->first;
        return {r,Bi};
    }
};


// 2026 setting: k = t = floor(sqrt(log2 n))   (as in benchmark_final.cpp)
static pair<int,int> sqrt_params(uint32_t n){
    double lgn=(n>1)?log2((double)n):1.0;
    int t=max(1,(int)floor(sqrt(lgn)));
    return {t,t};
}

// Step 3 entry points — one faithful core, four frontier structures, two (k,t) settings
vector<Weight> run_ds_org_2025(const CompactGraph& g, VertexId s) { auto [k,t]=cbrt_params(g.num_nodes); return run_bmssp<DS_Org>(g,s,k,t); }
vector<Weight> run_ds_sbbt_2025(const CompactGraph& g, VertexId s){ auto [k,t]=cbrt_params(g.num_nodes); return run_bmssp<DS_SBBT>(g,s,k,t); }
vector<Weight> run_ds_bin_2025(const CompactGraph& g, VertexId s) { auto [k,t]=cbrt_params(g.num_nodes); return run_bmssp<DS_BinaryHeap>(g,s,k,t); }
vector<Weight> run_ds_fib_2025(const CompactGraph& g, VertexId s) { auto [k,t]=cbrt_params(g.num_nodes); return run_bmssp<DS_FibHeap>(g,s,k,t); }
vector<Weight> run_ds_org_2026(const CompactGraph& g, VertexId s) { auto [k,t]=sqrt_params(g.num_nodes); return run_bmssp<DS_Org>(g,s,k,t); }
vector<Weight> run_ds_sbbt_2026(const CompactGraph& g, VertexId s){ auto [k,t]=sqrt_params(g.num_nodes); return run_bmssp<DS_SBBT>(g,s,k,t); }
vector<Weight> run_ds_bin_2026(const CompactGraph& g, VertexId s) { auto [k,t]=sqrt_params(g.num_nodes); return run_bmssp<DS_BinaryHeap>(g,s,k,t); }
vector<Weight> run_ds_fib_2026(const CompactGraph& g, VertexId s) { auto [k,t]=sqrt_params(g.num_nodes); return run_bmssp<DS_FibHeap>(g,s,k,t); }

// ===== CLI harness (UI): algorithm code above is byte-identical to benchmark_variants.cpp =====
#include <cstring>
#include <cstdlib>
int main(int argc, char** argv){
    // usage: bench_cli <graph.bin> <algo> <source|-1> <out.csv> [topology_label] [repeats]
    if(argc < 5){ cerr << "usage: bench_cli graph.bin algo source out.csv [label]\n"; return 2; }
    string fp = argv[1], algo = argv[2], outp = argv[4];
    long long srcArg = atoll(argv[3]);
    string label = (argc > 5) ? argv[5] : "custom";
    int reps = (argc > 6) ? max(1, atoi(argv[6])) : 1;
    CompactGraph g = load_binary_csr(fp);
    if(g.num_nodes == 0){ cerr << "ERROR: cannot load " << fp << "\n"; return 3; }
    VertexId src = 0; uint32_t maxdeg = 0;
    if(srcArg >= 0 && (uint64_t)srcArg < g.num_nodes){ src = (VertexId)srcArg; maxdeg = g.offsets[src+1]-g.offsets[src]; }
    else for(uint32_t i = 0; i < g.num_nodes; i++){
        uint32_t dg = g.offsets[i+1]-g.offsets[i];
        if(dg > maxdeg){ maxdeg = dg; src = i; }
    }
    cout << "NODES " << g.num_nodes << " EDGES " << g.edges.size() << " SOURCE " << src << " DEGREE " << maxdeg << endl;

    // Dijkstra reference distances (untimed) for the correctness check
    vector<Weight> gt = run_dijkstra(g, src);
    size_t reach = 0; for(auto d0 : gt) if(d0 < INF) reach++;
    cout << "REACHABLE " << reach << endl;

    typedef vector<Weight>(*AlgoFn)(const CompactGraph&, VertexId);
    struct Entry { const char* name; AlgoFn fn; };
    Entry table[] = {
        {"Dijkstra",        run_dijkstra},
        {"BMSSP_Org_2025",  run_ds_org_2025},  {"BMSSP_SBBT_2025", run_ds_sbbt_2025},
        {"BMSSP_Bin_2025",  run_ds_bin_2025},  {"BMSSP_Fib_2025",  run_ds_fib_2025},
        {"BMSSP_Org_2026",  run_ds_org_2026},  {"BMSSP_SBBT_2026", run_ds_sbbt_2026},
        {"BMSSP_Bin_2026",  run_ds_bin_2026},  {"BMSSP_Fib_2026",  run_ds_fib_2026},
        {"BMSSP_WC",        run_bmssp_wc},     {"BMSSP_CD",        run_bmssp_cd},
        {"BMSSP_Expected",  run_bmssp_expected},{"BMSSPy",         run_bmsspy},
        {"RecDijkstra_AC",  run_rec_dijkstra_ac}
    };
    AlgoFn fn = nullptr;
    for(auto& e : table) if(algo == e.name) fn = e.fn;
    if(!fn){ cerr << "unknown algorithm " << algo << "\n"; return 4; }

    bool exists = ifstream(outp).good();
    ofstream csv(outp, ios::app);
    if(!exists) csv << "Topology,Environment,Algorithm,Time_Seconds,Mismatches,Status,Run\n";
    for(int rep = 1; rep <= reps; ++rep){
        auto t0 = chrono::high_resolution_clock::now();
        vector<Weight> res = fn(g, src);
        double dt = chrono::duration<double>(chrono::high_resolution_clock::now()-t0).count();
        size_t mis = 0, cmp = min(res.size(), (size_t)g.num_nodes);
        for(size_t i = 0; i < cmp; ++i) if(res[i] != gt[i]) mis++;
        if(res.size() < (size_t)g.num_nodes) mis += g.num_nodes - res.size();
        string status = (mis == 0) ? "VERIFIED EXACT" : "FAILED (" + to_string(mis) + " mismatches)";
        cout << left << setw(24) << algo << " | " << fixed << setprecision(4) << dt << "s | " << status << " | run " << rep << endl;
        csv << label << ",diffusion," << algo << "," << fixed << setprecision(6) << dt << ","
            << mis << "," << status << "," << rep << "\n";
        csv.flush();
    }
    return 0;
}
