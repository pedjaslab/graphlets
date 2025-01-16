/**
 *
 * Jose Lugo-Martinez, jlugomar@indiana.edu
 * School of Informatics and Computing
 * Indiana University-Bloomington 
 *
 * May-09-2013
 *
 * Copyright (c) 2014 Jose Lugo-Martinez,
 * Vladimir Vacic, and Predrag Radivojac.
 *
 */

#ifndef __DIGKERNEL_H__
#define __DIGKERNEL_H__

#include "utils.h"
#include "simgraph.h"
#include <fstream>
#include <utility>
#include <list>
#include <map>
#include <vector>
using namespace std;

class DigraphKernel  {
public:
    DigraphKernel() : NORMALIZE(false), VERBOSE(false), SF(0.0), EM(0)  {}
    ~DigraphKernel()  {}

    /** Read a directed graph, node labels, and list of vertices of interest over input directed graph. */
    void read_graphs(string, string, const vector<unsigned> &);

    /** Read a probability similarity matrix for each vertex label as means to weight each label substitution. */
    void read_sim_matrix(string filename);    

    /** */
    inline void set_labels(const vector<int> &l)  { labels = l; }

    /** Compute cumulative random walk kernel matrix. */
    void compute_random_walk_cumulative_matrix(int steps, double restart);

    /** Compute random walk kernel matrix. */
    void compute_random_walk_matrix(int steps, double restart);

    /** Compute label substitutions kernel matrix. */
    void compute_label_mismatch_matrix();

    /** Compute edge indels kernel matrix. */
    void compute_edge_mismatch_matrix();

    /** Compute edit distance kernel matrix with 1 operation. */
    void compute_edit_distance_matrix();

    /** Compute edit distance kernel matrix with 2 operations. */
    void compute_edit_distance2_matrix();

    /** Writes distance kernel matrix in either:
    1:triangular binary form,
    2:triangular standard output form,
    3: squared-matrix standard output form.
    */
    void write_matrix(const char*);

    /** Writes vector of counts for label substitutions kernel on SVML^light format. */
    void write_sparse_svml_lm(const char*);

    /** Writes vector of counts for edge indels kernel on SVML^light format. */
    void write_sparse_svml_em(const char*);

    /** Writes vector of counts for edit distance kernel (1-operation) on SVML^light format. */
    void write_sparse_svml_ed(const char*);

    /** Writes vector of counts for edit distance kernel (2-operations) on SVML^light format. */
    void write_sparse_svml_ed2(const char*);

    /** Writes class labels for each example. */
    void write_labels(const char*);

    inline void set_normalize()  { NORMALIZE = true; } 

    inline void set_verbose()  { VERBOSE = true; }

    inline void set_number_label_mismatches(float fraction)  { SF = fraction; }
    
    inline void set_label_mismatches_alphabet(string alphabet)  { ALPHABET = alphabet; }

    inline void set_label_mismatches_root_alphabet(string alphabet)  { ALPHABET_ROOT = alphabet; }

	inline void set_number_edges_mismatches(unsigned edges_mismatches)  { EM = edges_mismatches; }
    
private:
	/** Returns the cumulative random walk kernel between two rooted neighborhoods. */
    float random_walk_cumulative(SimpleGraph &g, unsigned g1_root, unsigned g2_root, int steps, double restart);

	/** Returns the random walk kernel between two rooted neighborhoods. */
    float random_walk(SimpleGraph &g, unsigned g1_root, unsigned g2_root, int steps, double restart);

    /** Returns the counts of nonisomorphic labeled digraphlets on a rooted neighborhood. */
    vector<map<Key,MismatchInfo> > get_graphlets_counts(SimpleGraph &g, unsigned g_root);

    /** Adds the counts for inexact digraphlets based on vertex label substitutions. */
    void add_vertex_label_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, unsigned long g_type, int VLM, bool option);

    /** Updates the vector of count with inexact digraphlets based on vertex and edge label substitutions. */
    void update_label_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, unsigned long g_type, bool option, int VLM, bool eq);

    /** Adds the counts for inexact digraphlets based on edge insertions and deletions. */
    void add_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash);

    /** Adds the counts for inexact digraphlets based on 1-edge insertion and deletion. */
    void add_1_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash);

    /** Adds the counts for inexact digraphlets based on 2-edge insertions and deletions. */
    void add_2_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash);

    /** Normalizes the kernel matrix using the method for normalizing the spectral kernel matrix. */
    void normalize_spectral(map<Key,MismatchInfo> &, unsigned long g_type);

    /** Computes graphlet distance between two vector of digraphlet counts.*/
    float distance_hash_join(map<Key,MismatchInfo>, map<Key,MismatchInfo>, unsigned long g_type);

    /** Normalizes the kernel matrix using the method for normalizing the spectral kernel matrix. */
    void normalize_spectral(vector<map<Key,MismatchInfo> >&);

    /** Computes graphlet distance between two vector of digraphlet counts. */
    float distance_hash_join(vector<map<Key,MismatchInfo> >, vector<map<Key,MismatchInfo> >);

    // Data members.
    bool NORMALIZE, VERBOSE;
    float SF;
	unsigned EM;
    string ALPHABET;
    string ALPHABET_ROOT;

    vector<int> labels;
    SimpleGraph digraph;
    vector<unsigned>    roots;       // Vertices of interest.
    map<string,float>   sim_vlm_matrix;
    vector<vector<map<Key,MismatchInfo> > > hashes;
    vector<vector<float> >  kernel;
    map<Key, list<Key> >    vl_mismatch_neighborhood;
};

#endif

