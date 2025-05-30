/**
 *
 * Jose Lugo-Martinez, jlugomar@indiana.edu
 * School of Informatics and Computing
 * Indiana University-Bloomington 
 *
 * Aug-6-2012
 *
 * Copyright (c) 2014 Jose Lugo-Martinez,
 * Vladimir Vacic, and Predrag Radivojac.
 *
 */

#ifndef UTILS_H
#define	UTILS_H

#include "config.h"
#include <vector>
using namespace std;


float compare_labels(char label1, char label2);

int randint( int max); //[0.00-1.00)

double randdouble();

unsigned long get_graphlet_length(unsigned long g_type);

void compare_two(char &a, char &b, char &a1, char &b1);

unsigned int set_k(unsigned long g_type, float sf);

void insert_permutation(const Key &target, vector<Key> &mismatches_list);

string get_key(Key k);

string print_key(Key k);

Key make_key(char root, char a, char b, unsigned long g_type);

void initialize_vertices_labels(Key key, char &root, char &a, char &b);

Key get_feature_id(Key k, unsigned long g_type);

void increment_match_hash(map<Key,MismatchInfo> &hash, const Key &k, vector<Key> &mismatches_list);

float retrieve_exact_matches_count(map<Key,MismatchInfo> &hash, const Key &k);

float retrieve_edge_mismatch_count(map<Key,MismatchInfo> &hash, const Key &k);

float retrieve_label_mismatch_count(unsigned long g_type, map<Key,MismatchInfo> &hash, const Key key);

void insert_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, const Key &k);

void insert_mismatches_hash(map<Key,MismatchInfo> &mismatch_hash, const Key &k, vector<Key> mismatches_list);

Key create_permutations_subset(vector<Key> &mismatches, char root, char a, char b, unsigned long g_type);

void insert_graphlet_mismatch_neighborhood(list<Key> &neighborhood, Key k);

void generate_graphlet_mismatch_neighborhood_m1(list<Key> &neighborhood, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_vlm_matrix, unsigned long g_type, Key key);

void generate_graphlet_mismatch_neighborhood_m2(list<Key> &neighborhood, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_vlm_matrix, unsigned long g_type, Key key);

void generate_vertex_label_mismatch_graphlets(map<Key, list<Key> > &vl_mismatch_neighborhood, map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, Key key, unsigned long g_type, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_vlm_matrix, int VLM);

void increment_mismatch_count(map<Key,MismatchInfo> &hash, const Key &key, const Key &mismatch_key, float sim_score, float mult_factor);

float compare_graphlets(Key key1, Key key2, unsigned long g_type, map<string, float> &sim_vlm_matrix, float &sim_score);

void update_mismatch_count(map<Key,MismatchInfo> &hash, Key k, float mult_factor, unsigned long g_type, map<string, float> sim_vlm_matrix, int VLM, bool eq) ;

void increment_edge_mismatch_hash(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, const Key &k, float mult_factor, vector<Key> &mismatches_list);

void insert_edge_mismatch_graphlet(list<pair <unsigned long, Key> > &EM_set, pair <unsigned long, Key> graphlet, unsigned index);

void update_edge_mismatch_count(vector<list<pair <unsigned long, Key> > > &EM_set, Key k, unsigned g_type, unsigned EDGE_MISMATCHES_ALLOWED, unsigned vindex);

#endif

