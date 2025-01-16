#include "utils.h"
#include "string.h"
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <queue>
#include <iomanip>
using namespace std;


/************************ Auxiliary functions for Digraph Kernel class *************************/
float compare_labels(char label1, char label2)  {
    if (label1 == label2)
        return 1.0;
    return 0.0;
}
int randint(int max)  {
    if( max == 0)   {
		cerr << "ERROR: Check function randint()" << endl;
        return -1;
    }
    return rand() % max;
}

double randdouble()  {
    return randint(101) / 100.0;
}

unsigned long get_graphlet_length(unsigned long g_type)  {
    unsigned long g_length;
    if (g_type == 0)  {
        g_length = 1;
    }
    else if (g_type >= 1 && g_type <= 3)  {
        g_length = 2;
    }
    else if (g_type >= 4 && g_type <= 33)  {
        g_length = 3; 
    }
    else  {
        cerr << "ERROR: Digraphlet type: " << g_type << " is unsupported." << endl; exit(1);
    }
    return g_length;
}

void compare_two(char &a, char &b, char &a1, char &b1)  {
	if (a < b)  { // a,b
		a1 = a; b1 = b;
	} 
	else  { // b,a
		a1 = b; b1 = a;
	}
}

unsigned int set_k(unsigned long g_type, float sf)  {
    unsigned int k_mismatches = 0;   
    if ((g_type >= 0 && g_type < DIGRAPHLETS_TYPES) && (sf >= 0.0 && sf <= 1.0))  {
        k_mismatches = unsigned(float(get_graphlet_length(g_type)) * sf);
    }
    else  {
        if (g_type < 0 || g_type >= DIGRAPHLETS_TYPES)  {
            cerr << "ERROR: Digraphlet type: " << g_type << " is unsupported." << endl; exit(1);
        }
        else  {
            cerr << "ERROR: Scaling fraction: " << sf << " is invalid." << endl; exit(1);
        }
    }
    return k_mismatches;
}
   
/* This function first checks if graphlet permutation is already listed.
 * If not, it inserts it into a set of allowed permutations.
 */  
void insert_permutation(const Key &target, vector<Key> &mismatches_list)  {
	bool found(false);
	for (vector<Key>::iterator mismatch_key = mismatches_list.begin(); mismatch_key < mismatches_list.end(); mismatch_key++)  {
		if (*mismatch_key == target)  {
			found = true;
			break;
		}
	}
	if (!found)  {
	    mismatches_list.push_back(target);
	}
}

string get_key(Key k)  {	
	string temp(DIGRAPHLET_SIZE,ZERO_CHAR);
    for (unsigned i=DIGRAPHLET_SIZE-1; i > 0; i--)  {
        temp[i] = (k & ALPHABET_SIZE) + ZERO_CHAR; k = k >> LOG_ALPHABET_SIZE;
    }
    temp[0] = (k & ALPHABET_SIZE) + ZERO_CHAR; //root  
	
    return temp;
}

string print_key(Key k)  {
	string temp(DIGRAPHLET_SIZE,ZERO_CHAR);
    for (unsigned i=DIGRAPHLET_SIZE-1; i > 0; i--)  {
        temp[i] = (k & ALPHABET_SIZE) + ZERO_CHAR; k = k >> LOG_ALPHABET_SIZE;
    }
    temp[0] = (k & ALPHABET_SIZE) + ZERO_CHAR; //root  
	
    ostringstream s;
    s << temp ;
    return s.str();
}

Key make_key(char root, char a, char b, unsigned long g_type)  {
    Key curr_key(0);

	curr_key = ((unsigned long) (root-ZERO_CHAR)) << LOG_ALPHABET_SIZE;
    curr_key = (curr_key+a-ZERO_CHAR) << LOG_ALPHABET_SIZE;
	curr_key = (curr_key+b-ZERO_CHAR);
	
    return curr_key;
}

void initialize_vertices_labels(Key key, char &root, char &a, char &b)  {
    b = (key & ALPHABET_SIZE) + ZERO_CHAR;  key = key >> LOG_ALPHABET_SIZE;
    a = (key & ALPHABET_SIZE) + ZERO_CHAR;  key = key >> LOG_ALPHABET_SIZE;
    root = (key & ALPHABET_SIZE) + ZERO_CHAR;
}

Key get_feature_id(Key k, unsigned long g_type)  {
    Key feature_id(0);
    feature_id = (k << LOG_DIGRAPHLET_TYPES_SIZE) + g_type;
    return feature_id;
}

void increment_match_hash(map<Key,MismatchInfo> &hash, const Key &k, vector<Key> &mismatches_list)  {
    map<Key,MismatchInfo>::iterator it;
    if ((it = hash.find(k)) == hash.end())  {
        hash[k].matches = 1.0;
        hash[k].mismatches = 0.0;
        for (vector<Key>::iterator mismatch_key = mismatches_list.begin(); mismatch_key < mismatches_list.end(); mismatch_key++)  {
            hash[k].mismatchesGraph[*mismatch_key] = 0.0;
        }
    }
    else  {
        hash[k].matches = it->second.matches + 1.0;
    }
	mismatches_list.clear();
}

float retrieve_exact_matches_count(map<Key,MismatchInfo> &hash, const Key &k)  {
    float perfect_matches = 0.0;
    map<Key,MismatchInfo>::iterator it = hash.find(k);
    if (it != hash.end())  {
        perfect_matches = hash[k].matches ;
    }
    return perfect_matches;
}

float retrieve_edge_mismatch_count(map<Key,MismatchInfo> &hash, const Key &k)  {
    float counts = 0.0;
    map<Key,MismatchInfo>::iterator it = hash.find(k);
    if (it != hash.end())  {
        counts = hash[k].matches + hash[k].mismatches ;
    }
    return counts;
}

float retrieve_label_mismatch_count(unsigned long g_type, map<Key,MismatchInfo> &hash, const Key key)  {
    map<Key,MismatchInfo>::iterator it = hash.find(key);
    float counts = 0.0;
    if (it != hash.end())  {
	    counts = hash[key].matches + hash[key].mismatches;
	    for (map<Key,float>::iterator mismatches = hash[key].mismatchesGraph.begin(); mismatches != hash[key].mismatchesGraph.end(); mismatches++)  {
		    counts += mismatches->second;
	    }
    }
	return counts;
}

void insert_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, const Key &k)  {
    map<Key,MismatchInfo>::iterator it;    
    if ((it = hash.find(k)) == hash.end())  {
        hash[k].matches = mismatch_hash[k].matches;
        hash[k].mismatches = mismatch_hash[k].mismatches;
        for (map<Key,float>::iterator mismatch_key = mismatch_hash[k].mismatchesGraph.begin(); mismatch_key != mismatch_hash[k].mismatchesGraph.end(); mismatch_key++)  {
            hash[k].mismatchesGraph[mismatch_key->first] = mismatch_hash[k].mismatchesGraph[mismatch_key->first];
        }
    }
    else  {
        cerr << "Warning: While merging found a double counted digraphlet " << print_key(k) << endl;
    }
}

void insert_mismatches_hash(map<Key,MismatchInfo> &mismatch_hash, const Key &k, vector<Key> mismatches_list)  {
    map<Key,MismatchInfo>::iterator it;
	mismatch_hash[k].matches = 0.0;
    mismatch_hash[k].mismatches = 0.0;
	for (vector<Key>::iterator mismatch_key = mismatches_list.begin(); mismatch_key < mismatches_list.end(); mismatch_key++)  {
		mismatch_hash[k].mismatchesGraph[*mismatch_key] = 0.0;
    }
}

// Create all vertex-labeled graphlet permutations for each digraphlet found.
Key create_permutations_subset(vector<Key> &mismatches, char root, char a, char b, unsigned long g_type)  {
	Key k(0), curr_key(0); //For label mismatches labels.
	char a1, b1;

	mismatches.clear();

	switch (g_type)  {
		case 0:
			// Make key label
			// R (Root) 
			k = make_key(root, ZERO_CHAR, ZERO_CHAR, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			break;
			
		case 1:
        case 2:
        case 3:
			// Make key label
			// R-A 
			k = make_key(root, a, ZERO_CHAR, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			break;
			
		case 13:
		case 17:
        case 18:
        case 19:
        case 27:
        case 32:
			// Make key label
			// R-A1-A2
			compare_two(a, b, a1, b1);
			k = make_key(root, a1, b1, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// R-A2-A1
			curr_key = make_key(root, b1, a1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;
			
		case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 14:
        case 15:
        case 16:
        case 20:
        case 21:
        case 22:
        case 23:
        case 24:
        case 25:
        case 26:
        case 28:
        case 29:
        case 30:
        case 31:
        case 33:
			// Make key label
			// R-A-B                 
			k = make_key(root, a, b, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			break;
        default:
            cerr << "ERROR: Digraphlet type: " << g_type << " is unsupported." << endl; exit(1);
            break;		
	} //End of switch

	return k;
}

void insert_graphlet_mismatch_neighborhood(list<Key> &neighborhood, Key k)  {
    list<Key>::iterator list_it;
    bool found(false);
    for (list_it = neighborhood.begin(); list_it != neighborhood.end(); list_it++)  {
        if (k == *list_it)  {
            found = true;
            return;
        }
    }
    neighborhood.push_back(k);
}

inline float get_sim_score(char vlabel1, char vlabel2, map<string, float> &sim_vlm_matrix)  {
    string lookup;
    lookup += toupper(vlabel1);
    lookup += toupper(vlabel2);
    map<string,float>::iterator git = sim_vlm_matrix.find(lookup);
    if(git != sim_vlm_matrix.end())  {
        return git->second;
    }
    return 0.0;
}

void generate_graphlet_mismatch_neighborhood_m1(list<Key> &neighborhood, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_vlm_matrix, unsigned long g_type, Key key)  {
    Key k;
    char root, a, b;
	vector<Key> mismatches;
    float sim_score(0.0);

	initialize_vertices_labels(key, root, a, b);
    
    for (unsigned i=0; i<ALPHABET_ROOT.length(); i++)  {
        sim_score = get_sim_score(root, ALPHABET_ROOT[i], sim_vlm_matrix);
        if (ALPHABET_ROOT[i] != root && sim_score >= SIMILARITY_THRESHOLD)  {        
            k = create_permutations_subset(mismatches, ALPHABET_ROOT[i], a, b, g_type);
            insert_graphlet_mismatch_neighborhood(neighborhood, k);
        }
    }

    for (unsigned i=0; i<ALPHABET.length(); i++)  {
        sim_score = get_sim_score(a, ALPHABET[i], sim_vlm_matrix);
        if (ALPHABET[i] != a && sim_score >= SIMILARITY_THRESHOLD)  {    
            k = create_permutations_subset(mismatches, root, ALPHABET[i], b, g_type);
            insert_graphlet_mismatch_neighborhood(neighborhood, k);
        }        
        if (get_graphlet_length(g_type) > 2 && ALPHABET[i] != b)  {
            sim_score = get_sim_score(b, ALPHABET[i], sim_vlm_matrix);
            if (sim_score >= SIMILARITY_THRESHOLD)  {
                k = create_permutations_subset(mismatches, root, a, ALPHABET[i], g_type);
                insert_graphlet_mismatch_neighborhood(neighborhood, k);
            }
        }
    }
}

void generate_graphlet_mismatch_neighborhood_m2(list<Key> &neighborhood, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_vlm_matrix, unsigned long g_type, Key key)  {
    Key k;
    char root, a, b;
	vector<Key> mismatches;
    float sim_score1(0.0), sim_score2(0.0);

	initialize_vertices_labels(key, root, a, b);

	for (unsigned i=0; i<ALPHABET_ROOT.length(); i++)  {
        sim_score1 = get_sim_score(root, ALPHABET_ROOT[i], sim_vlm_matrix);
        if (sim_score1 >= SIMILARITY_THRESHOLD)  {
            for (unsigned j=0; j<ALPHABET.length(); j++)  {
                sim_score2 = get_sim_score(a, ALPHABET[j], sim_vlm_matrix);
                if (sim_score2 >= SIMILARITY_THRESHOLD)  {
                    k = create_permutations_subset(mismatches, ALPHABET_ROOT[i], ALPHABET[j], b, g_type);
                    insert_graphlet_mismatch_neighborhood(neighborhood, k);
                }
                if (get_graphlet_length(g_type) > 2)  {
                    sim_score2 = get_sim_score(b, ALPHABET[j], sim_vlm_matrix);
                    if (sim_score2 >= SIMILARITY_THRESHOLD)  {
                        k = create_permutations_subset(mismatches, ALPHABET_ROOT[i], a, ALPHABET[j], g_type);
                        insert_graphlet_mismatch_neighborhood(neighborhood, k);
                    }
                }
            }
        }
    }

    for (unsigned i=0; i<ALPHABET.length(); i++)  { 
        for (unsigned j=0; j<ALPHABET.length(); j++)  {
            if (get_graphlet_length(g_type) > 2)  {
                sim_score1 = get_sim_score(a, ALPHABET[i], sim_vlm_matrix);
                sim_score2 = get_sim_score(b, ALPHABET[j], sim_vlm_matrix);
                if (sim_score1 >= SIMILARITY_THRESHOLD && sim_score2 >= SIMILARITY_THRESHOLD)  {
                    k = create_permutations_subset(mismatches, root, ALPHABET[i], ALPHABET[j], g_type);
                    insert_graphlet_mismatch_neighborhood(neighborhood, k);
                }
            }
        }
    }
}

// Generate corresponding vertex label mismatch graphlets for each graphlet found.	 
void generate_vertex_label_mismatch_graphlets(map<Key, list<Key> > &vl_mismatch_neighborhood, map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, Key key, unsigned long g_type, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_vlm_matrix, int VLM)  {
    char root, a, b;
	list<Key>::iterator list_it;
	map<Key, list<Key> >::iterator lit;

	if ((lit = vl_mismatch_neighborhood.find(key)) == vl_mismatch_neighborhood.end())  {
		if (VLM == 1)  {
			generate_graphlet_mismatch_neighborhood_m1(vl_mismatch_neighborhood[key], ALPHABET_ROOT, ALPHABET, sim_vlm_matrix, g_type, key);
		}
		if (VLM == 2)  {
			generate_graphlet_mismatch_neighborhood_m2(vl_mismatch_neighborhood[key], ALPHABET_ROOT, ALPHABET, sim_vlm_matrix, g_type, key);
		}
    }

    for (list_it = vl_mismatch_neighborhood[key].begin(); list_it != vl_mismatch_neighborhood[key].end(); list_it++)  {
		map<Key,MismatchInfo>::iterator it;
		map<Key,MismatchInfo>::iterator mit;
		if (((it = hash.find(*list_it)) == hash.end()) && ((mit = mismatch_hash.find(*list_it)) == mismatch_hash.end()))  {
            vector<Key> mismatches;
			initialize_vertices_labels(*list_it, root, a, b);
			Key k = create_permutations_subset(mismatches, root, a, b, g_type);
			//Insert mismatch graphlets into hash
			insert_mismatches_hash(mismatch_hash, k, mismatches);
		}
	}		
}

void increment_mismatch_count(map<Key,MismatchInfo> &hash, const Key &key, const Key &mismatch_key, float sim_score, float mult_factor)  {
	map<Key,MismatchInfo>::iterator it = hash.find(key);
    if (sim_score >= SIMILARITY_THRESHOLD)  {
        if (it != hash.end())  {
		    map<Key,float>::iterator mismatch = hash[it->first].mismatchesGraph.find(mismatch_key);
		    if (mismatch != hash[key].mismatchesGraph.end())  {
			    hash[it->first].mismatchesGraph[mismatch_key] = hash[it->first].mismatchesGraph[mismatch_key] + (mult_factor * sim_score);
		    }
			else  {
				cerr << "ERROR: Mismatch graphlet permutation " << print_key(mismatch_key) << " cannot be found among subset of allowed permutations for mismatch graphlet " << print_key(key) << "." << endl; exit(1);
			}
	    }
		else  {
            cerr << "ERROR: Mismatch derived graphlet " << print_key(key) << " cannot be found." << endl; exit(1);
        }
    }
}

float compare_graphlets(Key key1, Key key2, unsigned long g_type, map<string, float> &sim_vlm_matrix, float &sim_score)  {
	string key1_str = get_key(key1), key2_str = get_key(key2);
	float distance(0.0);
    sim_score = 1.0;

	for (unsigned i=0; i<get_graphlet_length(g_type); i++)  {
        if(key1_str[i] != key2_str[i])  {
            distance = distance + 1.0;
            string lookup;
		    lookup += toupper(key2_str[i]);
		    lookup += toupper(key1_str[i]);
            map<string,float>::iterator git = sim_vlm_matrix.find(lookup);
		    if(git != sim_vlm_matrix.end())  {
                sim_score *= git->second;
		    }
		    else {
		        cerr << "ERROR: Amino acid substiution pair  " << lookup << " cannot be found in similarity matrix." << endl; exit(1);
		    }
        }
	}
	return distance;
}

void update_mismatch_count(map<Key,MismatchInfo> &hash, Key k, float mult_factor, unsigned long g_type, map<string, float> sim_vlm_matrix, int VLM, bool eq)  {
	float curr_dist, min_dist;
    float min_score, sim_score;
	Key min_key, min_graphlet;
	
	if (mult_factor > 0.0 && VLM > 0)  { 
        for (map<Key,MismatchInfo>::iterator it = hash.begin(); it != hash.end(); it++)  {
			if (it->first != k) {
                min_dist = float(get_graphlet_length(g_type)) + 1.0;
                min_score = 0.0;
				min_graphlet = make_key(ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, g_type);
				min_key = make_key(ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, g_type);
				for (map<Key,float>::iterator mismatches_hash = hash[(it->first)].mismatchesGraph.begin(); mismatches_hash != hash[(it->first)].mismatchesGraph.end(); mismatches_hash++)  {
					curr_dist = compare_graphlets(mismatches_hash->first, k, g_type, sim_vlm_matrix, sim_score);
					if (curr_dist < min_dist || (curr_dist == min_dist && sim_score > min_score))  {
						min_dist = curr_dist;
						min_key = it->first;
						min_graphlet = mismatches_hash->first;
						min_score = sim_score;
					}
				}
                if (eq)  {
			        if (min_dist == VLM)  {
				        increment_mismatch_count(hash, min_key, min_graphlet, min_score, mult_factor);
			        }
                }
                else  {
			        if (min_dist <= VLM)  {
				        increment_mismatch_count(hash, min_key, min_graphlet, min_score, mult_factor);
			        }
                }
            }
		}
    }
}

void increment_edge_mismatch_hash(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, const Key &k, float mult_factor, vector<Key> &mismatches_list)  {
    map<Key,MismatchInfo>::iterator it = hash.find(k);
    map<Key,MismatchInfo>::iterator mit = mismatch_hash.find(k);
    if (it != hash.end())  {
        hash[k].mismatches = it->second.mismatches + mult_factor;
    }
    else  {
        if (mit == mismatch_hash.end())  {
            mismatch_hash[k].matches = 0.0;
            mismatch_hash[k].mismatches = mult_factor;
            for (vector<Key>::iterator mismatch_key = mismatches_list.begin(); mismatch_key < mismatches_list.end(); mismatch_key++)  {
                mismatch_hash[k].mismatchesGraph[*mismatch_key] = 0.0;
            }
            mismatches_list.clear();
        }
        else  {
            mismatch_hash[k].mismatches = mit->second.mismatches + mult_factor;
        }
    }
}

void insert_edge_mismatch_graphlet(list<pair <unsigned long, Key> > &EM_set, pair <unsigned long, Key> graphlet, unsigned index)  {
    list<pair <unsigned long, Key> > :: iterator list_it;
    bool found(false);
    for (list_it = EM_set.begin(); list_it != EM_set.end(); list_it++)  {
        if (list_it->first == graphlet.first && list_it->second == graphlet.second)  {
            found = true; 
            break;
        }
    }
    if (!found)
        EM_set.push_back(graphlet);
    return;
}

// Update edge indels digraphlet counter.
void update_edge_mismatch_count(vector<list<pair <unsigned long, Key> > > &EM_set, Key k, unsigned g_type, unsigned EDGE_MISMATCHES_ALLOWED, unsigned vindex)  {
	Key key(k), curr_key(0);
	char root, a, b;
    char a1, b1;
    pair <unsigned long, Key> p (0,curr_key);

    if (EDGE_MISMATCHES_ALLOWED <= 0)
        return;    

	initialize_vertices_labels(key, root, a, b);

	switch (g_type)  {
		case 0: // R 
			// There are no edges addition/removal for this case.
			break;

		case 1: // R<->A 
			// R<-A (Edge Mismatch by removing edge from R to A)			
            p.first = 3;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 3, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A (Edge Mismatch by removing edge from A to R)
            p.first = 2;
            p.second = k;            
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 2, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));			
			break;
			
		case 2: // R->A 
			// R<->A (Edge Mismatch by adding edge from A to R)
            p.first = 1;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);        
            update_edge_mismatch_count(EM_set, k, 1, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 3: // R<-A 
			// R<->A (Edge mismatch by adding edge from R to A)
            p.first = 1;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 1, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 4: // R<->A<->B
			// R<->A<->B and R->B (Edge Mismatch by adding edge from R to B)
            p.first = 22;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A<->B and R<-B (Edge Mismatch by adding edge from B to R)
            p.first = 21;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A<->B (Edge Mismatch by removing edge from R to A)
            p.first = 10;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A<->B (Edge Mismatch by removing edge from A to R)
            p.first = 7;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 7, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A<-B (Edge Mismatch by removing edge from A to B)
            p.first = 6;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 6, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // R<->A->B (Edge Mismatch by removing edge from B to A)
            p.first = 5;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 5, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 5: // R<->A->B
			// R<->A<->B (Edge Mismatch by adding edge B to A)
            p.first = 4;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 4, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A->B and R->B (Edge Mismatch by adding edge from R to B)
            p.first = 25;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A->B and R<-B (Edge Mismatch by adding edge from B to R)
            p.first = 23;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 23, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A->B (Edge Mismatch by removing edge from R to A)
            p.first = 11;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A->B (Edge Mismatch by removing edge from A to R)
            p.first = 8;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 8, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 6: // R<->A<-B
			// R<->A<->B (Edge Mismatch by adding edge A to B)
            p.first = 4;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 4, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A<-B and R->B (Edge Mismatch by adding edge from R to B)
            p.first = 26;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 26, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A<-B and R<-B (Edge Mismatch by adding edge from B to R)
            p.first = 24;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A<-B (Edge Mismatch by removing edge from R to A)
            p.first = 12;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 12, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A<-B (Edge Mismatch by removing edge from A to R)
            p.first = 9;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 9, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 7: // R->A<->B 
			// R<->A<->B (Edge Mismatch by adding edge A to R)
            p.first = 4;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 4, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A<->B and R->B (Edge Mismatch by adding edge from R to B)
            compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 32;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A<->B and R<-B (Edge Mismatch by adding edge from B to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 29;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A<-B (Edge Mismatch by removing edge from A to B)
            p.first = 9;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 9, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A->B (Edge Mismatch by removing edge from B to A)
            p.first = 8;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 8, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 8: // R->A->B
			// R<->A->B (Edge Mismatch by adding edge A to R)
            p.first = 5;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 5, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A->B and R->B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 33;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A->B and R<-B (Edge Mismatch by adding edge from B to R)
            curr_key = make_key(root, a, b, g_type);
            p.first = 30;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A<->B (Edge Mismatch by adding edge B to A)
            p.first = 7;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 7, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 9: // R->A<-B 
			// R<->A<-B (Edge Mismatch by adding edge A to R)
            p.first = 6;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 6, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A<-B and R->B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, b, a, g_type);
            p.first = 33;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A<-B and R<-B (Edge Mismatch by adding edge from B to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 31;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A<->B (Edge Mismatch by adding edge A to B)
            p.first = 7;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 7, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 10: // R<-A<->B 
			// R<->A<->B (Edge Mismatch by adding edge R to A)
            p.first = 4;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 4, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A<->B and R->B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 29;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A<->B and R<-B (Edge Mismatch by adding edge from B to R)
            compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 27;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 27, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A<-B (Edge Mismatch by removing edge from A to B)
            p.first = 12;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 12, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A->B (Edge Mismatch by removing edge from B to A)
            p.first = 11;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, k, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 11: // R<-A->B
			// R<->A->B (Edge Mismatch by adding edge R to A)
            p.first = 5;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 5, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A->B and R->B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 31;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A->B and R<-B (Edge Mismatch by adding edge from B to R)
            curr_key = make_key(root, a, b, g_type);
            p.first = 28;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A<->B (Edge Mismatch by adding edge B to A)
            p.first = 10;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 12: // R<-A<-B
			// R<->A<-B (Edge Mismatch by adding edge R to A)
            p.first = 6;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 6, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A<-B and R->B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, b, a, g_type);
            p.first = 30;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A<-B and R<-B (Edge Mismatch by adding edge from B to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 28;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A<->B (Edge Mismatch by adding edge A to B)
            p.first = 10;
            p.second = k;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, k, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 13: // A1<->R<->A2 
			// A1<->R<->A2 and A1->A2 (Edge Mismatch by adding edge from A1 to A2)
            curr_key = make_key(root, a, b, g_type);
            p.first = 20;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<->R<->A2 and A1<-A2 (Edge Mismatch by adding edge from A2 to A1)
            curr_key = make_key(root, b, a, g_type);
            p.first = 20;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1->R<->A2 (Edge Mismatch by removing edge from R to A1)
            curr_key = make_key(root, b, a, g_type);
            p.first = 15;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 15, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<-R<->A2 (Edge Mismatch by removing edge from A1 to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 14;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<->R<-A2 (Edge Mismatch by removing edge from R to A2)
            curr_key = make_key(root, a, b, g_type);
            p.first = 15;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 15, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // A1<->R->A2 (Edge Mismatch by removing edge from A2 to R)
            curr_key = make_key(root, a, b, g_type);
            p.first = 14;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 14: // A<->R->B 
			// A<->R<->B (Edge Mismatch by adding edge from B to R)
            compare_two(a, b, a1, b1); 
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 13;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 13, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R->B and A->B (Edge Mismatch by adding edge from A to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 25;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R->B and A<-B (Edge Mismatch by adding edge from B to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 26;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 26, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R->B (Edge Mismatch by removing edge from R to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 16;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 16, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // A<-R->B (Edge Mismatch by removing edge from A to R)
            compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 17;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 15: // A<->R<-B
			// A<->R<->B (Edge Mismatch by adding edge from B to R)
            compare_two(a, b, a1, b1); 
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 13;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 13, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R<-B and A->B (Edge Mismatch by adding edge from A to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 23;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 23, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R<-B and A<-B (Edge Mismatch by adding edge from B to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 24;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R<-B (Edge Mismatch by removing edge from R to A)
            compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 18;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 18, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // A<-R<-B (Edge Mismatch by removing edge from A to R)            
            curr_key = make_key(root, b, a, g_type);
            p.first = 16;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 16, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
        
		case 16: // A->R->B
			// A<->R->B (Edge Mismatch by adding edge from R to A)            
            curr_key = make_key(root, a, b, g_type);
            p.first = 14;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R<->B (Edge Mismatch by adding edge from B to R)            
            curr_key = make_key(root, b, a, g_type);
            p.first = 15;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 15, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R->B and A->B (Edge Mismatch by adding edge from A to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 31;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R->B and A<-B (Edge Mismatch by adding edge from B to A)
            curr_key = make_key(root, b, a, g_type);
            p.first = 30;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;                

		case 17: // A1<-R->A2 
			// A1<->R->A2 (Edge Mismatch by adding edge from A1 to R)            
            curr_key = make_key(root, a, b, g_type);
            p.first = 14;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<-R<->A2 (Edge Mismatch by adding edge from A2 to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 14;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<-R->A2 and A1->A2 (Edge Mismatch by adding edge from A1 to A2)
            curr_key = make_key(root, a, b, g_type);
            p.first = 33;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<-R->A2 and A1<-A2 (Edge Mismatch by adding edge from A2 to A1)
            curr_key = make_key(root, b, a, g_type);
            p.first = 33;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;        

		case 18: // A1->R<-A2
			// A1<->R<-A2 (Edge Mismatch by adding edge from R to A1)
            curr_key = make_key(root, a, b, g_type);
            p.first = 15;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 15, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1->R<->A2 (Edge Mismatch by adding edge from R to A2)
            curr_key = make_key(root, b, a, g_type);
            p.first = 15;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 15, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1->R<-A2 and A1->A2 (Edge Mismatch by adding edge from A1 to A2)
            curr_key = make_key(root, a, b, g_type);
            p.first = 28;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1->R<-A2 and A1<-A2 (Edge Mismatch by adding edge from B to A)
            curr_key = make_key(root, b, a, g_type);
            p.first = 28;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 19: // A1<->R<->A2 and A1<->A2
			// A1->R<->A2 and A1<->A2 (Edge Mismatch by removing edge from R to A1)
            curr_key = make_key(root, b, a, g_type);
            p.first = 21;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<-R<->A2 and A1<->A2 (Edge Mismatch by removing edge from A1 to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 22;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1->R<->A2 and A1<->A2 (Edge Mismatch by removing edge from R to A2)
            curr_key = make_key(root, a, b, g_type);
            p.first = 21;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// A1<-R<->A2 and A1<->A2 (Edge Mismatch by removing edge from A2 to R)
            curr_key = make_key(root, a, b, g_type);
            p.first = 22;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<->R<->A2 and A1<-A2 (Edge Mismatch by removing edge from A1 to A2)
            curr_key = make_key(root, b, a, g_type);
            p.first = 20;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<->R<->A2 and A1->A2 (Edge Mismatch by removing edge from A2 to A1)
            curr_key = make_key(root, a, b, g_type);
            p.first = 20;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;
    
		case 20: // B<->R<->A and B<-A
			// B<->R<->A and B<->A (Edge Mismatch by adding edge from B to A)
			compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 19;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// B<->R<-A and B<-A (Edge Mismatch by removing edge from R to A)
            curr_key = make_key(root, b, a, g_type);
            p.first = 24;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// B<->R->A and B<-A (Edge Mismatch by removing edge from A to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 26;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 26, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// B->R<->A and B<-A (Edge Mismatch by removing edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 23;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 23, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// B<-R<->A and B<-A (Edge Mismatch by removing edge from B to R)
            curr_key = make_key(root, a, b, g_type);
            p.first = 25;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// B<->R<->A (Edge Mismatch by removing edge from A to B)
			compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 13;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 13, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 21: // A<->R<-B and A<->B
			// A<->R<->B and A<->B (Edge Mismatch by adding edge from R to B)
			compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 19;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R<-B and A<->B (Edge Mismatch by removing edge from R to A)
			compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 27;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 27, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R<-B and A<->B (Edge Mismatch by removing edge from A to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 29;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R<-B and A<-B (Edge Mismatch by removing edge from A to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 24;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// A<->R<-B and A->B (Edge Mismatch by removing edge from B to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 23;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 23, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A<->B (Edge Mismatch by removing edge from B to R)			
            curr_key = make_key(root, a, b, g_type);
            p.first = 4;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 4, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 22: // A<->R->B and A<->B
			// A<->R<->B and A<->B (Edge Mismatch by adding edge from R to B)
			compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 19;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R->B and A<->B (Edge Mismatch by removing edge from R to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 29;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R->B and A<->B (Edge Mismatch by removing edge from A to R)
			compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 32;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R->B and A<-B (Edge Mismatch by removing edge from A to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 26;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 26, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// A<->R->B and A->B (Edge Mismatch by removing edge from B to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 25;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A<->B (Edge Mismatch by removing edge from B to R)			
            curr_key = make_key(root, a, b, g_type);
            p.first = 4;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 4, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 23: // A<->R<-B and A->B
			// A<->R<->B and A->B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 20;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R<-B and A<->B (Edge Mismatch by adding edge from B to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 21;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R<-B and A->B (Edge Mismatch by removing edge from R to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 28;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R<-B and A->B (Edge Mismatch by removing edge from A to R)
            curr_key = make_key(root, a, b, g_type);
            p.first = 30;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A->B (Edge Mismatch by removing edge from B to R)			
            curr_key = make_key(root, a, b, g_type);
            p.first = 5;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 5, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R->B and A<-B (Edge Mismatch by removing edge from A to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 15;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 15, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 24: // A<->R<-B and A<-B
			// A<->R<->B and A<-B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, b, a, g_type);
            p.first = 20;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R<-B and A<->B (Edge Mismatch by adding edge from A to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 21;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R<-B and A<-B (Edge Mismatch by removing edge from R to A)
            curr_key = make_key(root, b, a, g_type);
            p.first = 28;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R<-B and A<-B (Edge Mismatch by removing edge from A to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 31;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A<-B (Edge Mismatch by removing edge from B to R)			
            curr_key = make_key(root, a, b, g_type);
            p.first = 6;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 6, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R->B and A<-B (Edge Mismatch by removing edge from B to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 15;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 15, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 25: // A<->R->B and A->B
			// A<->R<->B and A->B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 20;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R->B and A<->B (Edge Mismatch by adding edge from B to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 22;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R->B and A->B (Edge Mismatch by removing edge from R to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 31;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R->B and A->B (Edge Mismatch by removing edge from A to R)
            curr_key = make_key(root, a, b, g_type);
            p.first = 33;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A->B (Edge Mismatch by removing edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 5;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 5, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R->B and A<-B (Edge Mismatch by removing edge from A to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 14;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 26: // A<->R->B and A<-B
			// A<->R<->B and A<-B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, b, a, g_type);
            p.first = 20;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R->B and A<->B (Edge Mismatch by adding edge from A to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 22;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R->B and A<-B (Edge Mismatch by removing edge from R to A)
            curr_key = make_key(root, b, a, g_type);
            p.first = 30;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R->B and A<-B (Edge Mismatch by removing edge from A to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 33;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<->A<-B (Edge Mismatch by removing edge from R to B)			
            curr_key = make_key(root, a, b, g_type);
            p.first = 6;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 6, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<->R->B and A<-B (Edge Mismatch by removing edge from B to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 14;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 27: // A1->R<-A2 and A1<->A2
			// A1<->R<-A2 and A1<->A2 (Edge Mismatch by adding edge from R to A1)
            curr_key = make_key(root, a, b, g_type);
            p.first = 21;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1->R<->A2 and A1<->A2 (Edge Mismatch by adding edge from R to A2)
            curr_key = make_key(root, b, a, g_type);
            p.first = 21;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A2<->A1 (Edge Mismatch by removing edge from R to A1)
            curr_key = make_key(root, b, a, g_type);
            p.first = 10;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A1<->A2 (Edge Mismatch by removing edge from R to A2)
            curr_key = make_key(root, a, b, g_type);
            p.first = 10;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1->R<-A2 and A1<-A2 (Edge Mismatch by removing edge from A1 to A2)			
            curr_key = make_key(root, b, a, g_type);
            p.first = 28;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1->R<-A2 and A1->A2 (Edge Mismatch by removing edge from A2 to A1)
            curr_key = make_key(root, a, b, g_type);
            p.first = 28;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;
            
		case 28: // A->R<-B and A->B
			// A<->R<-B and A->B (Edge Mismatch by adding edge from R to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 23;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 23, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R<->B and A->B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, b, a, g_type);
            p.first = 24;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R<-B and A<->B (Edge Mismatch by adding edge from B to A)
			compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 27;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 27, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-B<-A (Edge Mismatch by removing edge from A to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 12;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 12, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A->B (Edge Mismatch by removing edge from B to R)
            curr_key = make_key(root, a, b, g_type);
            p.first = 11;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A->R<-B (Edge Mismatch by removing edge from A to B)
			compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 18;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 18, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;

		case 29: // B<-R<-A and B<->A
			// B<-R<->A and B->A (Edge Mismatch by adding edge from R to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 22;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// B<->R<-A and B->A (Edge Mismatch by adding edge from B to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 21;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->B<->A (Edge Mismatch by removing edge from A to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 7;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 7, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-A<->B (Edge Mismatch by removing edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 10;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// B<-R<-A and B->A (Edge Mismatch by removing edge from A to B)			
            curr_key = make_key(root, b, a, g_type);
            p.first = 30;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// B<-R<-A and B<-A (Edge Mismatch by removing edge from B to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 31;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 30: // A<-R<-B and A<-B
			// A<->R<-B and A<-B (Edge Mismatch by adding edge from A to R)
            curr_key = make_key(root, a, b, g_type);
            p.first = 23;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 23, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R<->B and A<-B (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, b, a, g_type);
            p.first = 26;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 26, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R<-B and A<->B (Edge Mismatch by adding edge from B to A)
            curr_key = make_key(root, b, a, g_type);
            p.first = 29;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R<-B<-A (Edge Mismatch by removing edge from R to A)
            curr_key = make_key(root, b, a, g_type);
            p.first = 12;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 12, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A->B (Edge Mismatch by removing edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 8;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 8, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R<-B (Edge Mismatch by removing edge from A to B)
            curr_key = make_key(root, b, a, g_type);
            p.first = 16;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 16, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 31: // B<-R<-A and B<-A
			// B<-R<->A and B<-A (Edge Mismatch by adding edge from R to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 25;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// B<->R<-A and B<-A (Edge Mismatch by adding edge from R to B)
            curr_key = make_key(root, b, a, g_type);
            p.first = 24;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// B<-R<-A and B<->A (Edge Mismatch by adding edge from B to A)
            curr_key = make_key(root, a, b, g_type);
            p.first = 29;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->B<-A (Edge Mismatch by removing edge from A to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 9;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 9, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// R<-A->B (Edge Mismatch by removing edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 11;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R<-B (Edge Mismatch by removing edge from A to B)			
            curr_key = make_key(root, a, b, g_type);
            p.first = 16;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 16, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 32: // A1<-R->A2 and A1<->A2
			// A1<->R->A2 and A1<->A2 (Edge Mismatch by adding edge from R to A1)
            curr_key = make_key(root, a, b, g_type);
            p.first = 22;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<-R<->A2 and A1->A2 (Edge Mismatch by adding edge from R to A2)
            curr_key = make_key(root, b, a, g_type);
            p.first = 22;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A2<-A1 (Edge Mismatch by removing edge from R to A1)
            curr_key = make_key(root, b, a, g_type);
            p.first = 7;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 7, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A1<->A2 (Edge Mismatch by removing edge from R to A2)
            curr_key = make_key(root, a, b, g_type);
            p.first = 7;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 7, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<-R->A2 and A1<-A2 (Edge Mismatch by removing edge from A1 to A2)
            curr_key = make_key(root, b, a, g_type);
            p.first = 33;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A1<-R->A2 and A1->A2 (Edge Mismatch by adding edge from A2 to A1)			
            curr_key = make_key(root, a, b, g_type);
            p.first = 33;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;
            
		case 33: // A<-R->B and A->B
			// A<->R->B and A->B (Edge Mismatch by adding edge from A to R)
            curr_key = make_key(root, a, b, g_type);
            p.first = 25;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R<->B and A->B (Edge Mismatch by adding edge from B to R)
            curr_key = make_key(root, b, a, g_type);
            p.first = 26;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 26, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R->B and A<->B (Edge Mismatch by adding edge from B to A)
			compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 32;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->B<-A (Edge Mismatch by removing edge from R to A)
            curr_key = make_key(root, b, a, g_type);
            p.first = 9;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 9, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// R->A->B (Edge Mismatch by removing edge from R to B)
            curr_key = make_key(root, a, b, g_type);
            p.first = 8;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 8, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// A<-R->B (Edge Mismatch by removing edge from A to B)
			compare_two(a, b, a1, b1);
            curr_key = make_key(root, a1, b1, g_type);
            p.first = 17;
            p.second = curr_key;
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, curr_key, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;
	}
}
