#include "digkernel.h"
#include "string.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <iomanip>


/*********************** DigraphKernel methods ***********************/
void DigraphKernel::read_graphs(string nlabels_file, string graph_file, const vector<unsigned> &vertices)  {
    if (VERBOSE)  cerr << "Reading input data ... ";

    digraph = SimpleGraph::read_digraph(nlabels_file.c_str(), graph_file.c_str());

    for (unsigned i=0; i<vertices.size(); i++)  {
        if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;

        roots.push_back(vertices[i]);
    }

    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::read_sim_matrix(string sim_matrix_file)  {
    if (VERBOSE)  cerr << "Reading probability similarity matrix for vertex labels file ... ";
    string column, row, key;
	
    // Read vertex labels similarity matrix.
    ifstream p(sim_matrix_file.c_str(), ios::in);
    if (p.fail()) {
        cerr << "ERROR: Vertex labels similarity matrix file " << sim_matrix_file << " cannot be opened." << endl; exit(1);
    }
    else  {
		if (VERBOSE)  cerr << sim_matrix_file.c_str();
        getline(p, column);
        while(getline(p, row))  {
            vector<string> tokens = split(row, '\t'); 
            for(unsigned i=0 ; i<column.size() ; i++)  {
                key = tokens[0] + column[i];
                map<string, float>::iterator it = sim_vlm_matrix.find(key);
                if(it == sim_vlm_matrix.end())  {
                    sim_vlm_matrix[key] = to_f(tokens[i+1]);
                }
            }
        }
    }
    p.close();

    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::compute_random_walk_cumulative_matrix(int steps, double restart)  {
    if (VERBOSE)  cerr << "Computing Cumulative Random Walk Digraph Kernel for #steps = " << steps << " restart prob = " << restart << "... ";

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);

    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;                                    

        kernel[i][i] = random_walk_cumulative(digraph, roots[i], roots[i], steps, restart);
        for (unsigned j=0; j<i; j++)  {			
            kernel[i][j] = random_walk_cumulative(digraph, roots[i], roots[j], steps, restart);
        }
    }
    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::compute_random_walk_matrix(int steps, double restart)  {
    if (VERBOSE)  cerr << "Computing Random Walk Digraph Kernel for #steps = " << steps << " restart prob = " << restart << "... ";

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);

    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;                                    

        kernel[i][i] = random_walk(digraph, roots[i], roots[i], steps, restart);
        for (unsigned j=0; j<i; j++)  {			
            kernel[i][j] = random_walk(digraph, roots[i], roots[j], steps, restart);
        }
    }
    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::compute_label_mismatch_matrix()  {
    if (VERBOSE)  {
        if (set_k((DIGRAPHLETS_TYPES-1), SF) > 0)
            cerr << "Computing Label Substitutions Digraphlet Kernel ... ";
        else
            cerr << "Computing Standard Digraphlet Kernel ... ";
    }

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);
    
	for (unsigned i=0; i<roots.size(); i++)  {
		for (unsigned j=0; j<=i; j++)
			kernel[i][j] = 0.0;
	}

    for (unsigned i=0; i<roots.size(); i++)  {
        hashes.push_back(get_graphlets_counts(digraph, roots[i]));
    }

	unsigned long g_type = 0;
	while (g_type < DIGRAPHLETS_TYPES)  {
		if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {        
            int VLM = set_k(g_type, SF);
            map<Key,MismatchInfo> mismatch_hash;

		    for (unsigned i=0; i<roots.size(); i++)  {
			    if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;

                add_vertex_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, VLM, false);

                update_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, false, VLM, false);

			    kernel[i][i] = kernel[i][i] + distance_hash_join(hashes[i][g_type], hashes[i][g_type], g_type);
			    for (unsigned j=0; j<i; j++)  {
				    kernel[i][j] = kernel[i][j] + distance_hash_join(hashes[i][g_type], hashes[j][g_type], g_type);
			    }
		    }
			vl_mismatch_neighborhood.clear();

			for (unsigned i=0; i<roots.size(); i++)  {
				hashes[i][g_type].clear();
            }
        }
        g_type = g_type + 1; 
	}

    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::compute_edge_mismatch_matrix()  {
    if (VERBOSE)   cerr << "Computing Edge Indels Digraphlet Kernel ... ";

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);

    for (unsigned i=0; i<roots.size(); i++)  {
        hashes.push_back(get_graphlets_counts(digraph, roots[i]));
        add_edge_mismatch_counts(hashes[i]);
    }

    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;
        
        kernel[i][i] = distance_hash_join(hashes[i], hashes[i]);
        for (unsigned j=0; j<i; j++)  {
            kernel[i][j] = distance_hash_join(hashes[i], hashes[j]);
        }
    }

    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::compute_edit_distance_matrix()  {
    if (VERBOSE)  cerr << "Computing Edit Distance Digraphlet Kernel (d=1) ... ";

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);
    
	for (unsigned i=0; i<roots.size(); i++)  {
		for (unsigned j=0; j<=i; j++)
			kernel[i][j] = 0.0;
	}

    for (unsigned i=0; i<roots.size(); i++)  {
        hashes.push_back(get_graphlets_counts(digraph, roots[i]));
        add_edge_mismatch_counts(hashes[i]);
    }

	unsigned long g_type = 0;
	while (g_type < DIGRAPHLETS_TYPES)  {
        if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
            int VLM = set_k(g_type, SF);
            map<Key,MismatchInfo> mismatch_hash;

		    for (unsigned i=0; i<roots.size(); i++)  {
			    if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;
                
                if (VLM >= 1)  {
                    add_vertex_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, 1, false);
                }

                update_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, false, 1, true);

			    kernel[i][i] = kernel[i][i] + distance_hash_join(hashes[i][g_type], hashes[i][g_type], g_type);
			    for (unsigned j=0; j<i; j++)  {
				    kernel[i][j] = kernel[i][j] + distance_hash_join(hashes[i][g_type], hashes[j][g_type], g_type);
			    }
		    }
			vl_mismatch_neighborhood.clear();

            for (unsigned i=0; i<roots.size(); i++)  {
                hashes[i][g_type].clear();
            }
        }
        g_type = g_type + 1; 
	}

    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::compute_edit_distance2_matrix()  {
    if (VERBOSE)  cerr << "Computing Edit Distance Digraphlet Kernel (d=2) ... ";

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);
    
	for (unsigned i=0; i<roots.size(); i++)  {
		for (unsigned j=0; j<=i; j++)
			kernel[i][j] = 0.0;
	}

    for (unsigned i=0; i<roots.size(); i++)  {
        hashes.push_back(get_graphlets_counts(digraph, roots[i]));
        add_1_edge_mismatch_counts(hashes[i]);
    }

	unsigned long g_type = 0;
	while (g_type < DIGRAPHLETS_TYPES)  {
        if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
            int VLM = set_k(g_type, SF);
            map<Key,MismatchInfo> mismatch_hash;

		    for (unsigned i=0; i<roots.size(); i++)  {
                if (VLM >= 1)  {
					add_vertex_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, 1, true);
                }
               
                update_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, true, 1, true);
            }
        }
		vl_mismatch_neighborhood.clear();
        g_type = g_type + 1;
    }
    for (unsigned i=0; i<roots.size(); i++)  {
        add_2_edge_mismatch_counts(hashes[i]);
    }

	g_type = 0;
	while (g_type < DIGRAPHLETS_TYPES)  {
        if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
            int VLM = set_k(g_type, SF);
            map<Key,MismatchInfo> mismatch_hash;

            for (unsigned i=0; i<roots.size(); i++)  {
                if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;

                if (VLM >= 2)  {
                    add_vertex_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, 2, false);
                }

                update_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, false, 2, true);

                kernel[i][i] = kernel[i][i] + distance_hash_join(hashes[i][g_type], hashes[i][g_type], g_type);
			    for (unsigned j=0; j<i; j++)  {
				    kernel[i][j] = kernel[i][j] + distance_hash_join(hashes[i][g_type], hashes[j][g_type], g_type);
			    }
		    }
			vl_mismatch_neighborhood.clear();

            for (unsigned i=0; i<roots.size(); i++)  {
                hashes[i][g_type].clear();
            }
        }
        g_type = g_type + 1; 
	}

    if (VERBOSE)  cerr << endl;
}

#if OUTPUT_FORMAT == 0
void DigraphKernel::write_matrix(const char *file)  {
    ofstream out(file, ios::out | ios::binary);
    
    unsigned g_size = roots.size();
    out.write((char*) &g_size, sizeof(unsigned)); 

    for (unsigned i=0; i<g_size; i++)
        for (unsigned j=0; j<=i; j++) 
            out.write((char*) &kernel[i][j], sizeof(float));
            
    out.close();
}
#elif OUTPUT_FORMAT == 1
void DigraphKernel::write_matrix(const char *file)  {
    ofstream out(file, ios::out );
    
    unsigned g_size = roots.size();
    
    for (unsigned i=0; i<g_size; i++)  {        
        for (unsigned j=0; j<=i; j++)  {
            out  << kernel[i][j] << "\t";
        }
        out  <<  endl;
    }
    out.close();
}
#else
void DigraphKernel::write_matrix(const char *file)  {
	ofstream out(file, ios::out );
	
	unsigned g_size = roots.size();
	
	for (unsigned i=0; i<g_size; i++)  {
		for (unsigned j=0; j<g_size; j++)   {
			if (j>i) {
                out  << setprecision (10) << kernel[j][i] <<  "\t"  ;
			} else  {
                out  << setprecision (10) << kernel[i][j]  <<  "\t"  ;
            }
        }
        out  <<  endl;
	}
	out.close();
}
#endif

void DigraphKernel::write_sparse_svml_lm(const char *file)  {
    if (VERBOSE)   { 
        if (set_k((DIGRAPHLETS_TYPES-1), SF) > 0)
            cerr << "Computing attributes for Label Substitutions Digraphlet Kernel ... ";
        else
            cerr << "Computing attributes for Standard Digraphlet Kernel ... ";
    }

    ofstream out(file, ios::out);
	unsigned g_type;
    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;

        out << labels[i];

        vector<map<Key,MismatchInfo> > g_hash = get_graphlets_counts(digraph, roots[i]);
        for (g_type=0; g_type<DIGRAPHLETS_TYPES; g_type++)  {
			if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
                int VLM= set_k(g_type, SF);
                map<Key,MismatchInfo> mismatch_hash;

                add_vertex_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, VLM, false);

                update_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, false, VLM, false);
                
				for (map<Key,MismatchInfo>::iterator it = g_hash[g_type].begin(); it != g_hash[g_type].end(); it++)
                    out << " " << get_feature_id(it->first, g_type) << ":" << retrieve_label_mismatch_count(g_type, g_hash[g_type], it->first);
			}
			vl_mismatch_neighborhood.clear();
        }
        out << " #" << i << endl;
    }
    out.close();

    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::write_sparse_svml_em(const char *file)  {
    if (VERBOSE)  cerr << "Computing attributes for Edge Indels Digraphlet Kernel ... ";

    ofstream out(file, ios::out);
    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;

        out << labels[i];

        vector<map<Key,MismatchInfo> > g_hash = get_graphlets_counts(digraph, roots[i]);
		add_edge_mismatch_counts(g_hash);
        for (unsigned g_type=0; g_type<DIGRAPHLETS_TYPES; g_type++)  {
            if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
                for (map<Key,MismatchInfo>::iterator it = g_hash[g_type].begin(); it != g_hash[g_type].end(); it++)  {
					out << " " << get_feature_id(it->first, g_type) << ":" << retrieve_edge_mismatch_count(g_hash[g_type], it->first);
                }
            }
        }
        out << " #" << i << endl;
    }
    out.close();

    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::write_sparse_svml_ed(const char *file)  {
    if (VERBOSE)  cerr << "Computing attributes for Edit Distance Digraphlet Kernel (d=1) ... ";

    ofstream out(file, ios::out);
    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;

        out << labels[i];

        vector<map<Key,MismatchInfo> > g_hash = get_graphlets_counts(digraph, roots[i]);
        add_edge_mismatch_counts(g_hash);
        for (unsigned g_type=0; g_type<DIGRAPHLETS_TYPES; g_type++)  {            
            if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
                int VLM = set_k(g_type, SF);
                map<Key,MismatchInfo> mismatch_hash;

                if(VLM >= 1)         
                    add_vertex_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, 1, false);
                    
                update_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, false, 1, true);

                for (map<Key,MismatchInfo>::iterator it = g_hash[g_type].begin(); it != g_hash[g_type].end(); it++)  {                    
                    if (retrieve_label_mismatch_count(g_type, g_hash[g_type], it->first) > 0.0)
                        out << " " << get_feature_id(it->first, g_type) << ":" << retrieve_label_mismatch_count(g_type, g_hash[g_type], it->first);
                }
            }
			vl_mismatch_neighborhood.clear();
        }
        out << " #" << i << endl;
    }
    out.close();

    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::write_sparse_svml_ed2(const char *file)  {
    if (VERBOSE)  cerr << "Computing attributes for Edit Distance Digraphlet Kernel (d=2) ... ";

    ofstream out(file, ios::out);
    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 1000 == 0)  cerr << " " << i+1;

        out << labels[i];

        vector<map<Key,MismatchInfo> > g_hash = get_graphlets_counts(digraph, roots[i]);
        add_1_edge_mismatch_counts(g_hash);
        for (unsigned g_type=0; g_type<DIGRAPHLETS_TYPES; g_type++)  {
            if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
                int VLM = set_k(g_type, SF);
                map<Key,MismatchInfo> mismatch_hash;

                if(VLM >= 1)                    
                    add_vertex_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, 1, true);
                                    
                update_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, true, 1, true);
                vl_mismatch_neighborhood.clear();

                add_2_edge_mismatch_counts(g_hash);

                if (VLM >= 2)
                    add_vertex_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, 2, false);

                update_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, false, 2, true);
                vl_mismatch_neighborhood.clear();

                for (map<Key,MismatchInfo>::iterator it = g_hash[g_type].begin(); it != g_hash[g_type].end(); it++)  {                    
                    if (retrieve_label_mismatch_count(g_type, g_hash[g_type], it->first) > 0.0)
                        out << " " << get_feature_id(it->first, g_type) << ":" << retrieve_label_mismatch_count(g_type, g_hash[g_type], it->first);
                }
            }
        }
        out << " #" << i << endl;
    }
    out.close();

    if (VERBOSE)  cerr << endl;
}

void DigraphKernel::write_labels(const char *file)  {
    unsigned num_pos(0), num_neg(0);

    ofstream out(file, ios::out);

    out << labels[0];
    labels[0]==1 ? num_pos++ : num_neg++;

    for (unsigned i=1; i<labels.size(); i++)  {
        out << "\n" << labels[i];
        labels[i]==1 ? num_pos++ : num_neg++;
    }
    out.close();

    if (VERBOSE)  {
        cerr << "Positive examples: " << num_pos << endl;
        cerr << "Negative examples: " << num_neg << endl;
    }
}

float DigraphKernel::random_walk_cumulative(SimpleGraph &g, unsigned g1_root, unsigned g2_root, int STEPS, double RESTART)  {
    unsigned i1, i2, j1, j2;
    unsigned g1_next_node, g2_next_node; 
	float value = compare_labels(g.nodes[g1_root], g.nodes[g2_root]);

    if(g.adj[g1_root].size() == 0 || g.adj[g2_root].size() == 0)  {
        return value;
    }

    i1 = g1_root;
    i2 = g2_root;
    int step(1);
    while(step < STEPS)  {
        if(g.adj[i1].size() == 0 || g.adj[i2].size() == 0)  {
            return value;
        }
        j1 = randint(g.adj[i1].size());
        j2 = randint(g.adj[i2].size());
        g1_next_node = g.adj[i1][j1];
        g2_next_node = g.adj[i2][j2];
        char node1_label = g.nodes[g1_next_node];
        char node2_label = g.nodes[g2_next_node];
        value = value + compare_labels(node1_label, node2_label);
        double prob = randdouble();
        // Restart random walk
        if(prob < RESTART)  {
			i1 = g1_root;
			i2 = g2_root;
        }
        // Continue random walk
        else  {
            i1 = g1_next_node;
            i2 = g2_next_node;
        }
        step += 1;
    }

    return value;
}

float DigraphKernel::random_walk(SimpleGraph &g, unsigned g1_root, unsigned g2_root, int STEPS, double RESTART)  {
    unsigned i1, i2, j1, j2;
    unsigned g1_next_node, g2_next_node; 
	float value(0.0);
    string g1_walk_seq, g2_walk_seq;

    if (g.nodes[g1_root] != g.nodes[g2_root] || g.adj[g1_root].size() == 0 || g.adj[g2_root].size() == 0)  {
        return value;
    }

    g1_walk_seq.push_back(g.nodes[g1_root]);
    g2_walk_seq.push_back(g.nodes[g2_root]);
    i1 = g1_root;
    i2 = g2_root;
    int step(1);
    while (step < STEPS)  {
        if (g.adj[i1].size() == 0 || g.adj[i2].size() == 0)  {
            return value;
        }
        j1 = randint(g.adj[i1].size());
        j2 = randint(g.adj[i2].size());
        g1_next_node = g.adj[i1][j1];
        g2_next_node = g.adj[i2][j2];
        g1_walk_seq.push_back(g.nodes[g1_next_node]);
        g2_walk_seq.push_back(g.nodes[g2_next_node]);
        double prob = randdouble();
        // Restart random walk
        if(prob < RESTART)  {
            i1 = g1_root;
            i2 = g2_root;
            if (g1_walk_seq.compare(g2_walk_seq) == 0) {
                value = value + 1.0;
            }
            g1_walk_seq.clear();
            g2_walk_seq.clear();
            g1_walk_seq.push_back(g.nodes[g1_root]);
            g2_walk_seq.push_back(g.nodes[g2_root]);
        }
        // Continue random walk
        else  {
            i1 = g1_next_node;
            i2 = g2_next_node;
        }
        step += 1;
    }

    return value;
}

// Count digraphlets starting from root.
vector<map<Key,MismatchInfo> > DigraphKernel::get_graphlets_counts(SimpleGraph &g, unsigned g_root)  {
    map<Key,MismatchInfo> T;
    vector<map<Key,MismatchInfo> > hash(DIGRAPHLETS_TYPES, T);
    vector<unsigned> dist = g.breadth_first_sort(g_root);
	vector<Key> mismatches;
	Key key;
    
	unsigned i, j;
    char root, a, b;

    root = g.nodes[g_root];
    
    if (DIGRAPHLETS_1)  {
		// 1-digraphlets, case 0, r
		key = create_permutations_subset(mismatches, root, ZERO_CHAR, ZERO_CHAR, 0);	
		increment_match_hash(hash[0], key, mismatches);
	}

    for (unsigned i_=0; i_<g.adj[g_root].size(); i_++)  {
        i = g.adj[g_root][i_];
        a = g.nodes[i];	
        bool found_ip(false);
        unsigned t(0);
        while (!found_ip && t<g.adj[i].size())  {
            found_ip = (g_root == g.adj[i][t++]);
        }
        // 2-digraphlets, cases r<->a or r->a
		if (DIGRAPHLETS_2)  {
            if (found_ip)  { //case 1, r<->a
				key = create_permutations_subset(mismatches, root, a, ZERO_CHAR, 1);
				increment_match_hash(hash[1], key, mismatches);
            }
            else  { //case 2, r->a
				key = create_permutations_subset(mismatches, root, a, ZERO_CHAR, 2);
				increment_match_hash(hash[2], key, mismatches);
            }
		}

	    if (DIGRAPHLETS_3)  {
            // 3-digraphlets, cases from a<->r<->b upto a<-r->b excluding all a->r or r<-b cases
            for (unsigned j_=0; j_<i_; j_++)  {
                j = g.adj[g_root][j_];
                b = g.nodes[j];
			    bool found_jp(false), found_ij(false), found_ji(false);
                unsigned t(0);
                while (!found_jp && t<g.adj[j].size())  {
                    found_jp = (g_root == g.adj[j][t++]);
                }

                t = 0;
				while (!found_ij && t<g.adj[i].size())  {
					found_ij = (j == g.adj[i][t++]);
                }

                t = 0;
				while (!found_ji && t<g.adj[j].size())  {
					found_ji = (i == g.adj[j][t++]);
				}

				if (found_ip)  {
                    if (found_jp)  {
                        if(found_ij)  {
                            if(found_ji)  { //case 19, a<->r<->b and a<->b
					            key = create_permutations_subset(mismatches, root, a, b, 19);
								increment_match_hash(hash[19], key, mismatches);
                            }
                            else  { //case 20, a<->r<->b and a->b
								key = create_permutations_subset(mismatches, root, a, b, 20);
								increment_match_hash(hash[20], key, mismatches);
                            }
                        }
                        else  {
                            if(found_ji)  { //case 20, a<->r<->b and a<-b
								key = create_permutations_subset(mismatches, root, b, a, 20);
								increment_match_hash(hash[20], key, mismatches); 
                            }
                            else  { //case 13, a<->r<->b
								key = create_permutations_subset(mismatches, root, a, b, 13);
								increment_match_hash(hash[13], key, mismatches);
                            }
                        }
                    }
                    else  {
                        if(found_ij)  {
                            if(found_ji)  { //case 22, a<->r->b and a<->b 
								key = create_permutations_subset(mismatches, root, a, b, 22);
								increment_match_hash(hash[22], key, mismatches);
                            }
                            else  { //case 25, a<->r->b and a->b
								key = create_permutations_subset(mismatches, root, a, b, 25);
								increment_match_hash(hash[25], key, mismatches);
                            }
                        }
                        else  {
                            if(found_ji)  { //case 26, a<->r->b and a<-b
								key = create_permutations_subset(mismatches, root, a, b, 26);
								increment_match_hash(hash[26], key, mismatches);                                
                            }
                            else  { //case 14, a<->r->b
								key = create_permutations_subset(mismatches, root, a, b, 14);
								increment_match_hash(hash[14], key, mismatches);
                            }
                        }
                    }
                }
				else  {
                    if (found_jp)  {
                        if(found_ij)  {
                            if(found_ji)  { //case 22, a<-r<->b and a<->b
								key = create_permutations_subset(mismatches, root, b, a, 22);
								increment_match_hash(hash[22], key, mismatches);
                            }
                            else  { //case 26, a<-r<->b and a->b
								key = create_permutations_subset(mismatches, root, b, a, 26);
								increment_match_hash(hash[26], key, mismatches);
                            }
                        }
                        else  {
                            if(found_ji)  { //case 25, a<-r<->b and a<-b
								key = create_permutations_subset(mismatches, root, b, a, 25);
								increment_match_hash(hash[25], key, mismatches);
                            }
                            else  {  //case 14, a<-r<->b                    
								key = create_permutations_subset(mismatches, root, b, a, 14);
								increment_match_hash(hash[14], key, mismatches);
                            }
                        }
                    }
                    else  {
                        if(found_ij)  {
                            if(found_ji)  { //case 32, a<-r->b and a<->b
								key = create_permutations_subset(mismatches, root, a, b, 32);
								increment_match_hash(hash[32], key, mismatches);
                            }
                            else  { //case 33, a<-r->b and a->b
								key = create_permutations_subset(mismatches, root, a, b, 33);
								increment_match_hash(hash[33], key, mismatches);
                            }
                        }
                        else  {
                            if(found_ji)  { //case 33, a<-r->b and a<-b
								key = create_permutations_subset(mismatches, root, b, a, 33);
								increment_match_hash(hash[33], key, mismatches);
                            }
                            else  { //case 17, a<-r->b
								key = create_permutations_subset(mismatches, root, a, b, 17);
								increment_match_hash(hash[17], key, mismatches);
                            }
                        }
                    }
                }
			}
        }

        if (DIGRAPHLETS_3)  {
            // 3-digraphlets, cases a<->r<-b or a<-r<-b
            for (j=0; j<g.nodes.size(); j++)  {
                for (unsigned j_=0; j_<g.adj[j].size(); j_++)  {
                    if (g.adj[j][j_] == g_root && j != i)  {
                        b = g.nodes[j];
                        bool found_pj(false), found_ij(false), found_ji(false);
                        unsigned t(0);
                        while (!found_pj && t<g.adj[g_root].size())  {
                            found_pj = (j == g.adj[g_root][t++]);
                        }
    
                        t = 0;
				        while (!found_ij && t<g.adj[i].size())  {
					        found_ij = (j == g.adj[i][t++]);
                        }

                        t = 0;
				        while (!found_ji && t<g.adj[j].size())  {
					        found_ji = (i == g.adj[j][t++]);
				        }

                        if(!found_pj)  {
                            if (found_ip)  {
                                if(found_ij)  {
                                    if(found_ji)  { //case 21, a<->r<-b and a<->b
										key = create_permutations_subset(mismatches, root, a, b, 21);
										increment_match_hash(hash[21], key, mismatches);
                                    }
                                    else  { //case 23, a<->r<-b and a->b
										key = create_permutations_subset(mismatches, root, a, b, 23);
										increment_match_hash(hash[23], key, mismatches);
                                    }
                                }
                                else  {
                                    if(found_ji)  { //case 24, a<->r<-b and a<-b
										key = create_permutations_subset(mismatches, root, a, b, 24);
										increment_match_hash(hash[24], key, mismatches);
                                    }
                                    else  { //case 15, a<->r<-b
										key = create_permutations_subset(mismatches, root, a, b, 15);
										increment_match_hash(hash[15], key, mismatches);
                                    }
                                }
                            }
                            else  {
                                if(found_ij)  {
                                    if(found_ji)  { //case 29, a<-r<-b and a<->b
										key = create_permutations_subset(mismatches, root, b, a, 29);
										increment_match_hash(hash[29], key, mismatches);
                                    }
                                    else  { //case 30, a<-r<-b and a->b
										key = create_permutations_subset(mismatches, root, a, b, 30);
										increment_match_hash(hash[30], key, mismatches);
                                    }
                                }
                                else  {
                                    if(found_ji)  { //case 31, a<-r<-b and a<-b
										key = create_permutations_subset(mismatches, root, b, a, 31);
										increment_match_hash(hash[31], key, mismatches);
                                    }
                                    else  { //case 16, a<-r<-b
										key = create_permutations_subset(mismatches, root, b, a, 16);
										increment_match_hash(hash[16], key, mismatches);
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
            }
        }

        if (DIGRAPHLETS_3)  {
		    // 3-digraphlets, cases r<->a<->b upto r->a->b excluding all r<-a and a<-b cases
		    for (unsigned j_=0; j_<g.adj[i].size(); j_++)  {
			    j = g.adj[i][j_];
			    b = g.nodes[j];
                if (j == g_root)  continue;
			    bool found_pj(false), found_jp(false), found_ji(false);
                unsigned t(0);
                while (!found_pj && t<g.adj[g_root].size())  {
                    found_pj = (j == g.adj[g_root][t++]);
                }
                
                t = 0;
                while (!found_jp && t<g.adj[j].size())  {
                    found_jp = (g_root == g.adj[j][t++]);
                }

                t = 0;
                while (!found_ji && t<g.adj[j].size())  {
                    found_ji = (i == g.adj[j][t++]);
                }

                if(!found_pj && !found_jp)  {
                    if (found_ip)  {
                        if (found_ji)  { //case 4, r<->a<->b
							key = create_permutations_subset(mismatches, root, a, b, 4);
							increment_match_hash(hash[4], key, mismatches);                    
                        }
                        else  { //case 5, r<->a->b
							key = create_permutations_subset(mismatches, root, a, b, 5);
							increment_match_hash(hash[5], key, mismatches);
                        }
                    }
                    else  {
                        if (found_ji)  { //case 7, r->a<->b
							key = create_permutations_subset(mismatches, root, a, b, 7);
							increment_match_hash(hash[7], key, mismatches);
                        }
                        else  { //case 8, r->a->b
							key = create_permutations_subset(mismatches, root, a, b, 8);
							increment_match_hash(hash[8], key, mismatches);
                        }
                    }
                }
		    }
        }
        
        if (DIGRAPHLETS_3)  {
            // 3-digraphlets, cases r<->a<-b or r->a<-b 
            for (j=0; j<g.nodes.size(); j++)  {
                for (unsigned j_=0; j_<g.adj[j].size(); j_++)  {
                    if (g.adj[j][j_] == i && 1 < dist[j])  {
                        b = g.nodes[j];
                        bool found_pj(false), found_jp(false), found_ij(false);
                        unsigned t(0);
                        while (!found_pj && t<g.adj[g_root].size())  {
                            found_pj = (j == g.adj[g_root][t++]);
                        }

                        t = 0;
                        while (!found_jp && t<g.adj[j].size())  {
                            found_jp = (g_root == g.adj[j][t++]);
                        }

                        t = 0;
                        while (!found_ij && t<g.adj[i].size())  {
                            found_ij = (j == g.adj[i][t++]);
                        }
                        
                        if (!found_pj && !found_jp && !found_ij)  {
                            if (found_ip)  { //case 6, r<->a<-b
								key = create_permutations_subset(mismatches, root, a, b, 6);
								increment_match_hash(hash[6], key, mismatches);
                            }
                            else  { //case 9, r->a<-b
								key = create_permutations_subset(mismatches, root, a, b, 9);
								increment_match_hash(hash[9], key, mismatches);
                            }
                        }
                        break;
                    }
                }
            }
        }
    }

    for (i=0; i<g.nodes.size(); i++)  {
        for (unsigned j_=0; j_<g.adj[i].size(); j_++)  {
            // 2-digraphlets, case r<-a
            if (g.adj[i][j_] == g_root)  {
                a = g.nodes[i];
                bool found_pi(false);
                unsigned t(0);
                while (!found_pi && t<g.adj[g_root].size())  {
                    found_pi = (i == g.adj[g_root][t++]);
                } 
                
		        if (!found_pi)  { //case 3, r<-a
                    if (DIGRAPHLETS_2)  {
						key = create_permutations_subset(mismatches, root, a, ZERO_CHAR, 3);
						increment_match_hash(hash[3], key, mismatches);
                    }
                    if (DIGRAPHLETS_3)  {
                        // 3-digraphlets, case r<-a<->b or r<-a->b
		                for (unsigned k_=0; k_<g.adj[i].size(); k_++)  {
			                j = g.adj[i][k_];
			                b = g.nodes[j];
                            if (k_ == j_ || j == g_root)  continue;
			                bool found_pj(false), found_jp(false), found_ji(false);
                            unsigned t(0);
                            while (!found_pj && t<g.adj[g_root].size())  {
                                found_pj = (j == g.adj[g_root][t++]);
                            }

                            t = 0;
                            while (!found_jp && t<g.adj[j].size())  {
                                found_jp = (g_root == g.adj[j][t++]);
                            }

                            t = 0;
                            while (!found_ji && t<g.adj[j].size())  {
                                found_ji = (i == g.adj[j][t++]);
                            }

                            if(!found_pj && !found_jp)  {
                                if (found_ji)  { //case 10, r<-a<->b
									key = create_permutations_subset(mismatches, root, a, b, 10);
									increment_match_hash(hash[10], key, mismatches);
                                }
                                else  { //case 11, r<-a->b
									key = create_permutations_subset(mismatches, root, a, b, 11);
									increment_match_hash(hash[11], key, mismatches);
                                }
                            }
                        }
                    }
                    if (DIGRAPHLETS_3)  { 
                        for (j=0; j<g.nodes.size(); j++)  {
                            for (unsigned k_=0; k_<g.adj[j].size(); k_++)  {
                                // 3-digraphlets, case r<-a<-b
                                if (g.adj[j][k_] == i && g.adj[j][k_] != g_root)  {
                                    b = g.nodes[j];
                                    bool found_pj(false), found_jp(false), found_ij(false);
                                    unsigned t(0);
                                    while (!found_pj && t<g.adj[g_root].size())  {
                                        found_pj = (j == g.adj[g_root][t++]);
                                    }

                                    t = 0;
                                    while (!found_jp && t<g.adj[j].size())  {
                                        found_jp = (g_root == g.adj[j][t++]);
                                    }
                                    t = 0;
                                    while (!found_ij && t<g.adj[i].size())  {
                                        found_ij = (j == g.adj[i][t++]);
                                    }

                                    if(!found_pj && !found_jp)  {
                                        if (!found_ij)  { //case 12, r<-a<-b
											key = create_permutations_subset(mismatches, root, a, b, 12);
											increment_match_hash(hash[12], key, mismatches);
                                        }
                                    }
		                            break;
                                }
                            }
                        }
                    }
                    if (DIGRAPHLETS_3)  {
                        for (j=0; j<i; j++)  {
                            for (unsigned j_=0; j_<g.adj[j].size(); j_++)  {
                                // 3-digraphlets, all cases a->r<-b
                                if (g.adj[j][j_] == g_root)  {
                                    b = g.nodes[j];
                                    bool found_pj(false), found_ij(false), found_ji(false);
                                    unsigned t(0);
                                    while (!found_pj && t<g.adj[g_root].size())  {
                                        found_pj = (j == g.adj[g_root][t++]);
                                    }

                                    t = 0;
                                    while (!found_ij && t<g.adj[i].size())  {
                                        found_ij = (j == g.adj[i][t++]);
                                    }

                                    t = 0;
                                    while (!found_ji && t<g.adj[j].size())  {
                                        found_ji = (i == g.adj[j][t++]);
                                    }

                                    if (!found_pj)  {
                                        if (found_ij)  {
                                            if (found_ji)  { //case 27, a->r<-b and a<->b
												key = create_permutations_subset(mismatches, root, a, b, 27);
												increment_match_hash(hash[27], key, mismatches);
                                            }
                                            else  { //case 28, a->r<-b and a->b
												key = create_permutations_subset(mismatches, root, a, b, 28);
												increment_match_hash(hash[28], key, mismatches);
                                            }
                                        }
                                        else  {
                                            if (found_ji)  { //case 28, a->r<-b and a<-b
												key = create_permutations_subset(mismatches, root, b, a, 28);
												increment_match_hash(hash[28], key, mismatches);
                                            }
                                            else  { //case 18, a->r<-b
												key = create_permutations_subset(mismatches, root, a, b, 18);
												increment_match_hash(hash[18], key, mismatches);
                                            }
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
                break;
            }
        }
    }

    if (NORMALIZE)
        normalize_spectral(hash);

    return hash;
}

// Add inexact graphlets by allowing vertex and edge label mismatches upto VLM.
void DigraphKernel::add_vertex_label_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, unsigned long g_type, int VLM, bool option)  {
	// Update counts to include vertex label mismacthes.
	if(VLM > 0 && ((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3)))  {
        // For each exact graphlet, generate all mismatch graphlets upto vertex label distance VLM.
        for (map<Key,MismatchInfo>::iterator git = hash.begin(); git != hash.end(); git++)  {
            if (hash[git->first].matches > 0.0)
                generate_vertex_label_mismatch_graphlets(vl_mismatch_neighborhood, hash, mismatch_hash, git->first, g_type, ALPHABET_ROOT, ALPHABET, sim_vlm_matrix, VLM);
            else  {
                if (option)
                    generate_vertex_label_mismatch_graphlets(vl_mismatch_neighborhood, hash, mismatch_hash, git->first, g_type, ALPHABET_ROOT, ALPHABET, sim_vlm_matrix, VLM);
            }
        }
    }
}

void DigraphKernel::update_label_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, unsigned long g_type, bool option, int VLM, bool eq)  { 
	if(VLM > 0 && ((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3)))  {
        // Update mismatch counts for newly created graphlets.
        for (map<Key,MismatchInfo>::iterator git = hash.begin(); git != hash.end(); git++)  {
            if (option)  {
                update_mismatch_count(hash, git->first, (hash[git->first].matches + hash[git->first].mismatches), g_type, sim_vlm_matrix, VLM, eq);
                update_mismatch_count(mismatch_hash, git->first, (hash[git->first].matches + hash[git->first].mismatches), g_type, sim_vlm_matrix, VLM, eq);
            }
            else  {
                update_mismatch_count(hash, git->first, hash[git->first].matches, g_type, sim_vlm_matrix, VLM, eq);
                update_mismatch_count(mismatch_hash, git->first, hash[git->first].matches, g_type, sim_vlm_matrix, VLM, eq);
            }
        }
        // Merge all non-zero graphlets into one hashable list. 
        for (map<Key,MismatchInfo>::iterator git = mismatch_hash.begin(); git != mismatch_hash.end(); git++)  {
            insert_mismatch_counts(hash, mismatch_hash, git->first);
        }
        mismatch_hash.clear();
    }
}

// Add inexact graphlets by allowing edge insertions and deletions upto EM.
void DigraphKernel::add_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash)  {
    map<Key,MismatchInfo> T;
    vector<map<Key,MismatchInfo> > mismatch_hash(DIGRAPHLETS_TYPES, T);
    unsigned vindex = 0;
    list<pair <unsigned long, Key> > L;
    float mult_factor;

    if (EM > 0)  {
        for (unsigned g_type=0; g_type<DIGRAPHLETS_TYPES; g_type++)  {
            if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
                for (map<Key,MismatchInfo>::iterator it = hash[g_type].begin(); it != hash[g_type].end(); it++)  { 
                    vindex = 0;
                    vector<list<pair <unsigned long, Key> > > EM_set(EM, L);
                    pair <unsigned long, Key> p (g_type, it->first); 
                    mult_factor = retrieve_exact_matches_count(hash[g_type], it->first);

                    // Generate edge mismacth neighborhood for this graphlet.
                    update_edge_mismatch_count(EM_set, it->first, g_type, EM, vindex);

                    list<pair <unsigned long, Key> > :: iterator list_it, search_it;
                    bool already_added(false);

                    // Add graphlets within edge mismacth neighborhood of current graphlet.
                    for (unsigned em = 0; em < EM; em++)  { 
                        EM_set[em].sort();
                        for (list_it = EM_set[em].begin(); list_it != EM_set[em].end(); list_it++)  {
                            already_added = false;
                            for (int previous = em-1; previous >= 0 ; previous--)  {
                                for (search_it = EM_set[previous].begin(); search_it != EM_set[previous].end(); search_it++)  {
                                    if ((list_it->first == p.first && list_it->second == p.second) || (list_it->first == search_it->first && list_it->second == search_it->second))  {
                                        already_added = true;
                                        break;
                                    }
                                }
                            }
                            if (!already_added)  {
                                vector<Key> mismatches;
                                char root, a, b;
			                    initialize_vertices_labels(list_it->second, root, a, b);
			                    Key k = create_permutations_subset(mismatches, root, a, b, list_it->first);
                                if (k != list_it->second)  {
                                    cerr << "ERROR: Digraphlets keys do not match " << print_key(k) << " vs " << print_key(list_it->second) << endl; exit(1);
                                }
                                increment_edge_mismatch_hash(hash[list_it->first], mismatch_hash[list_it->first], list_it->second, mult_factor, mismatches);
                            }
                        }
                    }
                }
            }
        }
        for (unsigned g_type=0; g_type<DIGRAPHLETS_TYPES; g_type++)  {
            if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
                for (map<Key,MismatchInfo>::iterator mit = mismatch_hash[g_type].begin(); mit != mismatch_hash[g_type].end(); mit++)  {
					insert_mismatch_counts(hash[g_type], mismatch_hash[g_type], mit->first);					
                }
            }
        }
        mismatch_hash.clear();
    }
}

// Add inexact graphlets by allowing 1-edge insertion and deletion.
void DigraphKernel::add_1_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash)  {
    map<Key,MismatchInfo> T;
    vector<map<Key,MismatchInfo> > mismatch_hash(DIGRAPHLETS_TYPES, T);
    unsigned vindex = 0;
    list<pair <unsigned long, Key> > L;
    float mult_factor;

    if (EM > 0)  {
        for (unsigned g_type=0; g_type<DIGRAPHLETS_TYPES; g_type++)  {
            if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
                for (map<Key,MismatchInfo>::iterator it = hash[g_type].begin(); it != hash[g_type].end(); it++)  { 
                    vindex = 0;
                    vector<list<pair <unsigned long, Key> > > EM_set(EM, L);
                    pair <unsigned long, Key> p (g_type, it->first); 
                    mult_factor = retrieve_exact_matches_count(hash[g_type], it->first);

                    // Generate edge mismacth neighborhood for this graphlet.
                    update_edge_mismatch_count(EM_set, it->first, g_type, 1, vindex);

                    list<pair <unsigned long, Key> > :: iterator list_it;

                    // Add graphlets within edge mismacth neighborhood of current graphlet.
                    EM_set[0].sort();
					for (list_it = EM_set[0].begin(); list_it != EM_set[0].end(); list_it++)  {
                        vector<Key> mismatches;
                        char root, a, b;
			            initialize_vertices_labels(list_it->second, root, a, b);
			            Key k = create_permutations_subset(mismatches, root, a, b, list_it->first);
                        if (k != list_it->second)  {
                            cerr << "ERROR: Digraphlets keys do not match " << print_key(k) << " vs " << print_key(list_it->second) << endl; exit(1);
                        }
                        increment_edge_mismatch_hash(hash[list_it->first], mismatch_hash[list_it->first], list_it->second, mult_factor, mismatches);                
                    }
                }
            }
        }
        for (unsigned g_type=0; g_type<DIGRAPHLETS_TYPES; g_type++)  {
            if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
                for (map<Key,MismatchInfo>::iterator mit = mismatch_hash[g_type].begin(); mit != mismatch_hash[g_type].end(); mit++)  {                    
					insert_mismatch_counts(hash[g_type], mismatch_hash[g_type], mit->first);
                }
            }
        }
        mismatch_hash.clear();
    }
}

// Add inexact graphlets by allowing 2-edge insertions and deletions.
void DigraphKernel::add_2_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash)  {
    map<Key,MismatchInfo> T;
    vector<map<Key,MismatchInfo> > mismatch_hash(DIGRAPHLETS_TYPES, T);
    unsigned vindex = 0;
    list<pair <unsigned long, Key> > L;
    float mult_factor;

    if (EM > 0)  {
        for (unsigned g_type=0; g_type<DIGRAPHLETS_TYPES; g_type++)  {
            if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
                for (map<Key,MismatchInfo>::iterator it = hash[g_type].begin(); it != hash[g_type].end(); it++)  { 
                    vindex = 0;
                    vector<list<pair <unsigned long, Key> > > EM_set(EM, L);
                    pair <unsigned long, Key> p (g_type, it->first); 
                    mult_factor = retrieve_exact_matches_count(hash[g_type], it->first);

                    // Generate edge mismacth neighborhood for this graphlet.
                    update_edge_mismatch_count(EM_set, it->first, g_type, 2, vindex);

                    list<pair <unsigned long, Key> > :: iterator list_it, search_it;
                    bool already_added(false);

                    // Add graphlets within edge mismacth neighborhood of current graphlet.
                    for (unsigned em = 1; em < EM; em++)  { 
                        EM_set[em].sort();
                        for (list_it = EM_set[em].begin(); list_it != EM_set[em].end(); list_it++)  {
                            already_added = false;
                            for (int previous = em-1; previous >= 0 ; previous--)  {
                                for (search_it = EM_set[previous].begin(); search_it != EM_set[previous].end(); search_it++)  {
                                    if ((list_it->first == p.first && list_it->second == p.second) || (list_it->first == search_it->first && list_it->second == search_it->second))  {
                                        already_added = true;
                                        break;
                                    }
                                }
                            }
                            if (!already_added)  {
                                vector<Key> mismatches;
                                char root, a, b;
			                    initialize_vertices_labels(list_it->second, root, a, b);
			                    Key k = create_permutations_subset(mismatches, root, a, b, list_it->first);
                                if (k != list_it->second)  {
                                    cerr << "ERROR: Digraphlets keys do not match " << print_key(k) << " vs " << print_key(list_it->second) << endl; exit(1);
                                }
                                increment_edge_mismatch_hash(hash[list_it->first], mismatch_hash[list_it->first], list_it->second, mult_factor, mismatches);
                            }
                        }
                    }
                }
            }
        }
        for (unsigned g_type=0; g_type<DIGRAPHLETS_TYPES; g_type++)  {
            if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
                for (map<Key,MismatchInfo>::iterator mit = mismatch_hash[g_type].begin(); mit != mismatch_hash[g_type].end(); mit++)  {
					insert_mismatch_counts(hash[g_type], mismatch_hash[g_type], mit->first);
                }
            }
        }
        mismatch_hash.clear();
    }
}


void DigraphKernel::normalize_spectral(map<Key,MismatchInfo> &hash, unsigned long g_type)  {
    float norm = sqrt(distance_hash_join(hash, hash, g_type));    
    if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {
		for (map<Key,MismatchInfo>::iterator it = hash.begin(); it != hash.end(); it++)  {
			it->second.matches = retrieve_label_mismatch_count(g_type, hash, it->first)/norm;
		}
	}
}

float DigraphKernel::distance_hash_join(map<Key,MismatchInfo> g_hash, map<Key,MismatchInfo> h_hash, unsigned long g_type)  {
    float sum(0);
	
	if((g_type == 0 && (DIGRAPHLETS_1)) || ((g_type >= 1 && g_type <= 3) && DIGRAPHLETS_2) || ((g_type >= 4 && g_type <= 33) && DIGRAPHLETS_3))  {	
		for (map<Key,MismatchInfo>::iterator git = g_hash.begin(); git != g_hash.end(); git++)  {
			map<Key,MismatchInfo>::iterator hit = h_hash.find(git->first);
			if (hit != h_hash.end())  {
				if (NORMALIZE)  {
					sum += git->second.matches * hit->second.matches;
                }
				else  {                
					sum += retrieve_label_mismatch_count(g_type, g_hash, git->first) * retrieve_label_mismatch_count(g_type, h_hash, hit->first);
                }
			}
		}
	}
    return sum;
}

void DigraphKernel::normalize_spectral(vector<map<Key,MismatchInfo> > &hash)  {
    float norm = sqrt(distance_hash_join(hash, hash));
    
    for (unsigned i=0; i<DIGRAPHLETS_TYPES; i++)  {
        for (map<Key,MismatchInfo>::iterator it = hash[i].begin(); it != hash[i].end(); it++)
            it->second.matches = retrieve_edge_mismatch_count(hash[i], it->first)/norm;
    }
}

float DigraphKernel::distance_hash_join(vector<map<Key,MismatchInfo> > g_hash, vector<map<Key,MismatchInfo> > h_hash)  {
    float sum(0);
    
    for (unsigned i=0; i<DIGRAPHLETS_TYPES; i++)  {
        if((i == 0 && (DIGRAPHLETS_1)) || ((i >= 1 && i <= 3) && DIGRAPHLETS_2) || ((i >= 4 && i <= 33) && DIGRAPHLETS_3))  {	
            for (map<Key,MismatchInfo>::iterator git = g_hash[i].begin(); git != g_hash[i].end(); git++)  {
                map<Key,MismatchInfo>::iterator hit = h_hash[i].find(git->first);
            
                if (hit != h_hash[i].end())  {
                    if (NORMALIZE)
                        sum += git->second.matches * hit->second.matches;
                    else
                        sum += retrieve_edge_mismatch_count(g_hash[i], git->first) * retrieve_edge_mismatch_count(h_hash[i], hit->first);
                }
            }
        }
    }
    return sum;
}

