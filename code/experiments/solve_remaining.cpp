/*
 * File:            solve_remaining.cpp
 * Description:     This file uses the depth_first_bnb algorithm in parallel to
 *                  compute folding statistics for the passed proteins.
 *                  The computed statistics are:
 *                      - Time it took to compute exact minimum of the score.
 *                      - Minimum found score.
 *                      - Number of placements for an amino acid.
 *                      - Number of full conformations checked.
 *
 * NOTE: Never run this experiment with more MPI ranks than proteins to solve!
 */

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <vector>

/* Include the algorithm from Prospr. */
#include "../prospr/prospr/core/src/depth_first_bnb.hpp"
#include "../prospr/prospr/core/src/protein.hpp"

/* Set namespaces. */
namespace fs = std::filesystem;
using namespace std;

/* Define constants. */
string RESFILE_PREFIX = "dfs_bnb";

/* Load protein IDs of already folded proteins as keys in a hashmap. */
map<tuple<string,string,string>, int> _load_finished_proteins(
        int cluster_size, string results_file) {
    ifstream cur_file;
    string line, word;
    stringstream linestream;
    map<tuple<string,string,string>, int> finished_proteins;

    for (int i = 0; i < cluster_size; i++) {
        /* Check if rank previously finished fully or partially. */
        if (fs::exists(results_file + "_r" + to_string(i) + ".csv")) {
            cur_file = ifstream(results_file + "_r" + to_string(i) + ".csv");
        } else if (fs::exists(results_file + "_r" + to_string(i) + "_TMP.csv")) {
            cur_file = ifstream(results_file + "_r" + to_string(i) + "_TMP.csv");
        }

        /* Process finished proteins, if there are any. */
        if (cur_file.is_open()) {
            /* Skip header and read protein IDs. */
            getline(cur_file, line);

            /* Add every line to the map of finished proteins. */
            while (std::getline(cur_file, line)) {
                linestream = stringstream(line);

                /* Read all columns, but use length,hratio,protein_id as key. */
                std::vector<std::string> cols;
                while (std::getline(linestream, word, ','))  {
                    cols.push_back(word);
                }
                std::tuple<string,string,string> key{ cols[0], cols[1], cols[2] };
                finished_proteins[key] = 1;
            }

            /* Close and clear filestream for next rank. */
            cur_file.close();
            cur_file.clear();
        }
    }

    return finished_proteins;
}

/* Load proteins from the CSV for this rank. */
vector<tuple<string,string,string,string>> load_proteins(string csv_path, 
        int cluster_size, int rank, fs::path results_path) {
    vector<tuple<string,string,string,string>> proteins;
    ifstream cur_file;
    string line, word;
    stringstream linestream;
    string protein_seq;
    string results_file = results_path.u8string() + "/" + RESFILE_PREFIX;

    /* Check if all proteins have already been finished. */
    if (fs::exists(results_file + ".csv")) {
        return proteins;
    }

    /* Get finished protein IDs in a map for efficient lookups. */
    map<tuple<string,string,string>, int> finished_proteins =
        _load_finished_proteins(cluster_size, results_file);

    /* Open input CSV file as stream and skip header. */
    cur_file = ifstream(csv_path, ios::in);
    getline(cur_file, line);

    while (std::getline(cur_file, line)) {
        linestream = stringstream(line);

        /* Read all columns and construct key for indexing. */
        std::vector<std::string> cols;
        while (std::getline(linestream, word, ','))  {
            cols.push_back(word);
        }
        std::tuple<string,string,string> key{ cols[0], cols[1], cols[2] };

        /* Verify that the protein has not yet been solved, and add to list. */
        if (auto it = finished_proteins.find(key); it == finished_proteins.end()) {
            proteins.push_back({ cols[0], cols[1], cols[2], cols[3] });
        }
    }

    cur_file.close();

    /* Determine what proteins to load for this specific node. */
    int begin = rank * (proteins.size() / cluster_size);
    int ranks_extra = proteins.size() % cluster_size;
    begin += min(rank, ranks_extra);

    int num_proteins = proteins.size() / cluster_size;
    if (rank < ranks_extra) {
        num_proteins++;
    }

    /* Slice and return proteins. */
    proteins = vector<tuple<string,string,string,string>>(
        proteins.begin() + begin, proteins.begin() + begin + num_proteins
    );

    return proteins;
}

/* Solves the given proteins sequentially using dfs_bnb. */
void dfs_bnb_parallel(string results_path, 
        vector<tuple<string,string,string,string>> proteins, int dim, int rank) {
    vector<int> min_conf;
    Protein *p;
    ofstream output_file;
    string results_file;

    chrono::time_point<chrono::high_resolution_clock> t1;
    chrono::time_point<chrono::high_resolution_clock> t2;
    chrono::duration<double> duration;

    /* Check if TMP results already exist for this job, append if so. 
     * Create new TMP file and write header otherwise.
     */
    results_file = results_path + "/" + RESFILE_PREFIX + "_r" + to_string(rank);
    if (fs::exists(results_file + "_TMP.csv")) {
        output_file.open(results_file + "_TMP.csv", ios::app);
    } else {
        output_file.open(results_file + "_TMP.csv", ios::out);
        output_file << "length,hratio,protein_id,algorithm,time,score,checked,"
                    << "placed,hash\n" << flush;
    }

    /* Fold each Protein and store resulting minimum energy conformation. */
    for (tuple<string,string,string,string> cur_p : proteins) {
        p = new Protein(get<3>(cur_p), dim, "HP");
        t1 = chrono::high_resolution_clock::now();
        depth_first_bnb(p, "reach_prune");
        t2 = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(t2 - t1);
        min_conf = p->hash_fold();

        /* Write statistics to output file. */
        output_file << get<0>(cur_p) << "," << get<1>(cur_p) << "," 
                    << get<2>(cur_p) << ",dfs_bnb_reach_prune," 
                    << duration.count() << "," << p->get_score() << "," 
                    << p->get_solutions_checked() 
                    << "," << p->get_aminos_placed() << "," << flush;

        /* Write conformation to file and close with newline. */
        output_file << "\"[";
        for (int b : vector<int>(min_conf.begin(), min_conf.end() - 1)) {
            output_file << b << ", " << flush;
        }
        output_file << min_conf.back() << "]\"\n" << flush;

        /* Deallocate the protein. */
        delete p;
    }

    /* When finished, close temporary rank files. */
    output_file.close();
    output_file.clear();

    /* Append TMP file to existing rank results, or simply rename file if none
     * exist yet. 
     */
    if (fs::exists(results_file + ".csv")) {
        ifstream tmp_results;
        string line;

        /* Open existing results and TMP file to read from. */
        output_file.open(results_file + ".csv", ios::app);
        tmp_results.open(results_file + "_TMP.csv", ios::in);

        /* Skip header and copy results from TMP to final results. */
        getline(tmp_results, line);
        output_file << tmp_results.rdbuf() << flush;

        /* Close file streams and remove TMP file. */
        tmp_results.close();
        output_file.close();
        remove((results_file + "_TMP.csv").c_str());
    } else {
        /* No existing results, so rename TMP file to rank results file. */
        rename((results_file + "_TMP.csv").c_str(),
            (results_file + ".csv").c_str());
    }
}

/* Merge results of all ranks in one final results file. */
void merge_results(fs::path results_path, int cluster_size) {
    string fin_results_path;
    string rank_results_path;
    ofstream fin_results_file;
    ifstream rank_results_file;
    struct stat buffer;
    string line;

    /* Open final results file and write header. */
    fin_results_path = results_path.u8string() + "/" + RESFILE_PREFIX + ".csv";
    fin_results_file.open(fin_results_path);
    fin_results_file << "length,hratio,protein_id,algorithm,time,score,"
                     << "checked,placed,hash\n"
                     << flush;

    /* Log path of merged results. */
    string node_debug = "[1/" + to_string(cluster_size) + 
        "]\tWriting results to:\t" + fin_results_path + "\n";
    cerr.write(node_debug.data(), node_debug.size());                       
    
    /* Fill final results file with results for each rank. */
    for (int rank = 0; rank < cluster_size; rank++) {
        /* Setup results filename according to whether it fully finished. */
        rank_results_path = results_path.u8string() + "/" +  RESFILE_PREFIX + 
            "_r" + to_string(rank);

        /* Use rank's TMP file if it still exists. */
        if (stat((rank_results_path + "_TMP.csv").c_str(), &buffer) == 0) {
            rank_results_path += "_TMP.csv";
        } else {
            rank_results_path += ".csv";
        }

        /* Read rank results, skip header, and write to final results. */
        rank_results_file.open(rank_results_path);
        getline(rank_results_file, line);
        fin_results_file << rank_results_file.rdbuf() << flush;
        rank_results_file.close();
        rank_results_file.clear();
    }

    /* Close final results file. */
    fin_results_file.close();
}

/* Entry point for experiments. */
int main(int argc, char *argv[]) {
    /* Experiment must be called with at least 4 arguments. */
    if (argc < 4) {
        cout << "Error: argc < 4..\n" << flush;
        exit(-1);
    }

    /* Determine rank and cluster size. */
    int cluster_size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &cluster_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Parse arguments. */
    fs::path results_path = fs::path(argv[1]);
    string csv_path = fs::path(argv[2]);
    int dim = atoi(argv[3]);

    /* Print general debug info if this is rank 0. */
    if (rank == 0) {
        cout << "Debug Info:\n";
        cout << "\tResults path:         " << results_path << "\n";
        cout << "\tInput CSV path:       " << csv_path << "\n";
        cout << "\tDimension of fold:    " << dim << "\n";
        cout << "\n\n";
    }

    /* Create results directory if it does not exist. */
    if (!fs::is_directory(results_path)) {
        fs::create_directories(results_path);
    }

    /* Notify creation of each rank. */
    string node_debug = "[" + to_string(rank + 1) + "/" +
        to_string(cluster_size) + "]\tI'm online..\n";
    cerr.write(node_debug.data(), node_debug.size());

    /* Run and time experiment with given arguments per specified length. */
    chrono::time_point<chrono::high_resolution_clock> t1;
    chrono::time_point<chrono::high_resolution_clock> t2;
    chrono::duration<double> duration;
    ofstream time_file;

    /* Load the proteins for this rank. */
    vector<tuple<string,string,string,string>> proteins = load_proteins(
        csv_path, cluster_size, rank, results_path
    );

    /* Time solving the set of proteins for this rank. */
    t1 = chrono::high_resolution_clock::now();
    dfs_bnb_parallel(results_path, proteins, dim, rank);
    t2 = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(t2 - t1);

    /* Each rank logs runtime when finished. */
    node_debug = "[" + to_string(rank + 1) + "/" + to_string(cluster_size) + 
        "]\tFinished my work in:\t" + to_string(duration.count()) + " us.\n";
    cerr.write(node_debug.data(), node_debug.size());

    /* Wait for everyone to finish, then rank 0 writes total runtime. */
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        t2 = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::microseconds>(t2 - t1);
        cout << "Finished " << " in " << duration.count() << "\n";
        time_file.open(results_path.u8string() + "/time_solve_remaining_us");
        time_file << duration.count() << "\n" << flush;
        time_file.close();
    }

    /* Close MPI processes. */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    node_debug = "[" + to_string(rank + 1) + "/" + to_string(cluster_size) +
                "]\tFinalized MPI..\n";
    cerr.write(node_debug.data(), node_debug.size());

    /* Merge results into one file per protein length. */
    if (rank == 0) {
        node_debug = "[" + to_string(rank + 1) + "/" + to_string(cluster_size) +
                "]\tMerging results..\n";
        cerr.write(node_debug.data(), node_debug.size());
        merge_results(results_path, cluster_size);
    }

    return 0;
}
