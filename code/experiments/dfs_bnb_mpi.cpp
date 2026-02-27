/*
 * File:            dfs_bnb_mpi.cpp
 * Description:     This file uses the depth_first_bnb algorithm in parallel to
 *                  compute folding statistics for the specified lengths of the 
 *                  specified dataset.
 *                  The computed statistics are:
 *                      - Time it took to compute exact minimum of the score.
 *                      - Minimum found score.
 *                      - Number of placements for an amino acid.
 *                      - Number of full conformations checked.
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

namespace fs = std::filesystem;
using namespace std;

/* Load protein IDs of already folded proteins as keys in a hashmap. */
map<int, int> _load_finished_proteins(int cluster_size, string results_file) {
  ifstream cur_file;
  string line, word;
  stringstream linestream;
  int protein_id;
  map<int, int> finished_ids;

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

      while (getline(cur_file, line)) {
        linestream = stringstream(line);
        getline(linestream, word, ',');
        protein_id = stoi(word);
        finished_ids[protein_id] = 1;
      }

      /* Close and clear filestream for next rank. */
      cur_file.close();
      cur_file.clear();
    }
  }

  return finished_ids;
}

/* Load proteins from dataset for each length. */
vector<pair<int, string>> load_proteins(string dataset_path, string setup, 
                                        int cluster_size, int rank, 
                                        fs::path results_path) {
  vector<pair<int, string>> proteins;
  ifstream cur_file;
  string line, word;
  stringstream linestream;
  int protein_id;
  string protein_seq;
  string results_file =
      results_path.u8string() + "/HP_" + setup + "_dfs_bnb";

  /* Check if all proteins have already been finished. */
  if (fs::exists(results_file + ".csv")) {
    return proteins;
  }

  /* Get finished protein IDs in a map.
   * Use underlying red-black tree for efficient lookups.
   */
  map<int, int> finished_ids =
      _load_finished_proteins(cluster_size, results_file);

  /* Open file as stream and skip header. */
  cur_file =
      ifstream(dataset_path + "/" + setup + ".csv", ios::in);
  getline(cur_file, line);

  /* Iterate line by line over the data. */
  while (getline(cur_file, line)) {
    linestream = stringstream(line);
    protein_id = -1;

    /* Iterate over each value of a line. */
    while (getline(linestream, word, ',')) {
      if (protein_id == -1)
        protein_id = stoi(word);
      else
        protein_seq = word;
    }

    /* Store protein data in vector, if not folded before. */
    if (finished_ids.count(protein_id) == 0) {
      proteins.push_back(make_pair(protein_id, protein_seq));
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
  proteins = vector<pair<int, string>>(proteins.begin() + begin,
                                       proteins.begin() + begin + num_proteins);

  return proteins;
}

/* Executes a serial session of depth_first_bnb searches over the dataset. */
void dfs_bnb_parallel(string results_path, vector<pair<int, string>> proteins,
                      int dim, string setup, int rank) {
  vector<int> min_conf;
  Protein *p;
  ofstream output_file;
  string results_file;

  chrono::time_point<chrono::high_resolution_clock> t1;
  chrono::time_point<chrono::high_resolution_clock> t2;
  chrono::duration<double> duration;

  /* Check if TMP results already exist for this job, append if so. Create
   * new TMP file and write header otherwise.
   */
  results_file = results_path + "/HP_" + setup + "_dfs_bnb_r" +
                 to_string(rank);
  if (fs::exists(results_file + "_TMP.csv")) {
    /* TMP file exists. */
    output_file.open(results_file + "_TMP.csv", ios::app);
  } else {
    /* TMP file does not exist. */
    output_file.open(results_file + "_TMP.csv", ios::out);
    output_file << "protein_id,algorithm,time,score,checked,placed,hash\n"
                << flush;
  }

  /* Fold each Protein and store resulting minimum energy conformation. */
  for (pair<int, string> cur_p : proteins) {
    p = new Protein(cur_p.second, dim, "HP");
    t1 = chrono::high_resolution_clock::now();
    depth_first_bnb(p, "reach_prune");
    t2 = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(t2 - t1);
    min_conf = p->hash_fold();

    /* Write statistics to output file. */
    output_file << cur_p.first << ",dfs_bnb_mpi," << duration.count() << ","
                << p->get_score() << "," << p->get_solutions_checked() << ","
                << p->get_aminos_placed() << "," << flush;

    /* Write conformation to file and close with newline. */
    output_file << "\"[";
    for (int b : vector<int>(min_conf.begin(), min_conf.end() - 1)) {
      output_file << b << ", " << flush;
    }
    output_file << min_conf.back() << "]\"\n" << flush;

    /* Deallocate the protein. */
    delete p;
  }

  /* If finished, close temporary file and rename to be the final file. */
  output_file.close();
  output_file.clear();

  /* Append TMP file to existing results, or simply rename file. */
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
    /* No existing results, so rename TMP file to results file. */
    rename((results_file + "_TMP.csv").c_str(),
           (results_file + ".csv").c_str());
  }
}

/* Merge results of all ranks in one final results file. */
void merge_results(fs::path results_path, int *lengths, int num_lengths,
                   string dataset_extra, int cluster_size) {
  string fin_results_path;
  string rank_results_path;
  ofstream fin_results_file;
  ifstream rank_results_file;
  struct stat buffer;
  string line;

  for (int l = 0; l < num_lengths; l++) {
    /* Open final results file and write header. */
    if (dataset_extra != "NULL") {
      fin_results_path = results_path.u8string() + "/HP_" +
                       to_string(lengths[l]) + "_" + dataset_extra + 
                       "_dfs_bnb.csv";
    } else {
      fin_results_path = results_path.u8string() + "/HP_" +
                       to_string(lengths[l]) + "_dfs_bnb.csv";
    }

    fin_results_file.open(fin_results_path);
    fin_results_file << "protein_id,algorithm,time,score,checked,"
                     << "placed,hash\n"
                     << flush;

    /* Fill final results file with results for each rank. */
    for (int rank = 0; rank < cluster_size; rank++) {
      if (dataset_extra != "NULL") {
        rank_results_path = results_path.u8string() + "/HP_" +
                          to_string(lengths[l]) + "_" + dataset_extra + 
                          "_dfs_bnb_r" + to_string(rank);
      } else {
        rank_results_path = results_path.u8string() + "/HP_" +
                          to_string(lengths[l]) + "_dfs_bnb_r" +
                          to_string(rank);
      }

      /* Check if rank is fully finished. */
      if (stat((rank_results_path + "_TMP.csv").c_str(), &buffer) == 0) {
        rank_results_path += "_TMP.csv";
      } else {
        rank_results_path += ".csv";
      }

      /* Open rank results, skip header, and write to final results. */
      rank_results_file.open(rank_results_path);
      getline(rank_results_file, line);
      fin_results_file << rank_results_file.rdbuf() << flush;
      rank_results_file.close();
      rank_results_file.clear();
    }

    /* Close final results file. */
    fin_results_file.close();
  }
}

/* Entry point for experiments. */
int main(int argc, char *argv[]) {
  /* Experiment must be called with at least 4 arguments. */
  if (argc < 6) {
    exit(-1);
  }

  /* Determine rank and cluster size. */
  int cluster_size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &cluster_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Parse arguments. */
  string job = argv[1];
  fs::path results_path = fs::path(argv[2]);
  string dataset_path = argv[3];
  string dataset_extra = argv[4];
  int dim = atoi(argv[5]);
  int NUM_FIXED_ARGS = 6;
  int num_lengths = argc - NUM_FIXED_ARGS;
  int lengths[num_lengths];

  for (int i = NUM_FIXED_ARGS; i < argc; i++) {
    lengths[i - NUM_FIXED_ARGS] = atoi(argv[i]);
  }

  /* Print general debug info if this is rank 0. */
  if (rank == 0) {
    cout << "Debug Info:\n";
    cout << "\tResults path:         " << results_path << "\n";
    cout << "\tDataset to use:       " << dataset_path << "\n";
    cout << "\tDataset extra:        " << dataset_extra << "\n";
    cout << "\tDimension of fold:    " << dim << "\n";
    cout << "\tLengths to fold:      ";
    for (int l : lengths) {
      cout << l << " ";
    }
    cout << "\n\n";
  }

  /* Append dataset_extra to the results_path, if it is set. */
  if (dataset_extra != "NULL") {
    results_path += "/" + dataset_extra;
  }

  /* Create results directory if it does not exist. */
  if (!fs::is_directory(results_path)) {
    fs::create_directories(results_path);
  }

  /* Notify creation of a rank. */
  string node_debug = "[" + to_string(rank + 1) + "/" +
                      to_string(cluster_size) + "]\tI'm online..\n";
  cerr.write(node_debug.data(), node_debug.size());

  /* Reload proteins for every length. */
  vector<pair<int, string>> proteins;

  /* Run and time experiment with given arguments per specified length. */
  chrono::time_point<chrono::high_resolution_clock> t1;
  chrono::time_point<chrono::high_resolution_clock> t2;
  chrono::duration<double> duration;
  ofstream time_file;

  /* Store the setup containing the length and potentially the dataset_extra. */
  string setup;

  for (int l : lengths) {
    /* Append dataset_extra to the length for storing the setup. */
    if (dataset_extra != "NULL") {
      setup = to_string(l) + "_" + dataset_extra;
    } else {
      setup = to_string(l);
    }

    /* Load the proteins for this rank. */
    proteins = load_proteins(dataset_path, setup, cluster_size, rank, 
                              results_path);

    /* Time execution and the wait for the other tasks. */
    t1 = chrono::high_resolution_clock::now();
    dfs_bnb_parallel(results_path, proteins, dim, setup, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    t2 = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(t2 - t1);

    /* Rank 0 reports runtime for this length. */
    if (rank == 0) {
      cout << "Finished " + setup + " in " << duration.count() << "\n";
      time_file.open(results_path.u8string() + "/time_HP_" + setup);
      time_file << duration.count() << "\n" << flush;
      time_file.close();
    }
  }

  /* Close MPI processes. */
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  node_debug = "[" + to_string(rank + 1) + "/" + to_string(cluster_size) +
               "]\tFinalized MPI..\n";
  cerr.write(node_debug.data(), node_debug.size());

  /* Merge results into one file per protein length. */
  if (rank == 0) {
    merge_results(results_path, lengths, num_lengths, dataset_extra, 
                  cluster_size);
  }

  return 0;
}
