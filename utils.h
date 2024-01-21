#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <stdlib.h>
#include <math.h>

#define INSERTION_INDEX 0
#define DELETION_INDEX 1
#define MATCH_MISMATCH_INDEX 2

#define NUM_STATES 3

#define SMOOTH_ALPHA 2
#define SMOOTH_BETA 2

#define STOPPING_THRESHOLD 0.0001

#define A_MATRIX_FILENAME "A_matrix.txt"
#define E_MATRIX_FILENAME "E_matrix.txt"
#define PI_MATRIX_FILENAME "PI_matrix.txt"

namespace fs = std::filesystem;
using namespace std;



string INPUT_FILEPATH = "/home/jan/Documents/bioinf2/HMM-pariwise-alignment/data/train/HIV1_ALL_2021_genome_DNA_reduced_preprocessed.fasta";
string OUTPUT_DIR = "/home/jan/Documents/bioinf2/HMM-pariwise-alignment/data/";


struct Pair {
    string first;
    string second;
};

//mapa koja parove opažanih stanja mapira u indekse za matricu E
vector<char> characters = {'-', 'A', 'C', 'G', 'T'};

vector<string> get_char_combinations(vector<char> characters) {
    vector<string> combinations;
    for (char firstChar : characters) {
        for (char secondChar : characters) {
            string combination = string(1, firstChar) + string(1, secondChar);
            combinations.push_back(combination);
        }
    }
    return combinations;
}

map<string, int> get_char_combination_map(vector<char> characters) {
    map<string, int> char_combination_to_idx;
    vector<string> char_combinations = get_char_combinations(characters);

    for (int i = 0; i < char_combinations.size(); ++i) {
        char_combination_to_idx[char_combinations[i]] = i;
    }

    return char_combination_to_idx;
}

map<string, int> char_combination_to_idx = get_char_combination_map(characters);

vector<Pair> load_pairs(string filename) {
    vector<Pair> pairs;
    fstream file(filename);
    string line;

    while( getline(file, line) ) {
        stringstream ss(line);
        Pair pair;
        getline(ss, pair.first, ' ');
        getline(ss, pair.second, ' ');
        // pairs.push_back(pair);
        pairs.emplace_back(move(pair)); // Ovo je riješilo problem s brisanjem vector<Pair> pairs iz memorije
    }
    return pairs;
}

int get_state_from_observation(char first, char second) {
    if (first == '-')
        return INSERTION_INDEX;
    else if (second == '-')
        return DELETION_INDEX;
    else
        return MATCH_MISMATCH_INDEX;
}

vector<int> get_states_from_pair(Pair pair) {
    vector<int> states;

    for (int i = 0; i < pair.first.size(); i++) {
        states.push_back(get_state_from_observation(pair.first[i], pair.second[i]));
    }

    return states;
}

void write_matrix_to_file(vector<vector<double>> matrix, string filepath) {
    ofstream output_file;
    output_file.open (filepath);

    for (const auto& row : matrix) {
        for (double element : row) {
            output_file << element << ' ';
        }
        output_file << '\n'; 
    }

    output_file.close();
}
void removeNan(vector<vector<double>> &matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[0].size(); j++) {
            if(isnan(matrix[i][j])) matrix[i][j] = 0.0;
        }
    }
}


struct HMM {
    vector<vector<double>> pi;
    vector<vector<double>> A;
    vector<vector<double>> E;
};