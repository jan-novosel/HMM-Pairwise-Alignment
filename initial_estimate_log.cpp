#include "utils.h"


void save_initial_estimate(vector<Pair> pairs, string dir) {
    string a_matrix_filename(A_MATRIX_FILENAME);
    string e_matrix_filename(E_MATRIX_FILENAME);
    string pi_matrix_filename(PI_MATRIX_FILENAME);

    fs::path a_matrix_filepath = dir;
    a_matrix_filepath /= a_matrix_filename;
    fs::path e_matrix_filepath = dir;
    e_matrix_filepath /= e_matrix_filename;
    fs::path pi_matrix_filepath = dir;
    pi_matrix_filepath /= pi_matrix_filename;

    double num_pairs = pairs.size();
    int num_char_combinations = char_combination_to_idx.size();

    vector<int> state_freqs(NUM_STATES, 0);
    vector<int> last_state_counter(NUM_STATES, 0);

    vector<double> pi_matrix(NUM_STATES, 0.0f);
    vector<vector<double>> a_matrix(NUM_STATES, vector<double>(NUM_STATES, 0.0f));
    vector<vector<double>> e_matrix(NUM_STATES, vector<double>(num_char_combinations, 0.0f));

    int progress = 0;

    for (Pair pair : pairs) {
        vector<int> states = get_states_from_pair(pair);

        // pi matrix
        pi_matrix[states[0]]++;

        for (int i = 0; i < states.size(); i++) {
            int state = states[i];
            state_freqs[state]++;

            // a_matrix
            if (i != states.size() - 1) {
                int following_state = states[i + 1];
                a_matrix[state][following_state]++;
            } else {
                last_state_counter[i]++;
            }

            // e_matrix
            string char_combination = string(1, pair.first[i]) + string(1, pair.second[i]);
            int combination_idx = char_combination_to_idx[char_combination];
            e_matrix[state][combination_idx]++;
        }

        // if (progress % 10 == 0) {
        //     cout << progress * 100 / num_pairs << "%" << endl;
        // }
        // ++progress;
    }

    // normalize pi_matrix
    for (int i = 0; i < pi_matrix.size(); i++) {
        pi_matrix[i] += (SMOOTH_ALPHA - 1);
        pi_matrix[i] /= (num_pairs + SMOOTH_ALPHA + SMOOTH_BETA - 2);
    }

    // normalize a_matrix, e_matrix
    for (int i = 0; i < a_matrix.size(); i++) {
        for (int j = 0; j < a_matrix[i].size(); j++) {
            a_matrix[i][j] += (SMOOTH_ALPHA - 1);
            a_matrix[i][j] /= (state_freqs[i] - last_state_counter[i] + SMOOTH_ALPHA + SMOOTH_BETA - 2);
        }
        for (int j = 0; j < e_matrix[i].size(); j++) {
            e_matrix[i][j] += (SMOOTH_ALPHA - 1);
            e_matrix[i][j] /= (state_freqs[i] + SMOOTH_ALPHA + SMOOTH_BETA - 2);
        }
    }

    // for(vector<double> v: a_matrix) {
    //     for(double i: v) {
    //         cout << i << " ";
    //     }
    //     cout << endl;
    // }

    //logaritmiraj matrice
    for (int i = 0; i < pi_matrix.size(); i++) {
        pi_matrix[i] = log(pi_matrix[i]);
    }

    for (int i = 0; i < a_matrix.size(); i++) {
        for (int j = 0; j < a_matrix[i].size(); j++) {
            a_matrix[i][j] = log(a_matrix[i][j]);
        }
        for (int j = 0; j < e_matrix[i].size(); j++) {
            e_matrix[i][j] = log(e_matrix[i][j]);
        }
    }

    vector<vector<double>> pi_matrix_2d = {pi_matrix};

    write_matrix_to_file(pi_matrix_2d, pi_matrix_filepath);
    write_matrix_to_file(a_matrix, a_matrix_filepath);
    write_matrix_to_file(e_matrix, e_matrix_filepath);

}


int main() {
    vector<Pair> pairs = load_pairs(INPUT_FILEPATH);
    save_initial_estimate(pairs, OUTPUT_DIR);
    
    return 0;
}
