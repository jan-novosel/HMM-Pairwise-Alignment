#include "utils.h"

#define MATCH 1
#define MISMATCH -1
#define GAP -2

int get_score(Pair alignment) {
    string x = alignment.first;
    string y = alignment.second;

    int score = 0;
    for(int i=0; i<x.size(); i++) {
        if(x[i] == y[i]) score += MATCH;
        else if(x[i] == '-' || y[i] == '-') score += GAP;
        else score += MISMATCH;
    }
    return score;
}

vector<int> viterbi(HMM &hmm, string &x, string &y) {
    int n = x.size();
    int m = y.size();

    // logaritmirane vjerojatnost
    vector<vector<double>> V_M(n, vector<double>(m, -1.0));
    vector<vector<double>> V_X(n, vector<double>(m, -1.0));
    vector<vector<double>> V_Y(n, vector<double>(m, -1.0));

    vector<vector<int>> V_M_bt(n, vector<int>(m, -1.0));
    vector<vector<int>> V_X_bt(n, vector<int>(m, -1.0));
    vector<vector<int>> V_Y_bt(n, vector<int>(m, -1.0));



    vector<int> backtrack;

    // inicijalizacija
    for (int i = 0; i < n; i++) {
        V_M[i][0] = -INFINITY;
        V_X[i][0] = -INFINITY;
        V_Y[i][0] = -INFINITY;

        V_M_bt[i][0] = DELETION_INDEX;
        V_X_bt[i][0] = DELETION_INDEX;
        V_Y_bt[i][0] = DELETION_INDEX;
    }

    for (int j = 0; j < m; j++) {
        V_M[0][j] = -INFINITY;
        V_X[0][j] = -INFINITY;
        V_Y[0][j] = -INFINITY;

        V_M_bt[0][j] = INSERTION_INDEX;
        V_X_bt[0][j] = INSERTION_INDEX;
        V_Y_bt[0][j] = INSERTION_INDEX;
    }

    // pocetak punjenja matrice
    V_M[0][0] = 0;
    
    // kraj backtrackinga
    V_M_bt[0][0] = -1;
    V_X_bt[0][0] = -1;
    V_Y_bt[0][0] = -1;

    string char_combination;
    vector<double> candidates;
    int argmax;
    double max_value;

    for (int i = 1; i < n; i++) {

        // if(i % 1000 == 0) cout << (float)(i)/n << endl;

        for (int j = 1; j < m; j++) {
            // V_M
            char_combination = string(1, x[i]) + string(1, y[j]);

            candidates = {
                V_X[i-1][j-1],
                V_Y[i-1][j-1],
                V_M[i-1][j-1]
            };
            argmax = get_argmax(candidates);
            max_value = candidates[argmax];

            V_M[i][j] = log(hmm.E[MATCH_MISMATCH_INDEX][char_combination_to_idx[char_combination]]) + max_value;
            V_M_bt[i][j] = argmax;

            // V_X
            char_combination = string(1, '-') + string(1, y[j]);

            candidates = {
                log(hmm.A[INSERTION_INDEX][INSERTION_INDEX]) + V_X[i][j-1],
                -INFINITY,
                log(hmm.A[MATCH_MISMATCH_INDEX][INSERTION_INDEX]) + V_M[i][j-1]
            };
            argmax = get_argmax(candidates);
            max_value = candidates[argmax];

            V_X[i][j] = log(hmm.E[INSERTION_INDEX][char_combination_to_idx[char_combination]]) + max_value;
            V_X_bt[i][j] = argmax;

            // V_Y
            char_combination = string(1, x[i]) + string(1, '-');

            candidates = {
                -INFINITY,
                log(hmm.A[DELETION_INDEX][DELETION_INDEX]) + V_Y[i-1][j],
                log(hmm.A[MATCH_MISMATCH_INDEX][DELETION_INDEX]) + V_M[i-1][j]
            };
            argmax = get_argmax(candidates);
            if(argmax == 0) argmax = 1;
            max_value = candidates[argmax];

            V_Y[i][j] = log(hmm.E[DELETION_INDEX][char_combination_to_idx[char_combination]]) + max_value;
            V_Y_bt[i][j] = argmax;
        }
    }

    int i = n-1, j = m-1;

    vector<double> final_candidates =  { V_X[i][j],V_Y[i][j],V_M[i][j] };

    int termination = get_argmax(final_candidates);


    vector<int> hidden_state_sequence;
    int current_state = termination;

    while (true) {
        hidden_state_sequence.push_back(current_state);

        switch (current_state) {
        case INSERTION_INDEX:
            j -= 1;
            current_state = V_X_bt[i][j];
            break;
        case DELETION_INDEX:
            i -= 1;
            current_state = V_Y_bt[i][j];
            break;
        case MATCH_MISMATCH_INDEX:
            i -= 1;
            j -= 1;
            current_state = V_M_bt[i][j];
            break;
        }

        if (current_state == -1) {
            break;
        }
    }

    reverse(hidden_state_sequence.begin(), hidden_state_sequence.end());

    return hidden_state_sequence;
}

Pair get_alignment_from_states(vector<int> &hidden_state_sequence, string &x, string &y) {
    string x_aligned = "";
    string y_aligned = "";
    int i = 0, j = 0;

    for (int state : hidden_state_sequence) {
        switch (state) {
            case INSERTION_INDEX:
                x_aligned += "-";
                y_aligned += y[j++];
                break;
            case DELETION_INDEX:
                x_aligned += x[i++];
                y_aligned += "-";
                break;
            case MATCH_MISMATCH_INDEX:
                x_aligned += x[i++];
                y_aligned += y[j++];
                break;
        }
    }

    Pair p = {x_aligned, y_aligned};

    return p;
}


int main() {
    struct HMM hmm;
    hmm.pi = load_matrix("data/PI_matrix.txt_BW");
    hmm.A = load_matrix("data/A_matrix.txt_BW");
    hmm.E = load_matrix("data/E_matrix.txt_BW");

    vector<Pair> pairs = load_pairs("data/test/HIV1_ALL_2021_genome_DNA_reduced_preprocessed.fasta");

    vector<int> scores;
    float cnt = 0;
    for(Pair p: pairs) {
        cout << cnt/pairs.size() << endl;
        vector<int> hidden_state_sequence = viterbi(hmm, p.first, p.second);
        Pair alignment = get_alignment_from_states(hidden_state_sequence, p.first, p.second);

        int score = get_score(alignment);
        scores.push_back(score);
        cnt++;
    }

    // Pair pair = {};
    // pair.first = "TCCCAACGAAGACAAGATATCCTTGATCTAAAAA";
    // pair.second = "GGATGGGTTAATTTACTCCCGGAAAAAGAGATCCTTGATCT";


    // auto hidden_state_sequence = viterbi(hmm, pair.first, pair.second);

    // // for(auto v: result) {
    // //     for(auto i: v) {
    // //         cout << i << " ";
    // //     }
    // //     cout << endl;
    // // }

    // for (auto i : hidden_state_sequence) {
    //     cout << i << " ";
    // }
    // cout << endl;

    // Pair alignment = get_alignment_from_states(hidden_state_sequence, pair.first, pair.second);
    // int score = get_score(alignment);

    // vector<int> scores;
    // scores.push_back(score);
    // scores.push_back(score);
    

    ofstream output_file;
    output_file.open("data/results.txt", ios::app);

    for(int i: scores) {
        output_file << i << ' ';
    }

    output_file << '\n';
    output_file.close();
    

    return 0;
}