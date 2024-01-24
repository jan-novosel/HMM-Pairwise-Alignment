#include "utils.h"

vector<vector<double>> viterbi(HMM &hmm, string &x, string &y) {
    int n = x.size();
    int m = y.size();
    vector<vector<double>> V_M(n, vector<double>(m, -1.0));
    vector<vector<double>> V_X(n, vector<double>(m, -1.0));
    vector<vector<double>> V_Y(n, vector<double>(m, -1.0));

    // inicijalizacija
    for (int i = 0; i < n; i++) {
        V_M[i][0] = 0;
        V_X[i][0] = 0;
        V_Y[i][0] = 0;
    }

    for (int j = 0; j < m; j++) {
        V_M[0][j] = 0;
        V_X[0][j] = 0;
        V_Y[0][j] = 0;
    }

    V_M[0][0] = 1;

    string char_combination;
    vector<double> candidates;
    int argmax;
    double max_value;

    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            // V_M
            char_combination = string(1, x[i]) + string(1, y[j]);

            candidates = {
                hmm.A[2][2] * V_M[i-1][j-1],
                hmm.A[0][2] * V_X[i-1][j-1],
                hmm.A[1][2] * V_Y[i-1][j-1]
            };
            argmax = get_argmax(candidates);
            max_value = candidates[argmax];

            V_M[i][j] = hmm.E[2][char_combination_to_idx[char_combination]] * max_value;

            // V_X
            char_combination = string(1, '-') + string(1, y[j]);

            candidates = {
                hmm.A[2][0] * V_M[i-1][j],
                hmm.A[0][0] * V_X[i-1][j]
            };
            argmax = get_argmax(candidates);
            max_value = candidates[argmax];

            V_M[i][j] = hmm.E[0][char_combination_to_idx[char_combination]] * max_value;

            // V_Y
            char_combination = string(1, x[i]) + string(1, '-');

            candidates = {
                hmm.A[2][1] * V_M[i-1][j],
                hmm.A[1][1] * V_Y[i-1][j]
            };
            argmax = get_argmax(candidates);
            max_value = candidates[argmax];

            V_M[i][j] = hmm.E[1][char_combination_to_idx[char_combination]] * max_value;
        }
    }
    return V_M;
}

int main() {
    struct HMM hmm;
    hmm.pi = load_matrix("data/PI_matrix.txt_BW");
    hmm.A = load_matrix("data/A_matrix.txt_BW");
    hmm.E = load_matrix("data/E_matrix.txt_BW");

    // vector<Pair> pair = load_pairs(TEST_FILE);
    // vector<int> seq = Viterbi(hmm, pair[0]);

    Pair pair = {};
    pair.first = "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCT";
    pair.second = "GGATGGGTTAATTTACTCCCGGAAAAGACAAGAGATCCTTGATCT";


    vector<vector<double>> result = viterbi(hmm, pair.first, pair.second);

    for(vector<double> v: result) {
        for(double i: v) {
            cout << i << " ";
        }
        cout << endl;
    }

    return 0;
}