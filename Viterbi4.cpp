#include "utils.h"


vector<vector<int>> viterbi(HMM &hmm, string &x, string &y) {
    int n = x.size();
    int m = y.size();

    vector<double> prev_row(m+1);
    vector<double> curr_row(m+1);

    vector<vector<int>> backtrack(n, vector<int>(m, -1));

    for(int i=0; i<prev_row.size(); i++) {
        prev_row[i] = log(0);
        curr_row[i] = -INFINITY;
    }

    prev_row[0] = log(1);
    curr_row[0] = log(0);

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            vector<double> max_values;
        
            //insert
            vector<double> v_Ix = {
                curr_row[j-1] + log(hmm.A[0][0]),
                -INFINITY,
                prev_row[j-1] + log(hmm.A[2][0])
            };
            int argmax_Ix = get_argmax(v_Ix);
            if(argmax_Ix == 1) cout << "err";
            max_values.push_back(v_Ix[argmax_Ix] + log(hmm.E[0][char_combination_to_idx[(1, '-') + string(1, y[j-1])]]));

            //delete
            vector<double> v_Iy = {
                -INFINITY,
                prev_row[j] + log(hmm.A[1][1]),
                prev_row[j-1] + log(hmm.A[2][1])
            };
            int argmax_Iy = get_argmax(v_Iy);
            max_values.push_back(v_Iy[argmax_Iy] + log(hmm.E[1][char_combination_to_idx[(1, x[i-1]) + string(1, '-')]]));
            
            //M
            vector<double> v_M = {
                curr_row[j-1],  //insert je prethodni
                prev_row[j],    //delete je prethodni
                prev_row[j-1]   //ostanak u M
            };
            int argmax_M = get_argmax(v_M);
            max_values.push_back(v_M[argmax_M] + log(hmm.E[2][char_combination_to_idx[(1, x[j-1]) + string(1, y[j-1])]]));


            vector<int> argmax_values = {argmax_Ix, argmax_Iy, argmax_M};

            //koje stanje stvara maksimum
            int argmax = get_argmax(max_values); 
            curr_row[j] = max_values[argmax];
            backtrack[i-1][j-1] = argmax_values[argmax];
        }
        prev_row = curr_row;
        curr_row[0] = log(0);
    }

    

    

    return backtrack;
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


    auto result = viterbi(hmm, pair.first, pair.second);

    for(auto v: result) {
        for(auto i: v) {
            cout << i << " ";
        }
        cout << endl;
    }

    // for (auto i : result) {
    //     cout << i << " ";
    // }
    

    return 0;
}