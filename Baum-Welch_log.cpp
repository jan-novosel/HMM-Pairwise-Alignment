#include "utils.h"


vector<vector<double>> forward(struct HMM &hmm, Pair pair) {
    int T = pair.first.size();
    vector<vector<double>> alpha(T, vector<double>(NUM_STATES, 0.0));
    for(int t=0; t<T; t++) {
        string char_combination = string(1, pair.first[t]) + string(1, pair.second[t]);
        for(int i=0; i<NUM_STATES; i++) {    //za svako od skrivenih stanja
            if(t == 0) {
                alpha[t][i] = hmm.pi[0][i] * hmm.E[i][char_combination_to_idx[char_combination]];
            } else {
                for(int j=0; j<NUM_STATES; j++) {
                    alpha[t][i] += exp(log(alpha[t-1][j]) + log(hmm.A[j][i]) + 
                        log(hmm.E[i][char_combination_to_idx[char_combination]]));
                    //if(alpha[t][i] == 0) alpha[t][i] = 0.0000000001;
                }
            }
        }
    }
    return alpha;
}

vector<vector<double>> backward(struct HMM &hmm, Pair pair) {
    int T = pair.first.size();
    vector<vector<double>> betha(T, vector<double>(NUM_STATES, 0.0));
    for(int t=T-1; t>=0; t--) {
        string char_combination;
        if(t != T-1) {
            char_combination = string(1, pair.first[t+1]) + string(1, pair.second[t+1]);
        }
        for(int i=0; i<NUM_STATES; i++) {
            if(t == T-1) {
                betha[t][i] = 1.0;
            } else {
                for(int j=0; j<NUM_STATES; j++) {
                    betha[t][i] += exp(log(betha[t+1][j]) + log(hmm.A[i][j]) + 
                        log(hmm.E[j][char_combination_to_idx[char_combination]]));
                    //if(betha[t][i] == 0) betha[t][i] = 0.0000000001;
                }
            }
        }
    }
    return betha;
}

vector<vector<double>> computeGamma(HMM &hmm, vector<vector<double>> &alpha, vector<vector<double>> &betha) {
    int T = alpha.size();
    vector<vector<double>> gamma(T, vector<double>(NUM_STATES, 0.0));

    for (int t = 0; t < T; t++) {
        double denominator = 0.0;
        for (int i = 0; i < NUM_STATES; i++) {
            gamma[t][i] = alpha[t][i] * betha[t][i];
            denominator += gamma[t][i];
        }
        for  (int i = 0; i < NUM_STATES; i++) {
            gamma[t][i] /= denominator;
        }
    }
    return gamma;
}


vector<vector<vector<double>>> computeXi(HMM &hmm, vector<vector<double>> &alpha, vector<vector<double>> &betha, Pair &pair) {
    int T = alpha.size();
    vector<vector<vector<double>>> xi(T-1, vector<vector<double>>(NUM_STATES, vector<double>(NUM_STATES, 0.0)));

    for (int t = 0; t < T-1; t++) {
        string char_combination = string(1, pair.first[t+1]) + string(1, pair.second[t+1]);
        double denominator = 0.0;

        for (int i = 0; i < NUM_STATES; i++) {
            for (int j = 0; j < NUM_STATES; j++) {
                xi[t][i][j] = exp(log(alpha[t][i]) + log(hmm.A[i][j] + betha[t+1][j]) +
                    log(hmm.E[j][char_combination_to_idx[char_combination]]));
				denominator += xi[t][i][j];
            }
        }

        for (int i = 0; i < NUM_STATES; i++) {
            for (int j = 0; j<NUM_STATES; j++) {
                xi[t][i][j] /= denominator;
            }
        }
    }
    return xi;
}

double matrix_1_norm(vector<vector<double>> matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    double norm = 0.0;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            norm += abs(matrix[i][j]);
        }
    }

    return norm;
}

bool stopping_criterion(HMM &hmm, HMM &prev_hmm, double threshold) {
    double hmm_1_norm = matrix_1_norm(hmm.pi) + matrix_1_norm(hmm.A) + matrix_1_norm(hmm.E);
    double prev_hmm_1_norm = matrix_1_norm(prev_hmm.pi) + matrix_1_norm(prev_hmm.A) + matrix_1_norm(prev_hmm.E);

    return abs(hmm_1_norm - prev_hmm_1_norm) < threshold;
}

void Baum_Welch(HMM &hmm, int max_iter, vector<Pair> &pairs) {
    int T;
    int N = pairs.size();

    HMM prev_hmm = hmm;

    for (int iter = 0; iter < max_iter; iter++) {
        
        vector<vector<vector<double>>> gammas;
        vector<vector<vector<vector<double>>>> xis;
        
        float cnt = 0.0f;

        for(Pair pair: pairs) {
            vector<vector<double>> alpha = forward(hmm, pair);
            vector<vector<double>> betha = backward(hmm, pair);

            gammas.push_back(computeGamma(hmm, alpha, betha));
            xis.push_back(computeXi(hmm, alpha, betha, pair));
        }

        //update pi
        for (int i = 0; i < NUM_STATES; i++) {
            double pi = 0.0;
            for (int n = 0; n < N; n++) {
                pi += gammas[n][0][i];
            }
            pi /= N;
            hmm.pi[0][i] = pi;
        }

        //update A
        for (int i = 0; i < NUM_STATES; i++) {
            double gamma_sum = 0.0;

            for (int n = 0; n < N; n++) {
                T = pairs[n].first.size();
                for (int t = 0; t < T; t++) {
                    gamma_sum += gammas[n][t][i];
                }
            }

            for (int j = 0; j < NUM_STATES; j++) {
                double xi_sum = 0.0;

                for (int n = 0; n < N; n++) {
                    T = pairs[n].first.size();
                    for (int t = 0; t < T-1; t++) {
                        xi_sum += xis[n][t][i][j];
                    }
                }
                hmm.A[i][j] = xi_sum / gamma_sum;
            }   
        }

        for (int i = 0; i < NUM_STATES - 1; i++) {
            for (int j = 0; j < NUM_STATES; j++) {
                if (i == j) {
                    hmm.A[i][j] -= 0.2;
                }
                else {
                    hmm.A[i][j] += 0.1;
                }
            }
        }

        //update E
        for (int i = 0; i < NUM_STATES; i++) {
            double denominator = 0.0;

            for (int n = 0; n < N; n++) {
                T = pairs[n].first.size();
                for (int t = 0; t < T; t++) {
                    denominator += gammas[n][t][i];
                }
            }

            for (int k = 0; k < hmm.E[0].size(); k++) {
                double numerator = 0.0;
                
                for (int n = 0; n < N; n++) {
                    for (int t = 0; t < T; t++) {
                        T = pairs[n].first.size();
                        string char_combination = string(1, pairs[n].first[t]) + string(1, pairs[n].second[t]);
                        if (char_combination_to_idx[char_combination] == k) {
                            numerator += gammas[n][t][i];
                        }
                    }
                }
                hmm.E[i][k] = numerator / denominator;
            }
        }

        if (stopping_criterion(hmm, prev_hmm, STOPPING_THRESHOLD)) {
            cout << "criterion " << iter << endl;
            break;
        }

    prev_hmm = hmm;

    }

}

void save_baum_welch_estimate(HMM &hmm, string dir) {
    string a_matrix_filename(A_MATRIX_FILENAME + string("_BW"));
    string e_matrix_filename(E_MATRIX_FILENAME + string("_BW"));
    string pi_matrix_filename(PI_MATRIX_FILENAME + string("_BW"));

    fs::path a_matrix_filepath = dir;
    a_matrix_filepath /= a_matrix_filename;
    fs::path e_matrix_filepath = dir;
    e_matrix_filepath /= e_matrix_filename;
    fs::path pi_matrix_filepath = dir;
    pi_matrix_filepath /= pi_matrix_filename;

    ofstream output_file;

    write_matrix_to_file(hmm.pi, pi_matrix_filepath);
    write_matrix_to_file(hmm.A, a_matrix_filepath);
    write_matrix_to_file(hmm.E, e_matrix_filepath);
}

int main() {
    struct HMM hmm;
    hmm.pi = load_matrix("data/PI_matrix.txt");
    hmm.A = load_matrix("data/A_matrix.txt");
    hmm.E = load_matrix("data/E_matrix.txt");


    // for(vector<double> v: hmm.pi) {
    //     for(double i: v) {
    //         cout << i << " ";
    //     }
    //     cout << endl;
    // }
    // for(vector<double> v: hmm.A) {
    //     for(double i: v) {
    //         cout << i << " ";
    //     }
    //     cout << endl;
    // }
    // for(vector<double> v: hmm.E) {
    //     for(double i: v) {
    //         cout << i << " ";
    //     }
    //     cout << endl;
    // }
    

    vector<Pair> pairs = load_pairs(INPUT_FILEPATH);
    Baum_Welch(hmm, 2, pairs);
    save_baum_welch_estimate(hmm, OUTPUT_DIR);

    return 0;
}