#include "utils.h"

int get_argmax(vector<double> v) {
    auto it = std::max_element(v.begin(), v.end());
    if(it != v.end())
    {
        return distance(v.begin(), it);
    }
    return -1;
}

vector<int> Viterbi(HMM &hmm, Pair pair) {
    vector<vector<double>> max_p(pair.first.size(), vector<double>(NUM_STATES, 1.0));
    vector<vector<double>> previous_states(pair.first.size()-1, vector<double>(NUM_STATES, 1.0));
    //vjerojatnost prvog stanja
    string char_combination = string(1, pair.first[0]) + string(1, pair.second[0]);
    for(int i = 0; i<NUM_STATES; i++) {
        max_p[0][i] *= hmm.pi[0][i] * hmm.E[i][char_combination_to_idx[char_combination]];

    }
    //vjerojatnosti ostalih stanja
    //unaprijedni prolaz
    for(int t = 1; t<pair.first.size(); t++) {
        char_combination = string(1, pair.first[t]) + string(1, pair.second[t]);
        for(int i = 0; i<NUM_STATES; i++) {
            vector<double> p(NUM_STATES, 1.0);
            for(int j = 0; j<NUM_STATES; j++) {
                p[j] *= max_p[t-1][j] * hmm.A[i][j] * hmm.E[j][char_combination_to_idx[char_combination]];
            }
            int argmax = get_argmax(p);
            max_p[t][i] = p[argmax];
            previous_states[t-1][i] = argmax;
        }
    }
    //unazadni prolaz za dobivanje sekvence
    vector<int> bestSeq;

    int argmax = get_argmax(max_p[pair.first.size()-1]);
    bestSeq.insert(bestSeq.begin(), argmax);
    int prev_state = previous_states[pair.first.size()-2][argmax];
    for(int t = pair.first.size()-2; t>=0; t--) {
        prev_state = previous_states[t][prev_state];
        bestSeq.insert(bestSeq.begin(), prev_state);
    }
    return bestSeq;
}

int main() {
    struct HMM hmm;
    hmm.pi = load_matrix("data/train/PI_matrix.txt");
    hmm.A = load_matrix("data/train/A_matrix.txt");
    hmm.E = load_matrix("data/train/E_matrix.txt");

    // vector<Pair> pair = load_pairs(TEST_FILE);
    // vector<int> seq = Viterbi(hmm, pair[0]);

    Pair pair = {};
    pair.first = "TGGAAGGGCTAATTCACTCCCAACGAAGACAA";
    pair.second = "GGA-TGGGTTAATTTACTCCCGGAAAAGACAA";

    vector<int> seq = Viterbi(hmm, pair);

    for(int i: seq) {
        cout << i << " ";
    }

    return 0;
}
