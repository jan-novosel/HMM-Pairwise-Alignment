import numpy as np
import matplotlib.pyplot as plt


with open('data/results.txt', 'r') as file:
    line = file.readline()
    viterbi = list(map(int, line.split()))

with open('data/results_nw.txt', 'r') as file:
    line = file.readline()
    nw = list(map(int, line.split()))

plt.plot(nw, label='Needleman-Wunsch', color='limegreen')
plt.plot(viterbi, label='Viterbi', color='pink')
plt.legend()
plt.show()