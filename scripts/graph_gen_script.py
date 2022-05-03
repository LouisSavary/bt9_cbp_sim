from cProfile import label
from statistics import mean
import matplotlib.pyplot as plt
import numpy as np
import asyncio
import sys


NUM_COL = 10
NUM_ROW = 42

f = open(sys.argv[1], "r")

f.readline()
f.readline()

headers = f.readline().split()
print(headers)

f.readline()

results = [[] for _ in range(NUM_ROW)]
for i in range(NUM_ROW):
  results[i] =  f.readline().split()
  # results[i,1:] = list(map(float, results[i,1:]))
f.close()
for row in results:
  row[1:] = [float(i) for i in row[1:]]

results_sorted = sorted(results, key=(lambda row: row[3]))
# usage in fuction of benchmark
graph = plt.bar(x=np.arange(NUM_ROW), height=([row[3] for row in results_sorted]))
plt.xticks(np.arange(NUM_ROW), labels=[row[0] for row in results_sorted],rotation = 90)
plt.hlines(mean([row[3] for row in results]),0,NUM_ROW,color='r',label="mean")
plt.title("Proportion of instruction executed in a predicted trace by workbench")
plt.tight_layout()
plt.show()
# print(data)


