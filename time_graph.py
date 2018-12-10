import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

with open('./test_ref/time_compare.txt') as file:
	lines = file.readlines()
	print(lines)
	cache = []
	normal = []
	for line in lines:
		cache.append(float(line.strip().split(' ')[0]))
		normal.append(float(line.strip().split(' ')[1]))
	cache = np.cumsum(np.array(cache))
	normal = np.cumsum(np.array(normal))
	df = pd.DataFrame({"number of ref": [int(i) for i in range(1,len(cache)+1)],
                   "cache alignment" : cache,
                   "ksw2" : normal })
df.plot(x="number of ref", y=["cache alignment", "ksw2"])
plt.show()