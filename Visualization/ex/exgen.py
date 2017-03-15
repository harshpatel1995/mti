import csv
import numpy as np

for i in range(20):
    with open(str(i) + '.gra', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([210007, 221988, 235279, 163164, 273119, 469008, 1335916, 1385755, 1435461, 511145])
        x = np.random.gumbel(0.5, 0.3, 10)
        x /= sum(x)
        if i % 3 == 0:
            x.sort()
        writer.writerow(x)
        writer.writerow(x * 0.342)
