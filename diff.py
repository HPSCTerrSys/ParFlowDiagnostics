import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

file1="serial.out"
file2="parallel.out"
printline=False
thresh = 0.000001
with open(file1, "r") as f:
    first = f.readlines()    
with open(file2, "r") as f:
    second = f.readlines()

diffs = []
for a,b in zip(first, second):
    a,b = a.strip(), b.strip()
    if printline:
        print(a, end='  ')
    for c,d in zip(a.split(), b.split()):
        try:
            e,f = float(c), float(d)
            if e != f:
                if printline:
                    print('other line:', b)
                diffs.append(abs(e-f))
        except:
            pass
diffs = np.asarray(diffs)
plot = sns.scatterplot(data=diffs)
plot.set_yscale("log")
plt.show()
