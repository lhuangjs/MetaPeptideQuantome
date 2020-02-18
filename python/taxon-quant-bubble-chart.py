import matplotlib.pyplot as plt
import json
import sys

taxonQuantJsonFile = sys.argv[1]
taxonQuantChartFile = sys.argv[2]

## parse json ##
file = open(taxonQuantJsonFile, 'r', encoding='utf-8')
taxonQuant = json.load(file)

# sample
x = []
# taxon
y = []
# quant
z = []
# color
color = []

# taxon
taxa = list(taxonQuant['taxon2Quants'].keys())
for i in range(len(taxa)):
    quantVals = taxonQuant['taxon2Quants'][taxa[i]]
    # sample
    for j in range(len(taxonQuant['samples'])):
        x.append(j)
        y.append(i)
        if quantVals[j] == "-Infinity":
            z.append(0)
        else:
            z.append(quantVals[j])

# use the scatter function
fig = plt.figure(figsize=(6, 15))
plt.scatter(x, y, z, c=x)
plt.xticks(list(range(len(taxonQuant['samples']))), taxonQuant['samples'])
plt.yticks(list(range(len(taxa))), taxa)
plt.xticks(fontsize=12, rotation=90)
plt.yticks(fontsize=12)
plt.title("Taxa distribution on " + taxonQuant['rank'] + ' level', fontweight='bold', fontsize='large',
          position=[0.5, 1.05])
fig.tight_layout()
plt.savefig(taxonQuantChartFile)
