import matplotlib.pyplot as plt
import json
import sys

peptideCount2TaxonCountJsonFile = sys.argv[1]
peptideCount2TaxonCountChart = sys.argv[2]

ranks = ['forma', 'varietas', 'subspecies', 'species', 'species_subgroup', 'species_group', 'subgenus', 'genus',
         'subtribe',
         'tribe', 'subfamily', 'family', 'superfamily', 'parvorder', 'infraorder', 'suborder', 'order', 'superorder',
         'infraclass', 'subclass', 'class', 'superclass', 'subphylum', 'phylum', 'superphylum', 'subkingdom', 'kingdom',
         'superkingdom']

file = open(peptideCount2TaxonCountJsonFile, 'r', encoding='utf-8')
# rank => {x-unique peptide => taxon count}
data = json.load(file)

# plot
ranksLen = len(ranks)
cols = 3
rows = -1
if ranksLen % cols == 0:
    rows = ranksLen // cols
else:
    rows = ranksLen // cols + 1


def sortFun(key):
    return int(key)


plt.figure(figsize=(25, 80))
for i in range(ranksLen):
    plt.subplot(rows, cols, i + 1)
    peptideCount2TaxonCount = data[ranks[i]]
    peptideCountList = sorted(peptideCount2TaxonCount.keys(), key=sortFun)
    if len(peptideCountList) > 0:
        taxonCountList = []
        for p in peptideCountList:
            taxonCountList.append(peptideCount2TaxonCount[p])
        plt.bar(peptideCountList, taxonCountList, color="#75b79e")
    else:
        plt.xticks([])
        plt.yticks([])
    plt.xlabel("Peptide count", fontweight='bold')
    plt.ylabel("Taxon Count", fontweight='bold')
    plt.title("Peptide distribution on " + ranks[i], fontweight='bold', fontsize='large', position=[0.5, 1.05])
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.1)

plt.savefig(peptideCount2TaxonCountChart)
# plt.show()
