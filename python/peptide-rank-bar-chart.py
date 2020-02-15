import matplotlib.pyplot as plt
import json
import sys

peptide2RankJsonFile = sys.argv[1]
peptide2RankChart = sys.argv[2]

ranks = ['forma', 'varietas', 'subspecies', 'species', 'species_subgroup', 'species_group', 'subgenus', 'genus',
         'subtribe',
         'tribe', 'subfamily', 'family', 'superfamily', 'parvorder', 'infraorder', 'suborder', 'order', 'superorder',
         'infraclass', 'subclass', 'class', 'superclass', 'subphylum', 'phylum', 'superphylum', 'subkingdom', 'kingdom',
         'superkingdom']

file = open(peptide2RankJsonFile, 'r', encoding='utf-8')
data = json.load(file)
peptideCountForSubrankList = []
peptideCountForRankList = []

for rank in ranks:
    peptideCountForSubrankList.append(data[rank]['peptideCountForSubrank'])
    peptideCountForRankList.append(data[rank]['peptideCountForRank'])

# plot
plt.figure(figsize=(15, 10))
plt.barh(ranks, peptideCountForSubrankList, color="#75b79e")
plt.barh(ranks, peptideCountForRankList, left=peptideCountForSubrankList, color="#fb8d62")

plt.xlabel("Peptide count", fontweight='bold')
plt.ylabel("Rank", fontweight='bold')
plt.title("Peptide distribution on different rank", fontweight='bold', fontsize ='large', position=[0.5, 1.05])
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.1)
plt.legend(['peptide count specified to this rank', 'peptide count specified to subranks'],
           loc ='lower right')
plt.savefig(peptide2RankChart)
# plt.show()


