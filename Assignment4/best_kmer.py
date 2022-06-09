import pandas as pd
import shutil

df = pd.read_csv('/students/2021-2022/master/Firdaws_DSLS/output/kmers.csv', header = None)

kmers = list(df[0].values)
n50 = list(df[1].values)


index = 0 
best_index = 0 
best_n50 = 0
for i in n50:

    if i > best_n50:
        best_index = index
        best_n50 = i
    index += 1
    
best_kmer = kmers[best_index]



original = '/students/2021-2022/master/Firdaws_DSLS/output/Firdaws_{}/contigs.fa'.format(best_kmer)
target = '/homes/fabadmus/programming_3/Programming3/Assignment4/output/contigs.fa'
csv = '/students/2021-2022/master/Firdaws_DSLS/output/kmers.csv'
target2 = '/homes/fabadmus/programming_3/Programming3/Assignment4/output/kmers.csv'

shutil.move(original, target)
shutil.move(csv, target2)