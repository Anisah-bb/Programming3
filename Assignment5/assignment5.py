import pyspark
from pyspark.sql import SQLContext
from pyspark import SparkContext
import csv
import pyspark.sql.functions as func
from pyspark.sql.functions import desc, avg, split, explode
from pyspark.sql.types import StructType
from pyspark.sql.types import StructField
from pyspark.sql.types import StringType, IntegerType
import pandas as pd

sc=SparkContext('local[16]')
# define schema
customSchema = StructType([
    StructField("Protein_accession", StringType(), True),        
    StructField("MD5", StringType(), True),
    StructField("Seq_len", IntegerType(), True),
    StructField("Analysis", StringType(), True),        
    StructField("Signature_accession", StringType(), True),
    StructField("Signature_description", StringType(), True),
    StructField("Start", IntegerType(), True),        
    StructField("Stop", IntegerType(), True),
    StructField("Score", StringType(), True),
    StructField("Status", StringType(), True),        
    StructField("Date", StringType(), True),
    StructField("InterPro_accession", StringType(), True),
    StructField("InterPro_discription", StringType(), True),
    StructField("GO_annotations", StringType(), True),        
    StructField("Pathways", StringType(), True),
])

#first read from csv file and create a dataframe with created Schema
df = SQLContext(sc).read.csv('/data/dataprocessing/interproscan/all_bacilli.tsv', sep=r'\t',
                            schema=customSchema,
                            inferSchema= True)

# answer 1
a = df.select("InterPro_accession").distinct()
answer11 = a.count()
answer12 =  a._sc._jvm.PythonSQLUtils.explainString(a._jdf.queryExecution(), 'simple')

# answer 2
df = df.filter(df.InterPro_accession!="-")
df = df.filter(df.Protein_accession!="-")
c_count = df.groupBy("Protein_accession").agg(func.count("InterPro_accession"))
c_avg = c_count.agg(func.mean("count(InterPro_accession)"))
answer21 = c_avg.collect()[0][0]
answer22 = c_avg._sc._jvm.PythonSQLUtils.explainString(c_avg._jdf.queryExecution(), 'simple')
# answer 3
b = df.withColumn("go", explode(split(df.GO_annotations, "[|]")))
b = b.filter(b.go!="-")
b = b.groupBy("go").agg(func.count("go"))
b = b.sort(desc("count(go)"))
answer31 = b.collect()[0][0]
answer32 = b._sc._jvm.PythonSQLUtils.explainString(b._jdf.queryExecution(), 'simple')

# answer 4
c1 = df.withColumn('size', (df['Stop'] - df['Start']))
c2 = c1.agg(avg('size'))
answer41 = c2.collect()[0][0]
answer42 = c1._sc._jvm.PythonSQLUtils.explainString(c1._jdf.queryExecution(), 'simple')

# answer 5
d = df.filter(df.InterPro_accession!="-").groupBy("InterPro_accession").count()
d1 = d.sort(func.desc("count"))
d2 = d1.select("InterPro_accession").head(10)
answer51=[i[0] for i in d2]
answer52 = d1._sc._jvm.PythonSQLUtils.explainString(d1._jdf.queryExecution(), 'simple')

# answer 6
e1 = c1.where (((c1.size / c1.Seq_len) * 100) >= 90)
e1 = e1.sort(desc("size"))
e2 = e1.select("InterPro_accession").head(200)
e3=[i[0] for i in e2]
answer61 = []
for i in e3:
    if not i in answer61:
        answer61.append(i)
answer61 = answer61[:10]
answer62 = e1._sc._jvm.PythonSQLUtils.explainString(d1._jdf.queryExecution(), 'simple')
# answer 7
g = df.filter(df.InterPro_discription!="-")
g = g.withColumn("discription_word", explode(split(df.InterPro_discription, '\s|,\s')))
g = g.groupBy("discription_word").agg(func.count("discription_word"))
g = g.sort(desc("count(discription_word)"))
g1 = g.head(10)
answer71=[i[0] for i in g1]
answer72 = g._sc._jvm.PythonSQLUtils.explainString(g._jdf.queryExecution(), 'simple')

# answer 8
g1 = g.tail(10)
answer81=[i[0] for i in g1]
answer82 = answer72

# answer 9
h  = c1.where (((c1.size / c1.Seq_len) * 100) >= 90)
k = h.filter(df.InterPro_discription!="-")
k = k.withColumn("discription_word", explode(split(df.InterPro_discription, '\s|,\s')))
k = k.groupBy("discription_word").agg(func.count("discription_word"))
k = k.sort(desc("count(discription_word)"))
k1 = k.head(10)
answer91=[i[0] for i in k1]
answer92 = k._sc._jvm.PythonSQLUtils.explainString(k._jdf.queryExecution(), 'simple')

# answer 10
m = df.filter(df.InterPro_accession !="-").groupBy('Protein_accession','Seq_len')
m = m.agg(func.count("Seq_len"))
answer101 = m.corr('Seq_len','count(Seq_len)')**2
answer102 = m._sc._jvm.PythonSQLUtils.explainString(m._jdf.queryExecution(), 'simple')

output = {'number' : range(1,11),
 'answer1' : [answer11, answer21, answer31, answer41, answer51, answer61, answer71, answer81, answer91, answer101],
 'answer2' : [answer12, answer22, answer32, answer42, answer52, answer62, answer72, answer82, answer92, answer102],
}
output = pd.DataFrame(output)
output.to_csv('output.csv', sep=',', index=False)