''' This is the main script that runs the data preparationing and traing scripts
usage
python assignment6.py -i /data/dataprocessing/interproscan/all_bacilli.tsv -o /homes/fabadmus/programming_3/Programming3/Assignment6/results
'''
import os
import warnings
import argparse as ap
import pandas as pd
from dask.distributed import Client
from data_cleaner import PrepareData
from classifier import Train

warnings.filterwarnings('ignore')

argparser = ap.ArgumentParser(
                            description= "Script that Prepares and classifies funtions of protein")
argparser.add_argument("--IN_FILE_PATH", "-i", action="store",  type = str,
                            help="TSV file of InterPROscan annotations")
argparser.add_argument("--OUT_FILE_PATH", "-o", action="store",  type = str,
                            help="Path to save results")
args = argparser.parse_args()
print("Getting data: ", args.IN_FILE_PATH)

in_path = args.IN_FILE_PATH

result_dir = args.OUT_FILE_PATH
if not os.path.isdir(result_dir):
    os.mkdir(result_dir) 
print("Result directory created")

# '/data/dataprocessing/interproscan/all_bacilli.tsv'
'''
Prepare data
'''
MODELDATA_PATH = f'{result_dir}/model_data'
my_data = PrepareData(in_path, MODELDATA_PATH)
loaded_data = my_data.load_data()
filtered_data = my_data.filter_data(loaded_data)
labeled_data = my_data.label_size(filtered_data)
clean_data = my_data.clean_data(labeled_data)
my_data.get_model_df(clean_data)
'''
Do machine learning
'''
data = pd.read_csv(MODELDATA_PATH)
REPORT_PATH = f'{result_dir}/model_report'
MODEL_PATH = f'{result_dir}/randomforest_model.sav'
model = Train(data, REPORT_PATH , MODEL_PATH)
# train_features, test_features, train_labels, test_labels = model.train_prep(data)
client = Client(processes=False, threads_per_worker=4,
            n_workers=1, memory_limit='128GB')
model.train_model()
