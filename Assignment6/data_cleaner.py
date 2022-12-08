''' This script prepares the data for machine learning'''
import pandas as pd
import dask.array as da
import dask.dataframe as dd
from sklearn.preprocessing import LabelEncoder

class PrepareData():
    '''
    class to perform data wrangling
    '''
    def __init__(self, in_path, modeldata_path):
        self.inpath = in_path
        self.outpath = modeldata_path
    def load_data(self):
        '''
        function to load data into a datafram from a given path
        returns: data as a dataframe
        '''
        headers = ["Protein_accession", "MD5", "Seq_len", "Analysis",
                   "Signature_accession","second_layer"
                    "Signature_description", "Start",
                    "Stop", "Score", "Status", "Date",
                    "InterPro_accession", "InterPro_discription",
                    "GO_annotations", "Pathways"]
        protein_data = dd.read_csv(self.inpath, sep=r'\t', engine='python', names=headers,
                        dtype={"Protein_accession": str,
                                "MD5": str,
                                "Seq_len": int,
                                "Analysis": str,
                                "Signature_accession": str,
                                "Signature_description": str,
                                "Start": int,
                                "Stop": int,
                                "Score": str,
                                "Status": str,
                                "Date": str,
                                "InterPro_accession": str,
                                "InterPro_discription": str,
                                "GO_annotations": str,
                                "Pathways": str
                        })
        protein_data = protein_data.repartition(npartitions=32).reset_index()
        return protein_data

    def filter_data(self, protein_data):
        ''' function to filter the data '''
        # remove rows with no interPro_accession and Protein_accession values
        protein_data = protein_data[(protein_data.InterPro_accession!="-")
                                    & (protein_data.Protein_accession!="-")]
        return protein_data[['Protein_accession', 'Seq_len', 'Start', 'Stop', 'InterPro_accession']]

    def label_size(self, protein_data):
        '''
        function to label the data according to size
        '''
        # calculate size
        protein_data['Size'] = ((protein_data['Stop'] - protein_data['Start'])
                                / protein_data['Seq_len']) * 100
        sizes = protein_data['Size'].to_dask_array()
        protein_data['Features'] = \
        da.where(sizes > 90, 'large','small')
        protein_data['Features'].compute()
        return protein_data

    def clean_data(self, protein_data):
        '''
        function to clean the data and remove unwanted rows
        '''
        # grouby the proteins that have 2 unique feature types (large and small)
        grouped = protein_data.groupby('Protein_accession')['Features'].nunique() == 2
        # map the grouping to protein accession id
        groups = protein_data['Protein_accession'].map(grouped)
        # filter out proteins that do not have 2 unique feature types
        protein_data = protein_data[groups == True]
        return protein_data.reset_index(drop=True)

    def get_model_df(self, protein_data):
        '''
        funtion to create input and label data
        '''
        # extract the label data as a dataframe
        label_df = protein_data.groupby(['Protein_accession'])['InterPro_accession',
                                                             'Size'].max().compute()
        # extract the small feautures as a dataframe
        small = protein_data[protein_data.Features == 'small']
        # aggregate by count of small features per protein
        input_df = small.groupby(['Protein_accession','InterPro_accession']).agg(
                            {'InterPro_accession': 'count'}).compute()
        # make into a dataframe
        input_df = pd.DataFrame(input_df)
        # rename columns
        input_df.rename(columns={'InterPro_accession':'Count'}, inplace=True)
        # pivot the dataframe
        input_df = input_df.pivot_table(index ='Protein_accession',
                                            columns ='InterPro_accession',
                                            values='Count',
                                            fill_value = 0)

        # merge X and y features into a model data fram
        model_df = input_df.merge(label_df, on='Protein_accession')
        # change the datatype of the label
        model_df['InterPro_accession'] = model_df['InterPro_accession'].astype('category')
        # drop size
        model_df.drop('Size', axis=1, inplace=True)
        # creating instance of labelencoder
        labelencoder = LabelEncoder()
        # Assigning numerical values
        model_df['InterPro_accession'] = labelencoder.fit_transform(model_df['InterPro_accession'])
        pd.DataFrame.to_csv(model_df, self.outpath)
        