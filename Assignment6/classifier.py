from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
import pandas as pd
import numpy as np
import dask
import dask.delayed
import time
import pickle
from dask_ml.compose import ColumnTransformer
from dask_ml.impute import SimpleImputer
from dask_ml.preprocessing import StandardScaler, DummyEncoder
from distributed import wait
from dask_ml.metrics import accuracy_score
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from dask_ml.model_selection import train_test_split
from sklearn.pipeline import Pipeline
import dask.array as da
import joblib
from dask.distributed import Client, progress
import warnings
warnings.filterwarnings('ignore')



class train():
    '''
    train model and predict
    '''
    def __init__(self):
        pass
    
    def train_prep(self, data):
        '''
        funtion to prepare the data for training
        '''
        # select features
        features = data.iloc[:,1:-1].values
        # Saving feature names for later use
        #feature_list = list(features.columns)
        # Convert to numpy array
        # features = da.array(features)
        
        # select labels
        labels = data.iloc[:, -1].values
        

        #split data
        train_features, test_features, train_labels, test_labels = train_test_split(features,
                                                                labels, test_size = 0.3, random_state=1)
        return train_features, test_features, train_labels, test_labels
    
    @staticmethod
    def train_model(train_features, test_features, train_labels, test_labels):
        '''
        function to perform random forest classification
        '''
        
      
        with joblib.parallel_backend("dask", scatter=[train_features, train_labels]):
            model = RandomForestClassifier(n_estimators=1000, random_state = 42)
            print('Now training')
            model.fit(train_features, train_labels)
            #predict and evaluate
            pred_labels = model.predict(test_features)
            print("Accuracy for Random Forest: ",accuracy_score(test_labels, pred_labels))
            # report = classification_report(test_labels, pred_labels, output_dict=True)
            
            # report_df = pd.DataFrame(report).transpose()
            # pd.DataFrame.to_csv(report_df, '/homes/fabadmus/programming_3/Programming3/Assignment6/model_report2')
            # save the model 
            filename = 'rf_model.sav'
            pickle.dump(model, open(filename, 'wb'))
       
    # @staticmethod    
    # def model_train(train_features, test_features, train_labels, test_labels):
        
    #     with joblib.parallel_backend('dask'):
            
    #         # model.fit(train_features, train_labels)
    #         param_grid = { 
    #                 'n_estimators': [200, 500],
    #                 'max_features': ['auto', 'sqrt', 'log2'],
    #                 'max_depth' : [4,5,6,7,8],
    #                 'criterion' :['gini', 'entropy']
    #             }
    #         model = RandomForestClassifier(n_estimators=1000, random_state = 42)
    #         CV_rfc = GridSearchCV(estimator=model, param_grid=param_grid, cv= 5)
    #         print('Now training')
            
    #         CV_rfc.fit(train_features, train_labels)
    #         print('Done training')
    
        


        
       
    #     pred_labels = CV_rfc.predict(test_features)
        
    #     print("Accuracy for Random Forest on CV data: ",accuracy_score(test_labels, pred_labels))
        
    #     report = classification_report(test_labels, pred_labels, output_dict=True)
        
    #     report_df = pd.DataFrame(report).transpose()
    #     pd.DataFrame.to_csv(report_df, '/homes/fabadmus/programming_3/Programming3/Assignment6/model_report')
        
    #     # save the model 
    #     filename = 'rf_model.sav'
    #     pickle.dump(CV_rfc, open(filename, 'wb'))




        
def main():

    data = pd.read_csv('/homes/fabadmus/programming_3/Programming3/Assignment6/model_data', )
    model = train()
    train_features, test_features, train_labels, test_labels = model.train_prep(data)
    # print(train_labels)
    client = Client(processes=False, threads_per_worker=4,
                n_workers=1, memory_limit='128GB')
    client
        
    model.train_model(train_features, test_features, train_labels, test_labels )
    
    

if __name__ == '__main__':
    main()