''' This ia class that performs machine learning on the model_data file.
To run this, please use the main script assignment6.py
'''
import pickle
import warnings
import joblib
from dask_ml.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from dask_ml.model_selection import train_test_split

warnings.filterwarnings('ignore')

class Train():
    '''
    class to train model and make predictions
    '''
    def __init__(self, clean_data, report_path, model_path):
        self.data = clean_data
        self.report_path = report_path
        self.model_path = model_path
        self.train_features, self.test_features, self.train_labels, self.test_labels = self.train_prep()
    def train_prep(self):
        '''
        funtion to prepare the data for training
        '''
        # select features
        features = self.data.iloc[:,1:-1].values
        labels = self.data.iloc[:, -1].values
        #split data
        train_features, test_features, train_labels, test_labels = train_test_split(
            features,labels, test_size = 0.3, random_state=1)
        return train_features, test_features, train_labels, test_labels

    def train_model(self):
        '''
        function to perform random forest classification
        '''
        with joblib.parallel_backend("dask",
                                     scatter=[self.train_features, self.train_labels]):
            model = RandomForestClassifier( random_state = 42)
            print('Now training')
            model.fit(self.train_features, self.train_labels)
            #predict and evaluate
            pred_labels = model.predict(self.test_features)
            report = f"Accuracy for Random Forest: {self.test_labels, pred_labels}"
            with open(self.report_path, 'a', encoding="utf-8") as out_file:
                out_file.write(report)
                print(report)
            # save model
            filepath = self.model_path
            with open(filepath, 'wb', encoding="utf-8"):
                pickle.dump(model, filepath)
