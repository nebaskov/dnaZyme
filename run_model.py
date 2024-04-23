import pickle
import pandas as pd

from utils import get_full_descriptors_set, sel_features

test_data = pd.read_csv('data/test_sample.csv', index_col=0)
seq_column_name = 'e'
buffer_column_name = 'buffer'
cofactor_column_name = 'metal_ions'
temperature_column_name = 'temperature'

X_test_full = get_full_descriptors_set(
    df=test_data,
    seq_column_name=seq_column_name,
    buffer_column_name=buffer_column_name,
    cofactor_column_name=cofactor_column_name,
    temperature_column_name=temperature_column_name
)

pickled_model = pickle.load(open('weights/model_v2.pkl', 'rb'))
predictions = pickled_model.predict(X_test_full[sel_features])
