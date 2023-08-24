
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
import seaborn as sns

df = pd.read_csv('sequences.csv', index_col='id')
df['split_enzyme'] = df['e'].apply(list)
seq_data = df['split_enzyme'].apply(lambda x: ' '.join(x).strip()).tolist()

vocabulary = ['A', 'C', 'G', 'T', 'U']
vectorizer = CountVectorizer(vocabulary=vocabulary)
seq_stats = vectorizer.fit_transform(seq_data)

features = vectorizer.get_feature_names_out()


