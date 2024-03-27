import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from sklearn import metrics
from sklearn.model_selection import (
    StratifiedKFold,
    GroupShuffleSplit,
    train_test_split,
    cross_val_score
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, KBinsDiscretizer
from utils import (
    get_full_descriptors_set,
    sel_features
)


# load data
data = pd.read_csv('../data/db_dnazymes.csv', index_col=0)
seq_column_name = 'e'
buffer_column_name = 'buffer'
cofactor_column_name = 'metal_ions'
temperature_column_name = 'temperature'

# encode data
X_full = get_full_descriptors_set(
    df=data,
    seq_column_name=seq_column_name,
    buffer_column_name=buffer_column_name,
    cofactor_column_name=cofactor_column_name,
    temperature_column_name=temperature_column_name
)
y_full = np.log10(data['kobs'])

# train weights
train_inds, val_inds = next(
    GroupShuffleSplit(
        test_size=.15,
        n_splits=2,
        random_state=10
    ).split(X_full, groups=data[seq_column_name])
)

X, y = X_full.iloc[train_inds], y_full[train_inds]
X_val, y_val = X_full.iloc[val_inds], y_full[val_inds]
y = np.array(y)
y_discretized = KBinsDiscretizer(
    n_bins=5,
    encode='ordinal',
    strategy='uniform',
    random_state=0
).fit_transform(y.reshape(-1, 1))

X_train, X_test, y_train, y_test = train_test_split(
    X, y.ravel(),
    test_size=0.2,
    random_state=0,
    stratify=y_discretized
)

y_train_discretized = KBinsDiscretizer(
    n_bins=5,
    encode='ordinal',
    strategy='uniform', random_state=0
).fit_transform(y_train.reshape(-1, 1))

model = Pipeline(
    [
        ('scaler', StandardScaler()),
        ('reg', LGBMRegressor(random_state=0, n_estimators=200, max_depth=6, learning_rate=0.1))
    ]
)

model.fit(X_train[sel_features], y_train, reg__eval_set=(X_test[sel_features], y_test), reg__eval_metric='r2')

# evaluation of the weights
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=0).split(X=X_train, y=y_train_discretized)
scores_cv = cross_val_score(
    model,
    X_train[sel_features],
    y_train,
    scoring='r2',
    cv=skf
)
y_pred = model.predict(X_test[sel_features])
y_val_pred = model.predict(X_val[sel_features])
y_train_pred = model.predict(X_train[sel_features])

r2_test = metrics.r2_score(y_test, y_pred)
r2_val = metrics.r2_score(y_val, y_val_pred)
rmse_test = np.sqrt(metrics.mean_squared_error(y_test, y_pred))
rmse_val = np.sqrt(metrics.mean_squared_error(y_val, y_val_pred))

scores_final = {
    'Q2': scores_cv.mean(),
    'R2_test': r2_test,
    'RMSE_test': rmse_test,
    'R2_val': r2_val,
    'RMSE_val': rmse_val
}
