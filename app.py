import streamlit as st
import pandas as pd
import numpy as np
# import altair as alt
# import pydeck as pdk

# SETTING PAGE CONFIG TO WIDE MODE
# st.beta_set_page_config(layout="wide")

st.title('Predicting covid test results')
st.subheader('Updated on 2021/5/12 by Ming Tang')

# input variables
st.subheader('Check if ture:')
x1 = st.checkbox('Cough')
x2 = st.checkbox('Fever')
x3 = st.checkbox('Sore throat')
x4 = st.checkbox('Shortness of breath')
x5 = st.checkbox('Headache')
x6 = st.checkbox('Age > 60')
x7 = st.checkbox('Male')
x8 = st.checkbox('Contact with confirmed cases')
x9 = st.checkbox('Travel aboard recently')

x1 = 1 if x1 else 0
x2 = 1 if x2 else 0
x3 = 1 if x3 else 0
x4 = 1 if x4 else 0
x5 = 1 if x5 else 0
x6 = 1 if x6 else 0
x7 = 1 if x7 else 0
x8 = 1 if x8 else 0
x9 = 1 if x9 else 0

@st.cache
def load_data(url):
    df = pd.read_pickle(url+ '?raw=true')
    return df

url = 'https://github.com/tangming2008/Datasets_for_projects/blob/master/df_CovidTest_10000rows.pkl'
df = load_data(url+ '?raw=true')


X = df[list(df.columns)[0:-1]]
y = df[list(df.columns)[-1:]]

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# classification model
from sklearn.linear_model import LogisticRegression
model = LogisticRegression(C=10000) # C large enough to remove regularization
model.fit(X_train, y_train)
# y_pred = model.predict(X_test)
y_pred = model.predict_proba(X_test)[:,1] >= 0.079

# prediction on the input data
onedata = [x1, x2, x3, x4, x5, x6, x7, x8, x9]
df_onedata = pd.DataFrame([onedata])
probability = model.predict_proba(df_onedata)[:,1]

predict_result = model.predict(df_onedata)[0]
predict_result = 'positive' if predict_result > 0.5 else 'negative'
probability_positive = "{:.0%}".format(probability[0])

st.subheader('Show results:')
'Class prediction (hard classification) :  ', predict_result
'Positive probability (soft classification) : ' , probability_positive

st.subheader('Threshold tuning:')

threshold = st.slider('threshold', min_value = 0.0, max_value = 1.0, step=0.02, value = 0.5)

# adjusted model
y_pred = model.predict_proba(X_test)[:, 1] >= threshold

from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score
accuracy = round(accuracy_score(y_test, y_pred),2)
recall = round(recall_score(y_test, y_pred),2)
precision = round(precision_score(y_test, y_pred),2)
f1 = round(f1_score(y_test, y_pred),2)

df = pd.DataFrame()
df['metric'] = ['accuracy', 'recall', 'precision', 'f1']
df['value'] = [accuracy, recall, precision, f1]

st.subheader("Metric scores")
df

st.subheader("Confusion Matrix")
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
from mlxtend.plotting import plot_confusion_matrix

cm = confusion_matrix(y_test, y_pred)
fig, ax = plot_confusion_matrix(conf_mat=cm, figsize=(3, 3))
st.pyplot(fig)

'data', cm
