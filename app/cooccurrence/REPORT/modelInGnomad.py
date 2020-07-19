from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn import tree
import pandas as pd
from sklearn import preprocessing


df1 = pd.read_csv('/content/gdrive/My Drive/RESEARCH/TOPMED/brca1DF.tsv', sep='\t')
df2 = pd.read_csv('/content/gdrive/My Drive/RESEARCH/TOPMED/brca2DF.tsv', sep='\t')

# transform boolean to numerica
le = preprocessing.LabelEncoder()
# brca1
df1['class'] = le.fit_transform(df1['class'])
# brca2
df2['class'] = le.fit_transform(df2['class'])

# define features and labels
# brca1
#X1 = df1[['class', 'popFreq', 'cohortFreq', 'aa', 'Aa', 'AA', 'hail_hw_pvalue']]
X1 = df1[['class', 'popFreq', 'cohortFreq', 'hail_hw_pvalue']]
#X1 = df1[['class', 'popFreq', 'cohortFreq']]
Y1 = df1[['inGnomad']]
# brca2
X2 = df2[['class', 'popFreq', 'cohortFreq', 'aa', 'Aa', 'AA', 'hail_hw_pvalue']]
Y2 = df2[['inGnomad']]


# find fraction of true false in data set
# brca1
Y1_true = len(Y1[Y1['inGnomad'] == True])
Y1_false = len(Y1[Y1['inGnomad'] == False])
# brca2
Y2_true = len(Y2[Y2['inGnomad'] == True])
Y2_false = len(Y2[Y2['inGnomad'] == False])
print("fraction of false in df1 = " + str(Y1_false/(Y1_true + Y1_false)))
print("fraction of false in df2 = " + str(Y2_false/(Y2_true + Y2_false)))

# split data
# brca1
X1_train, X1_test, y1_train, y1_test = train_test_split(X1, Y1, test_size=0.01)
# brca2
X2_train, X2_test, y2_train, y2_test = train_test_split(X2, Y2, test_size=0.01)

# build decision tree
# brca1
model1_1 = tree.DecisionTreeClassifier(max_depth=2)
model1_1.fit(X1_train, y1_train)
y1_predict_train = model1_1.predict(X1_train)
print('unnormalized dt train accuracy for brca1 = ' + str(accuracy_score(y1_train, y1_predict_train)))
y1_predict_test = model1_1.predict(X1_test)
print('unnormalized dt test accuracy for brca1 = ' + str(accuracy_score(y1_test, y1_predict_test)))
# brca2
model2_1 = tree.DecisionTreeClassifier(max_depth=2)
model2_1.fit(X2_train, y2_train)
y2_predict_train = model2_1.predict(X2_train)
print('unnormalized dt train accuracy for brca2 = ' + str(accuracy_score(y2_train, y2_predict_train)))
y2_predict_test = model2_1.predict(X2_test)
print('unnormalized dt test accuracy for brca2 = ' + str(accuracy_score(y2_test, y2_predict_test)))

# plot and save trees
import matplotlib.pyplot as plt
# brca1
plt.figure()
tree.plot_tree(model1_1,filled=True)
plt.savefig('/content/gdrive/My Drive/RESEARCH/TOPMED/brca1-dtree.eps',format='eps',bbox_inches = "tight")
# brca2
plt.figure()
tree.plot_tree(model2_1,filled=True)
plt.savefig('/content/gdrive/My Drive/RESEARCH/TOPMED/brca2-dtree.eps',format='eps',bbox_inches = "tight")

# build random forest
# brca1
model1_2 = RandomForestClassifier(n_estimators=100)
model1_2.fit(X1_train, y1_train)
y1_2_predict_train = model1_2.predict(X1_train)
print('unnormalized rf1 train accuracy = ' + str(accuracy_score(y1_train, y1_2_predict_train)))
y1_2_predict = model1_2.predict(X1_test)
print('unnormalized rf1 test accuracy = ' + str(accuracy_score(y1_test, y1_2_predict)))
print(model1_2.feature_importances_)
# brca2
model2_2 = RandomForestClassifier(n_estimators=100)
model2_2.fit(X2_train, y2_train)
y2_2_predict_train = model2_2.predict(X2_train)
print('unnormalized rf2 train accuracy = ' + str(accuracy_score(y2_train, y2_2_predict_train)))
y2_2_predict = model2_2.predict(X2_test)
print('unnormalized rf2 test accuracy = ' + str(accuracy_score(y2_test, y2_2_predict)))
print(model2_2.feature_importances_)
