from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn import tree
import pandas as pd
from sklearn import preprocessing
import pickle
import sys
import matplotlib.pyplot as plt
import sklearn
from sklearn.svm import SVC


def main():
    if len(sys.argv) != 5:
        print('input-dir output-dir test_pct normalize')
        sys.exit(1)

    inputDir =sys.argv[1]
    outputDir = sys.argv[2]
    test_pct = float(sys.argv[3])
    normalize = bool(sys.argv[4])

    le = preprocessing.LabelEncoder()

    print('inputDir = ' + inputDir)
    print('outputdir = ' + outputDir)
    print('test % = ' + str(test_pct))

    # build model based on f5 copd
    df1 = pd.read_csv(inputDir + '/F5/REPORT/f5_chr17_brca1_copdhmb_report.tsv', sep='\t')
    df2 = pd.read_csv(inputDir + '/F5/REPORT/f5_chr13_brca2_copdhmb_report.tsv', sep='\t')
    dfAll = pd.concat([df1, df2], axis=0)

    #features = ['popFreq', 'AF', 'hail_hweafp', 'chisquare', 'class']
    features = ['popFreq', 'cohortFreq', 'class', 'hail_hweafp', 'chisquare']


    #model_1_dt = buildModel(df1, features, tree.DecisionTreeClassifier(max_depth=2), test_pct)
    #model_2_dt = buildModel(df2, features, tree.DecisionTreeClassifier(max_depth=2), test_pct)

    #runMe(features, le, df1, normalize, test_pct)



    model_1_rf = buildModel(df1, features, RandomForestClassifier(n_estimators=10), test_pct, normalize, le)
    model_2_rf = buildModel(df2, features, RandomForestClassifier(n_estimators=10), test_pct, normalize, le)
    model_all_rf = buildModel(dfAll, features, RandomForestClassifier(n_estimators=10), test_pct, normalize, le)

    # save RF models
    f = open(outputDir + '/brca1_rf_model', 'wb')
    pickle.dump(model_1_rf, f)
    f.close()
    f = open(outputDir + '/brca2_rf_model', 'wb')
    pickle.dump(model_2_rf, f)
    f.close()
    f = open(outputDir + '/brca_all_model', 'wb')
    pickle.dump(model_all_rf, f)
    f.close()

    # predict for f8 gru + hmb
    brca1_f8_report_DF = pd.read_csv(inputDir + '/F8/f8_chr17_brca1_gruhmb_report.tsv', header=0, sep='\t')
    brca2_f8_report_DF = pd.read_csv(inputDir + '/F8/f8_chr13_brca2_gruhmb_report.tsv', header=0, sep='\t')
    brca_all_f8_report_DF = pd.concat([brca1_f8_report_DF, brca2_f8_report_DF], axis=0)

    brca1_f8_predictions = getPredictions(brca1_f8_report_DF, model_1_rf, features, normalize, le)
    brca2_f8_predictions = getPredictions(brca2_f8_report_DF, model_2_rf, features, normalize, le)
    brca_all_f8_predictions = getPredictions(brca_all_f8_report_DF, model_2_rf, features, normalize, le)


    # save to disk
    brca1_f8_predictions.to_csv(outputDir + '/F8/f8_chr17_brca1_predictions_report.tsv', sep='\t')
    brca2_f8_predictions.to_csv(outputDir + '/F8/f8_chr13_brca2_predictions_report.tsv', sep='\t')
    brca_all_f8_predictions.to_csv(outputDir + '/F8/f8_all_predictions_report.tsv', sep='\t')

    brca1_f8_predictions = pd.read_csv(inputDir + '/F8/f8_chr17_brca1_predictions_report.tsv', header=0, sep='\t')
    brca2_f8_predictions = pd.read_csv(inputDir + '/F8/f8_chr13_brca2_predictions_report.tsv', header=0, sep='\t')
    brca_all_f8_predictions = pd.read_csv(inputDir + '/F8/f8_all_predictions_report.tsv', header=0, sep='\t')


    # shave off variants from df that are predicted to be in gnomad and save to disk
    brca1_in = brca1_f8_predictions.loc[(brca1_f8_predictions['gnomadPrediction'] == True)]['variant']
    brca2_in = brca2_f8_predictions.loc[(brca2_f8_predictions['gnomadPrediction'] == True)]['variant']
    brca_all_in = brca_all_f8_predictions.loc[(brca_all_f8_predictions['gnomadPrediction'] == True)]['variant']

    brca1_in.to_csv(outputDir + '/F8/brca1_in.txt', index=False, header=0)
    brca2_in.to_csv(outputDir + '/F8/brca2_in.txt', index=False, header=0)
    brca_all_in.to_csv(outputDir + '/F8/brca_all_in.txt', index=False, header=0)


def runMe(features, le, df, normalize, testPctg):
    for f in features:
        if isinstance(df[f].iloc[0], str):
            df[f] = le.fit_transform(df[f])

    # define features and labels
    if normalize:
        X = preprocessing.normalize(df[features])
        Y = df[['inGnomad']]

    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=testPctg)

    for c in [0.01, 0.1, 1, 10, 100, 1000]:
        model = trainMe(cValue=c, gValue=1, x=x_train, y=y_train, kernel='rbf')
        y_predict_test = model.predict(x_test)
        print('rbf test accuracy for c = ' + str(c) + ' = ' + str(accuracy_score(y_test, y_predict_test)))
        print('score = ' + str(model.score(x_test, y_test)))

    for c in [0.01, 0.1, 1, 10, 100, 1000]:
        model = trainMe(cValue=c, gValue=1, x=x_train, y=y_train, kernel='linear')
        y_predict_test = model.predict(x_test)
        print('linear test accuracy for c = ' + str(c) + ' = ' + str(accuracy_score(y_test, y_predict_test)))

    for g in [0.01, 0.1, 1, 10, 100, 1000]:
        model = trainMe(cValue=1, gValue=g, x=x_train, y=y_train, kernel='rbf')
        y_predict_test = model.predict(x_test)
        print('rbf test accuracy for g = ' + str(g) + ' = ' + str(accuracy_score(y_test, y_predict_test)))
        print('score = ' + str(model.score(x_test, y_test)))

def trainMe(cValue, gValue, x, y, kernel):
    model = SVC(C = cValue, gamma = gValue, kernel=kernel)
    model.fit(x, y)
    return model


def getPredictions(df, model, features, normalize, le):
    predictions = list()

    transformedFeatures = list()
    for f in features:
        if isinstance(df[f].iloc[0], str):
            df[f] = le.fit_transform(df[f])
            transformedFeatures.append(f)

    if normalize:
        X = preprocessing.normalize(df[features])

    for i in range(len(X)):
        variant = df.iloc[i]
        predicted = model.predict([list(variant[features])])
        predictions.append(predicted[0])
    df['gnomadPrediction'] = predictions

    for f in transformedFeatures:
        df[f] = le.inverse_transform(df[f])

    return df

def buildModel(df, features, model, testPctg, normalize, le):
    print('building model type ' + str(type(model)))

    # transform boolean to numerica
    for f in features:
        if isinstance(df[f].iloc[0], str):
            df[f] = le.fit_transform(df[f])

    # define features and labels
    if normalize:
        X = preprocessing.normalize(df[features])
        Y = df[['inGnomad']]

    # find fraction of true false in data set
    Y_true = len(Y[Y['inGnomad'] == True])
    Y_false = len(Y[Y['inGnomad'] == False])
    print("fraction of false in df1 = " + str(Y_false/(Y_true + Y_false)))

    # split data
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=testPctg)

    # build model
    model.fit(X_train, y_train)

    # get predictions
    y_predict_train = model.predict(X_train)

    # get overall test accuracy
    print('normalize = ' + str(normalize) + ' train accuracy = ' + str(accuracy_score(y_train, y_predict_train)))
    y_predict_test = model.predict(X_test)
    print('normalize = ' + str(normalize) + ' test accuracy = ' + str(accuracy_score(y_test, y_predict_test)))

    # get false positive rate (sensitivity) and false positive rate (specificity)
    fp = 0
    fn = 0
    numCorrect = 0
    yActualList = y_test['inGnomad'].tolist()
    for i in range(len(yActualList)):
        if yActualList[i] == y_predict_test[i]:
            numCorrect += 1
        elif yActualList[i] == True and y_predict_test[i] == False:
            fn += 1
        else:
            fp +=1

    accuracy = (1.0 * numCorrect) / (1.0 * len(yActualList))
    fpr = 1 - (1.0 * fp) / (1.0 * len(yActualList))
    fnr = 1 - (1.0 * fn) / (1.0 * len(yActualList))
    print('accuracy = ' + str(accuracy))
    print('fpr = ' + str(fpr))
    print('fnr = ' + str(fnr))
    print('sum acc + fpr + fnr = ' + str(accuracy + fpr + fnr))


    if isinstance(model, sklearn.tree.DecisionTreeClassifier):
        plt.figure()
        tree.plot_tree(model,filled=True)
        plt.show()
        #plt.savefig('/content/gdrive/My Drive/RESEARCH/TOPMED/brca1-dtree.eps',format='eps',bbox_inches = "tight")

    elif isinstance(model, sklearn.ensemble.forest.RandomForestClassifier):
        for i in range(len(features)):
            print('feature importance: ')
            print(features[i] + ' : ' + str(model.feature_importances_[i]))

    return model

if __name__ == "__main__":
    main()