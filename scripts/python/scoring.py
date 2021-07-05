#!/usr/bin/env python3

# SCORING & SCORING HELPER FUNCTIONS
# 2021 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

"""
DESCRIPTION
"""

import warnings
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# get relevant features from a dataframe of features (threshold based)
def get_relevant_features(features,
                          diff_threshold = 0.05,
                          active_threshold = 0.1,
                          inactive_threshold = 0.1):

    """
    -- DESCRIPTION --
    """

    f = features.loc[(features["DIFFERENCE"] >= diff_threshold) & (features["ACTIVE_FREQUENCY"] >= active_threshold) & (features["INACTIVE_FREQUENCY"] >= inactive_threshold)]
    return f

# get feature impact e.g. if feature is more often present in actives or inactives
def get_feature_impact(features):

    """
    -- DESCRIPTION --
    """

    positives = []
    negatives = []

    for index, row in features.iterrows():
        if row["ACTIVE_FREQUENCY"] > row["INACTIVE_FREQUENCY"]:
            positives.append(row["INTERACTION"])
        else:
            negatives.append(row["INTERACTION"])

    if (len(positives) + len(negatives)) == features.shape[0]:
        return [positives, negatives]
    else:
        warnings.warn("Not all interactions could be assigned!", UserWarning)
        return [positives, negatives]

# scoring function (has to be applied row-wise to a dataframe)
def score(data_row,
          positives,
          negatives,
          strategy = "+"):

    """
    -- DESCRIPTION --
    """

    score = 0
    # strategy 1
    if strategy == "+":
        for interaction in positives:
            if data_row[interaction] > 0:
                score = score + 1
    # strategy 2
    elif strategy == "++":
        for interaction in positives:
            score = score + data_row[interaction]
    # strategy 3
    elif strategy == "+-":
        for interaction in positives:
            if data_row[interaction] > 0:
                score = score + 1
        for interaction in negatives:
            if data_row[interaction] > 0:
                score = score - 1
    # strategy 4
    elif strategy == "++--":
        for interaction in positives:
            score = score + data_row[interaction]
        for interaction in negatives:
            score = score - data_row[interaction]
    else:
        pass

    return score

# get optimal cutoff value based on accuracy
def get_cutoff(labels,
               scores):

    """
    -- DESCRIPTION --
    """

    best_cutoff = 0
    best_accuracy = 0

    for i in range(0, max(scores) + 1, 1):
        predicted = ["active" if score >= i else "inactive" for score in scores]
        agreement = 0
        for j, label in enumerate(labels):
            if label == predicted[j]:
                agreement = agreement + 1
        accuracy = agreement / len(labels)
        if accuracy >= best_accuracy:
            best_accuracy = accuracy
            best_cutoff = i

    return [best_cutoff, best_accuracy]

# calculate metrics for specific cutoff
def test_cutoff(labels,
                scores,
                cutoff,
                metric = "accuracy"):

    """
    -- DESCRIPTION --
    """

    predicted = ["active" if score >= cutoff else "inactive" for score in scores]
    agreement = 0
    FP = 0
    N = 0

    for j, label in enumerate(labels):
        if label == predicted[j]:
            agreement = agreement + 1
        if label == "inactive":
            N = N + 1
        if predicted[j] == "active":
            if label == "inactive":
                FP = FP + 1

    accuracy = agreement / len(labels)
    fpr = FP / N

    if metric == "accuracy":
        return accuracy
    elif metric == "fpr":
        return fpr
    else:
        return None

# get quality metrics for a dataframe with score and labels
def get_metrics(dataframe,
                cutoff):

    """
    -- DESCRIPTION --
    Ya: Yield of actives
    EF: Enrichment Factor
    REF: Relative Enrichment Factor
        Source:
        https://jcheminf.biomedcentral.com/track/pdf/10.1186/s13321-016-0189-4.pdf
        mapping: N = N
                 Ns = n
                 n = A
                 ns = TP
    """

    TP, n, A, N = 0, 0, 0, 0

    for index, row in dataframe.iterrows():
        if row["SCORE"] >= cutoff:
            n = n + 1
        if (row["SCORE"] >= cutoff) and (row["LABEL"] == "active"):
            TP = TP + 1
        if row["LABEL"] == "active":
            A = A + 1
        N = N + 1

    return {"Ya": TP / n, "EF": (TP / n) / (A / N), "REF": (100 * TP) / min(n, A)}

# plot confusion matrix
def plot_confusion_matrix(cm,
                          classes,
                          normalize = False,
                          title = "Confusion matrix",
                          cmap = plt.cm.Blues):

    """
    -- DESCRIPTION --
    """

    plt.imshow(cm, interpolation="nearest", cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation = 45)
    plt.yticks(tick_marks, classes)

    if normalize:
        cm = cm.astype("float") / cm.sum(axis = 1)[:, np.newaxis]
        print("Normalized confusion matrix")

    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, cm[i, j],
                 horizontalalignment = "center",
                 color = "white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylim(len(cm) - 0.5, -0.5)
    plt.ylabel("True label")
    plt.xlabel("Predicted label")

    return cm
