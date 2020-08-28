# python3

# environment:
# conda install numpy scikit-learn pandas scipy
# conda install -c conda-forge matplotlib

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import cycle

import sklearn.metrics as metrics
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score


def get_ROC_coordinates(y_test, y_predict):
    fpr, tpr, ROC_threshold = metrics.roc_curve(y_test, y_predict)
    roc_auc = metrics.auc(fpr, tpr)
    return(fpr, tpr, ROC_threshold, roc_auc)


def get_PRC_coordinates(y_test, y_predict):
    assert isinstance(y_test, list)
    precision, recall, PRC_threshold = precision_recall_curve(y_test, y_predict)
    average_precision = average_precision_score(y_test, y_predict)
    return(precision, recall, PRC_threshold, average_precision)


def main():
    fpr_dict = dict()
    tpr_dict = dict()
    ROC_threshold_dict = dict()
    roc_auc_dict = dict()
    precision_dict = dict()
    recall_dict = dict()
    PRC_threshold_dict = dict()
    average_precision_dict = dict()
    labels_ROC = []
    labels_PRC = []

# read data from pandas
    folder_path = './' # '/ur_folder/path'
    # data_path = folder_path+"Supplementary Table S2.csv"

    # 4 novel gRNAs only:
    data_path = folder_path+"Supplementary Table S3.csv"

    # Replace cells with no scores by 0
    # Skipped the first two rows("Supplementary Table S2/3" and the following blank line
    pd_data0 = pd.read_csv(data_path, skiprows=2)
    pd_data = pd_data0.replace(np.nan, 0, regex=True)

    # read algorithm names from file or provide a list
    algorithm_list = ['elevation', 'CRISTA', 'predictCRISPR',
                      'CFD', 'CNN_std', 'MIT', 'CCTOP', 'Hsu', 'COSMID', 'Cropit']

    # From the CRISPOR paper:
    # the MIT site gives us no score for many off-targets
    # so we're setting it to 0.0 for these
    # it's not entirely correct, as we should somehow treat these as "missing data"
    # this leads to a diagonal line in the plot... not sure how to avoid that

    # Here we treat the missing data the same as the CRISPOR paper. (using 0)

    for name in algorithm_list:
        y_test = pd_data['TrueOT'].tolist()
        y_predict = pd_data[name].tolist()
        fpr_dict[name], tpr_dict[name], ROC_threshold_dict[name], roc_auc_dict[name] = \
            get_ROC_coordinates(y_test, y_predict)
        precision_dict[name], recall_dict[name], PRC_threshold_dict[name], average_precision_dict[name] = \
            get_PRC_coordinates(y_test, y_predict)

    # For high color contrasts
    # elevation-sky blue  (0.35,0.7,0.9)
    # crista-orange (0.9, 0.6, 0)
    # predictCRISPR-black (0,0,0)
    # CFD-reddish purple (0.8, 0.6, 0.7)
    # CNN-STD Vermilion (0.8, 0.4, 0)
    # MIT-Blue (0, 0.45, 0.7)
    # CCTOP-bluish Green (0, 0.6, 0.5)
    # Dash lines: Hsu-blue; COSMID-orange, CROPit-black

    # plot ROC
    lw = 1.5
    color_CUD_ROC = cycle([(0.35,0.7,0.9), (0.9, 0.6, 0), (0,0,0), (0.8, 0.6, 0.7),
                       (0.8, 0.4, 0), (0, 0.45, 0.7), (0, 0.6, 0.5)])
    color_CUD_PRC = cycle([(0.35,0.7,0.9), (0.9, 0.6, 0), (0,0,0), (0.8, 0.6, 0.7),
                       (0.8, 0.4, 0), (0, 0.45, 0.7), (0, 0.6, 0.5)])
    linestyles = cycle(['-','-','-','-','-','-','-','--','--','--'])

    f = plt.figure(1)
    for i, color, linestyle in zip(algorithm_list, color_CUD_ROC, linestyles):
        plt.plot(fpr_dict[i], tpr_dict[i], color=color, lw=lw, linestyle=linestyle, alpha=0.85)
        labels_ROC.append('{0} (AUC = {1:0.2f})'''.format(i, roc_auc_dict[i]))

    plt.title('Receiver-Operating curve')
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r:')
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.legend(labels_ROC, loc=(0.65, 0.01), prop=dict(size=7))

    plt.savefig("ROC.svg", format='svg')
    f.show()

    g = plt.figure(2)

    for i, color, linestyle in zip(algorithm_list, color_CUD_PRC, linestyles):
        plt.plot(recall_dict[i], precision_dict[i], color=color, lw=lw, linestyle=linestyle, alpha=0.85)
        labels_PRC.append('{0} (AUC = {1:0.2f})'''.format(i, average_precision_dict[i]))
    plt.xlim([0, 1])
    # plt.ylim([0, 1])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall curve')
    plt.legend(labels_PRC, loc=(0.65, 0.60), prop=dict(size=7))
    plt.savefig("PRC.svg", format='svg')
    g.show()

    # get AUC table
    print('roc_auc_dict', roc_auc_dict)
    print('average_precision_dict', average_precision_dict)


if __name__ == '__main__':
    main()

