import numpy as np
import pandas as pd
import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

RUN_IDS_FILE = 'runs_val.txt'
PREDICTION_PARQUET_FOLDER = 'new_parquets'
LABEL_PARQUET_FILE = '~/labels.parquet'
POS_CLASSES = ['NSCLC', 'LIVER_CARCINOMA']
POS_LABELS = ['NSCLC']
NEG_LABELS = ['WBC']


def fpr(fp, tn):
  if (fp + tn > 0):
    fpr = fp/(fp + tn)
  else:
    fpr = -1

  return fpr


def recall(tp, fn):
  if tp + fn > 0:
    recall = tp/(tp+fn)
  else:
    recall = -1

  return recall


def create_pos(df, pos):
    
    def extract_class(class_list, discover_classes):
        probs = []
        for discover_class in discover_classes:
            for class_dict in class_list:
                if class_dict['cell_class'] == discover_class:
                    probs.append(class_dict['probability'])
                    
        if len(probs) == 0:
            raise ValueError
            
        return max(probs)
            
    df['pos'] = [extract_class(class_list, pos) for class_list in df['predictions']]
    return df


def read_predictions(run_ids):
    predictions = pd.DataFrame()
    for run_id in tqdm.tqdm(run_ids):
        df = pd.read_parquet(f'{PREDICTION_PARQUET_FOLDER}/{run_id}.parquet')
        predictions = pd.concat([predictions, df]).reset_index()
    predictions = create_pos(predictions, POS_CLASSES)
    return predictions


def read_labels(run_ids):
    df = pd.read_parquet(LABEL_PARQUET_FILE)
    df = df[df['run_id'].isin(run_ids)] 
    return df 

  
def get_label(label, sample_type):
  if sample_type == 'ADULT_BLOOD':
      return 'WBC'
  else:
      return label

    
def process_metrics(df, thresh, pos_labels, neg_labels):
    
    def calculate_metrics(row, thresh, pos_labels, neg_labels):
        FP, TP, FN, TN = 0, 0, 0, 0
        
        if row['label'] in pos_labels:
            if row['pos'] >= thresh:
                TP = 1
            else:
                FN = 1
                
        if row['label'] in neg_labels:
            if row['pos'] >= thresh:
                FP = 1
            else:
                TN = 1
                
        return (FP, TP, FN, TN)
    
    df['metrics'] = [calculate_metrics(row, thresh, pos_labels, neg_labels) for index, row in df.iterrows()]
    df['FP'] = [row[0] for row in df['metrics']]
    df['TP'] = [row[1] for row in df['metrics']]
    df['FN'] = [row[2] for row in df['metrics']]
    df['TN'] = [row[3] for row in df['metrics']]
    
    return (thresh, fpr(df['FP'].sum(), df['TN'].sum()), recall(df['TP'].sum(), df['FN'].sum()))


thresholds = [0, 0.00001, 0.0001, 0.001, 0.01, 0.02, 0.03,0.04, 0.999, 0.9999, 0.99995, 0.99999]
thresholds2 = np.arange(0.05, 0.95, 0.05, dtype=float)
thresholds3 = np.arange(0.96, 0.99, 0.01, dtype=float)

for elem in thresholds2: 
    thresholds.append(elem)
    
for elem in thresholds3:
    thresholds.append(elem)

thresholds = [round(threshold, 5) for threshold in thresholds]
thresholds = sorted(thresholds)


run_ids = open(RUN_IDS_FILE).read().splitlines()
predictions = read_predictions(run_ids)
labels = read_labels(run_ids)
metrics_df = pd.merge(predictions, labels, on=['cell_id'], how='left')
metrics_df['label'] = [get_label(cell_id['label'], cell_id['sample_type']) for _, cell_id in metrics_df.iterrows()]
metrics_df = metrics_df.dropna(subset=['label'])

metrics = []
for thresh in thresholds:
    metrics.append(process_metrics(metrics_df, thresh, POS_LABELS, NEG_LABELS))

fpr = pd.Series([metric[1] for metric in metrics])
recall = pd.Series([metric[2] for metric in metrics])
final_df = pd.DataFrame({'threshold': thresholds, 'fpr': fpr, 'recall': recall})

plt.title("ROC Curve", fontsize=18)
plt.xlabel("FPR", fontsize=18)
plt.ylabel("TPR", fontsize=18)
plt.plot(final_df["fpr"], final_df["recall"], label=f"{POS_CLASSES}")
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=18)
