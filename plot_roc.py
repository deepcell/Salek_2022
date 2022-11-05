"""Script to take predictions and labels and plot an ROC Curve.

Usage:

python plot_roc.py --run_ids_file runs.txt --predictions_folder preds_parquets/ --labels_file labels.parquet \
  --pos_classes NSCLC LIVER_CARCINOMA --pos_labels NSCLC --neg_labels NSCLC
"""

import argparse
import logging
import os
import sys
import tqdm

from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import modin.pandas as pd
import ray


os.environ["MODIN_ENGINE"] = "ray"
os.environ["RAY_verbose_spill_logs"] = "0"
ray.init(runtime_env={'env_vars': {'__MODIN_AUTOIMPORT_PANDAS__': '1'}})
logging.basicConfig(stream=sys.stdout, level=logging.INFO)


def read_predictions(run_ids: list[str], predictions_folder_dir: str, pos_classes: list[str]) -> pd.DataFrame:
    """Read and process inputted run predictions parquet files.

    Args:
      run_ids: List of RunIDs strings.
      predictions_folder_dir: Folder directory that contains prediction parquet files.
      pos_classes: List of positive model output classes.

    Returns:
      Pandas DataFrame that contains predictions.

    """

    def extract_max_pos_class(class_list: list[dict[str, int]], pos_classes: list[str]) -> float:
        """Extract max positive class probablity from a list of cell class prediction probabilities.

        Args:
          class_list: List of dictionaries that contain cell class to prediction probabilities.
          pos_classes: List of positive model output classes.

        Returns:
          Float that represents the max postive class.

        """

        probs = []
        for discover_class in pos_classes:
            for class_dict in class_list:
                if class_dict["cell_class"] == discover_class:
                    probs.append(class_dict["probability"])

        if len(probs) == 0:
            raise ValueError

        return max(probs)

    predictions = pd.DataFrame()
    for run_id in tqdm.tqdm(run_ids, unit="RunID"):
        source_df = pd.read_parquet(f"{predictions_folder_dir}/{run_id}.parquet")
        predictions = pd.concat([predictions, source_df]).reset_index(drop=True)

    predictions["pos"] = [extract_max_pos_class(class_list, pos_classes) for class_list in predictions["predictions"]]
    return predictions


def read_labels(run_ids: list[str], labels_file: str) -> pd.DataFrame:
    """Read the labels parquet file and filter by run_ids.

    Args:
      run_ids: List of RunIDs to read labels from.
      labels_file: Location of the labels parquet file.

    Returns:
      Pandas dataframe of labels

    """
    labels_df = pd.read_parquet(labels_file)
    return labels_df[labels_df["run_id"].isin(run_ids)]


def get_label(label: str, sample_type: str):
    """Helper function for getting an implied label if sample type is ADULT_BLOOD.

    Args:
        label: String representing label.
        sample_type: String representing sample type.

    """
    if sample_type == "ADULT_BLOOD":
        return "WBC"
    else:
        return label


def process_metrics(
    merged_df: pd.DataFrame, thresh: float, pos_labels: list[str], neg_labels: list[str]
) -> tuple[float, float, float]:
    """Process metrics for a given dataframe and thresholds.

    Args:
      df: Merged dataframe that contains predictions and labels.
      thresh: Threshold for plotting ROC Curve.
      pos_labels: List of positive labels.
      neg_labels: List of negative labels.

    Returns:
      Tuple of the threshold, false positive rate and recall.

    """

    def fpr(fp: float, tn: float) -> float:
        """Calculate False Positive Rate.

        Args:
          fp: Sum of all False Positives.
          tn: Sum of all True Negatives.

        Returns:
          Float representing False Positive Rate.

        """
        if fp + tn > 0:
            fpr = fp / (fp + tn)
        else:
            fpr = -1
        return fpr

    def recall(tp: float, fn: float) -> float:
        """Calculate Recall.

        Args:
          tp: Sum of all True Positives.
          fn: Sum of all False Negatives.

        Returns:
          Float representing Recall.

        """
        if tp + fn > 0:
            recall = tp / (tp + fn)
        else:
            recall = -1
        return recall

    def calculate_metrics(
        row: pd.Series,
        thresh: float,
        pos_labels: list[str],
        neg_labels: list[str],
    ) -> tuple[int, int, int, int]:
        """Calulate total amount of false postives, true positives, false negatives and true negatives.

        Args:
          row: Pandas series that represents a row of the merged DataFrame.
          thresh: Threshold.
          pos_labels: List of positive labels.
          neg_labels: List of negative labels.

        Returns:
          Tuple of total false positives, true positives, false negatives and true negatives.

        """
        FP, TP, FN, TN = 0, 0, 0, 0

        if row["label"] in pos_labels:
            if row["pos"] >= thresh:
                TP = 1
            else:
                FN = 1

        if row["label"] in neg_labels:
            if row["pos"] >= thresh:
                FP = 1
            else:
                TN = 1

        return (FP, TP, FN, TN)

    merged_df["metrics"] = [calculate_metrics(row, thresh, pos_labels, neg_labels) for _, row in merged_df.iterrows()]
    merged_df["FP"] = [row[0] for row in merged_df["metrics"]]
    merged_df["TP"] = [row[1] for row in merged_df["metrics"]]
    merged_df["FN"] = [row[2] for row in merged_df["metrics"]]
    merged_df["TN"] = [row[3] for row in merged_df["metrics"]]

    return (
        thresh,
        fpr(merged_df["FP"].sum(), merged_df["TN"].sum()),
        recall(merged_df["TP"].sum(), merged_df["FN"].sum()),
    )


def plot_data(roc_df: pd.DataFrame, pos_classes: list[str], image_filename: str = "ROC_Curve.png") -> None:
    """Computes and saves ROC Curve as an image.

    Args:
      roc_df: Dataframe that contains columns for FPR, Recall and Threshold.
      pos_classes: List of positive model output classes.

    """
    plt.title("ROC Curve", fontsize=18)
    plt.xlabel("FPR", fontsize=18)
    plt.ylabel("TPR", fontsize=18)
    plt.plot(roc_df["fpr"], roc_df["recall"], label=f"{pos_classes}")
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=18)
    plt.savefig(image_filename)


def compute_roc(
    run_ids_file: str,
    predictions_folder: str,
    labels_file: str,
    pos_classes: list[str],
    pos_labels: list[str],
    neg_labels: list[str],
) -> None:
    """Computes and saves ROC Curve.

    Args:
      run_ids_file: Text file location for Run IDs to be used for plotting.
      predictions_folder: Folder that contains the parquet files that contain predictions for each run.
      pos_classes: List of positive model output classes.
      pos_labels: List of positive labels.
      neg_labels: List of negative labels.

    """

    range_0 = [0, 0.00001, 0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.999, 0.9999, 0.99995, 0.99999]
    thresholds = np.sort(
        np.append(
            range_0, np.append(np.arange(0.05, 0.95, 0.05, dtype=float), np.arange(0.96, 0.99, 0.01, dtype=float))
        )
    )

    logging.info("Reading Run IDs")
    run_ids = open(run_ids_file).read().splitlines()

    logging.info("Reading Predictions")
    predictions = read_predictions(run_ids, predictions_folder, pos_classes)

    logging.info("Reading Labels")
    labels = read_labels(run_ids, labels_file)

    logging.info("Merging Predictions & Labels")
    merged_df = pd.merge(predictions, labels, on=["cell_id"], how="left")
    merged_df["label"] = [get_label(cell_id["label"], cell_id["sample_type"]) for _, cell_id in merged_df.iterrows()]
    merged_df = merged_df.dropna(subset=["label"])

    logging.info("Processing Metrics")
    metrics = []
    for thresh in tqdm.tqdm(thresholds, unit="threshold"):
        metrics.append(process_metrics(merged_df, thresh, pos_labels, neg_labels))

    fpr = pd.Series([metric[1] for metric in metrics])
    recall = pd.Series([metric[2] for metric in metrics])
    roc_df = pd.DataFrame({"threshold": thresholds, "fpr": fpr, "recall": recall})

    logging.info("Plotting Curve")
    plot_data(roc_df, pos_classes)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot ROC Curve from input data")
    parser.add_argument(
        "--run_ids_file", type=str, default="runs.txt", help="file directory which contains the run ids"
    )
    parser.add_argument(
        "--predictions_folder",
        type=str,
        default="pred_parquets",
        help="folder directory which contains the predictions",
    )
    parser.add_argument(
        "--labels_file", type=str, default="labels.parquet", help="file directory which contains the labels"
    )
    parser.add_argument("--pos_classes", nargs="+", type=str, default=["NSCLC", "LIVER_CARCINOMA"], help="")
    parser.add_argument("--pos_labels", nargs="+", type=str, default=["NSCLC"], help="")
    parser.add_argument("--neg_labels", nargs="+", type=str, default=["NSCLC"], help="")

    args = parser.parse_args()

    compute_roc(
        args.run_ids_file,
        args.predictions_folder,
        args.labels_file,
        args.pos_classes,
        args.pos_labels,
        args.neg_labels,
    )
