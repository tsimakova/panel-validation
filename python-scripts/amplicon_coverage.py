import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


def create_table(input_file: pathlib.PosixPath) -> pd.DataFrame:
    """
    :param input_file: The path to VariFind or Oncoscope coverage analysis results
    :return: Sorted pd.DataFrame
    """
    # Make a pandas dataframe from input data:
    data = pd.read_csv(input_file, sep="\t")
    # Count the total number of reads, relative number of reads and sort:
    total_amp_number = data["total_reads"].sum()
    data["amp_proc"] = data["total_reads"] / total_amp_number
    data_sorted = data.sort_values(by=['amp_proc'])
    data_sorted['amp_serial_num'] = np.arange(data.shape[0])
    return data_sorted


def lin_regression(sample_name: str, data_sorted: pd.DataFrame, threshold: int = 0.85) -> np.ndarray:
    """
    :param sample_name: use sample name in case of low R2 coefficient
    :param data_sorted: Sorted pd.Dataframe
    :param threshold: Threshold for R2 in linear regression
    :return: np.array containing predicted number of amplicon reads
    """
    # Sklearn linear regression:
    x = np.array(data_sorted['amp_serial_num']).reshape(-1, 1)
    y = np.array(data_sorted['amp_proc']).reshape(-1, 1)
    lin_reg = LinearRegression()
    lin_reg.fit(x, y)
    y_predict = lin_reg.predict(x)
    r2 = r2_score(y, y_predict)
    # Raise an error and stop if the observed R2 score is less than R2 threshold:
    if r2 >= threshold:
        return y_predict
    else:
        print(f"R2 for {sample_name} coverage results is less than threshold {threshold}")
        return "stop the loop"


def add_prediction(data_sorted: pd.DataFrame, y_predict: np.ndarray, under_ratio: int = 0.5, over_ratio: int = 1.3):
    """
    :param data_sorted: Sorted pd.DataFrame
    :param y_predict: Predicted number of amplicon reads
    :param under_ratio: The ratio of observed to predicted number of reads for undercovered amplicons
    :param over_ratio: The ratio of observed to predicted number of reads for overcovered amplicons
    :return: sorted pd.Dataframe and tables containing the undercovered and overcovered amplicons
    """
    # Add temporal columns for predicted number of amplicon reads and ratio
    data_sorted["amp_proc_predict"] = y_predict
    data_sorted["ratio"] = abs(data_sorted["amp_proc"] / data_sorted["amp_proc_predict"])
    under_amplicons = data_sorted.loc[data_sorted["ratio"] < under_ratio]
    over_amplicons = data_sorted.loc[data_sorted["ratio"] > over_ratio]
    return data_sorted, under_amplicons, over_amplicons


def amp_scatterplot(sample_name: str, data_sorted: pd.DataFrame, under_amplicons: pd.DataFrame,
                    over_amplicons: pd.DataFrame, output_dir: pathlib.PosixPath, figure_width: float,
                    figure_height: float):
    """
    :param sample_name: sample name for output files
    :param data_sorted: sorted pd.DataFrame
    :param under_amplicons: pd.DataFrame with undercovered amplicons
    :param over_amplicons: pd.DataFrame with overcovered amplicons
    :param output_dir: The path to output files
    :param figure_width: The width of the linear regression plot
    :param figure_height: The height of the linear regression plot
    """
    # Save amplicon coverage scatter plot:
    sns.set_style("darkgrid")
    sns.set(rc={'figure.figsize': (figure_width, figure_height)})
    sns.scatterplot(data=data_sorted, x="amp_serial_num", y="amp_proc", color='grey', linewidth=0, size=3)
    sns.scatterplot(data=under_amplicons, x='amp_serial_num', y='amp_proc', color='red', linewidth=0, size=3)
    sns.scatterplot(data=over_amplicons, x='amp_serial_num', y='amp_proc', color='green', linewidth=0, size=3)
    sns.lineplot(data=data_sorted, x="amp_serial_num", y="amp_proc_predict", color="purple")
    plt.title(f"Profile of amplicon coverage ({sample_name} results)", size=20)
    plt.xlabel('Amplicon number', size=15)
    plt.ylabel('Amplicon coverage', size=15)
    plt.legend([], [], frameon=False)
    plt.savefig(f"{output_dir}/{sample_name}_amplicon_coverage_scatterplot.png")
    plt.clf()


def create_output_table(sample_name: str, under_amplicons: pd.DataFrame, over_amplicons: pd.DataFrame,
                        output_dir: pathlib.PosixPath):
    """
    :param sample_name: sample name for output files
    :param under_amplicons: pd.DataFrame with undercovered amplicons
    :param over_amplicons: pd.DataFrame with overcovered amplicons
    :param output_dir: The path to the output files directory
    """
    # Drop the temporal columns:
    under_amplicons = under_amplicons.drop(["amp_proc", "amp_serial_num", "amp_proc_predict", "ratio"], axis=1)
    under_amplicons = under_amplicons.sort_index()
    over_amplicons = over_amplicons.drop(["amp_proc", "amp_serial_num", "amp_proc_predict", "ratio"], axis=1)
    over_amplicons = over_amplicons.sort_index()
    # Save tables with under- and overcovered amplicons:
    under_amplicons.to_csv(f"{output_dir}/{sample_name}_undercovered_amplicons.txt", sep="\t", index=False)
    over_amplicons.to_csv(f"{output_dir}/{sample_name}_overcovered_amplicons.txt", sep="\t", index=False)


def main(input_files, threshold, under_ratio, over_ratio, output_dir, figure_width, figure_height):
    for el in input_files:
        # Get an input file name for output file names
        sample_name = str(el).split("/")[-1].split(".tsv")[0]
        table = create_table(el)
        predictions = lin_regression(sample_name, table, threshold)
        if str(predictions) != "stop the loop":
            table, under, over = add_prediction(table, predictions, under_ratio, over_ratio)
            amp_scatterplot(sample_name, table, under, over, output_dir, figure_width, figure_height)
            create_output_table(sample_name, under, over, output_dir)
        else:
            continue


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for searching the undercovered amplicons")
    parser.add_argument("-t", "--threshold", type=float, default=0.85, help="Threshold for R2 in linear regression")
    parser.add_argument("-u", "--under_ratio", type=float, default=0.5,
                        help="The ratio of observed to predicted number of reads for undercovered amplicons")
    parser.add_argument("-o", "--over_ratio", type=float, default=1.3,
                        help="The ratio of observed to predicted number of reads for overcovered amplicons")
    parser.add_argument("-i", "--input_files", nargs="+", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to input files")
    parser.add_argument("-d", "--output_dir", type=lambda p: pathlib.Path(p).absolute(),
                        help="The path to the output files directory")
    parser.add_argument("-w", "--figure_width", type=float, default=10,
                        help="The width of the linear regression plot")
    parser.add_argument("-e", "--figure_height", type=float, default=6,
                        help="The height of the linear regression plot")
    args = parser.parse_args()
    main(input_files=args.input_files, threshold=args.threshold, under_ratio=args.under_ratio,
         over_ratio=args.over_ratio, output_dir=args.output_dir, figure_width=args.figure_width,
         figure_height=args.figure_height)
