import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass, field
from tqdm import tqdm
import os, sys
from pathlib import Path

# NOTE: Currently data{num_threads}.txt is structured as follow:
#   *File:* {test_case_id}*_*seqA.txt
#   SCORES: -m {m_score}, -x {x_score}, -o {o_score}, -e {e_score}
#   --------------------------------------------------------------
#   Function                       *[The last blank space before Score should be at Hyperparams.ALGO_NAME_IDX_IDENTIFIER]*Score      Consistent AVG Elapsed Time (s)   STD DEV    AVG Performance (GFLOPS/s)  STD DEV
#   *SW* {algo_name}                                                                                                      {score:int}{is_consistent: str}   {time:f}   {tstd_d:f} {perf:f}                    {pstd_d:f}
#   *SW* ...
#   *{score:int}* <- Here it's mandatory that the first char is a number-like. No other lines can have a first char that is a number-like
# IT IS MANDATORY THAT CHARS INSIDE ** ARE PERFECTLY MATCHED.

class Hyperparams:
    DATA_FILES_PATH = "./raw"
    PLOTS_DIR = "./plots"
    TEST_CASES_LINE_IDENTIFIER = "File"
    ALGO_NAME_LINE_IDENTIFIER = "SW"
    ALGO_NAME_IDX_IDENTIFIER = 30  # that is, the number used to pad prinf()

@dataclass
class PlotMetaData:
    num_threads: int = 0
    test_case: int = 0
    title: str = ""
    plot_dst_dir: str = ""

    algos_name: list[str] = field(default_factory=list)
    algos_time: list[float] = field(default_factory=list)
    algos_time_sd: list[float] = field(default_factory=list)
    algos_perf: list[float] = field(default_factory=list)
    algos_perf_sd: list[float] = field(default_factory=list)


def plot_data(alogs_name: list[str], data_array: list[float], std_dev: list[float], plot_title: str, xlabel: str, plot_dst_dir: str):
    plt.figure(figsize=(6, 6))
    plt.barh(alogs_name, data_array, height=0.9, align='center', color='skyblue', edgecolor='black')
    plt.errorbar(data_array, alogs_name, xerr=std_dev, fmt='none', ecolor='black', capsize=5)

    plt.yticks(rotation=45)
    plt.xlabel(xlabel)
    plt.title(plot_title)
    plt.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout()
    plt.savefig(f"{plot_dst_dir}/{plot_title}.png", dpi=400)
    plt.close()

def check_pmd(pmd: PlotMetaData, data_file_path: str):
    if pmd.num_threads == 0:
        print(f"[ERROR]: Failed to set num_threads (currently is {pmd.num_threads = }) while parsing {data_file_path}. Exiting...")
        sys.exit(1)
    if pmd.test_case == 0:
        print(f"[ERROR]: Failed to set test_case (currently is {pmd.test_case = }) while parsing {data_file_path}. Exiting...")
        sys.exit(1)
    if len(pmd.algos_name) == 0:
        print(f"[ERROR]: The array pmd.algos_name was not set properly {pmd.algos_name = }. Exiting")
        sys.exit(1)
    if len(pmd.algos_time) == 0:
        print(f"[ERROR]: The array pmd.algos_time was not set properly {pmd.algos_time = }. Exiting")
        sys.exit(1)
    if len(pmd.algos_time_sd) == 0:
        print(f"[ERROR]: The array pmd.algos_time_sd was not set properly {pmd.algos_time_sd = }. Exiting")
        sys.exit(1)
    if len(pmd.algos_perf) == 0:
        print(f"[ERROR]: The array pmd.algos_perf was not set properly {pmd.algos_perf = }. Exiting")
        sys.exit(1)
    if len(pmd.algos_perf_sd) == 0:
        print(f"[ERROR]: The array pmd.algos_perf_sd was not set properly {pmd.algos_perf_sd = }. Exiting")
        sys.exit(1)

def parse_data_file(data_file_path: str, filename: str) -> list[PlotMetaData]:
    data_file_path = f"{data_file_path}/{filename}"
    # get num_threads from path
    num_threads: str = filename.split("data")[-1].split(".")[0]
    if not num_threads.isdigit():
        print(f"[ERROR]: Got {num_threads = } while parsing {data_file_path = } | expected number-like. Exiting...")
        sys.exit(1)

    pmd_container: list[PlotMetaData] = []
    pmd = PlotMetaData()
    with open(data_file_path, "r") as df:
        lines_container = df.readlines()

    test_case = 0
    for line in lines_container:
        if line.startswith(Hyperparams.TEST_CASES_LINE_IDENTIFIER):
            test_case: str = line.split(f"{Hyperparams.TEST_CASES_LINE_IDENTIFIER}: ")[1].split("_")[0]
            if not test_case.isdigit():
                print(f"[ERROR]: Got {test_case = } while parsing {line = } | expectend number-like. Exiting...")
                sys.exit(1)
            else:
                pmd.test_case = test_case

        if line.startswith(Hyperparams.ALGO_NAME_LINE_IDENTIFIER):
            algo_name: str = line[:Hyperparams.ALGO_NAME_IDX_IDENTIFIER].strip()
            algo_name = algo_name.replace("O(^2)",r"$\mathcal{O}(n^2)$")
            pmd.algos_name.append(algo_name)

            # A1: assuming the current line is loaded *SW* {algo_name} {score:int} {is_consistent: str} {time:f} {tstd_d:f} {perf:f} {pstd_d:f}
            # starting at idx Hyperparams.ALGO_NAME_IDX_IDENTIFIER will result in the following substring {score:int} {is_consistent: str} {time:f} {tstd_d:f} {perf:f} {pstd_d:f}
            algo_raw_data: list[str] = line[Hyperparams.ALGO_NAME_IDX_IDENTIFIER:].split()
            # if A1 holds, algo_raw_data stores [{score:int} {is_consistent: str} {time:f} {tstd_d:f} {perf:f} {pstd_d:f}]
            pmd.algos_time.append(float(algo_raw_data[2]))
            pmd.algos_time_sd.append(float(algo_raw_data[3]))
            pmd.algos_perf.append(float(algo_raw_data[4]))
            pmd.algos_perf_sd.append(float(algo_raw_data[5]))

        if line[0].isdigit():
            pmd.num_threads = num_threads
            pmd.title = f"CASE: {pmd.test_case}"
            pmd.plot_dst_dir = f"{Hyperparams.PLOTS_DIR}/{filename.replace('.txt', '')}"
            pmd_container.append(pmd)
            check_pmd(pmd, data_file_path)
            del pmd
            pmd = PlotMetaData()
    return pmd_container

def speedup_plot(data:list[PlotMetaData]):
    data.sort(key=lambda pmd:pmd.num_threads)
    algs_data = {name : [] for name in data[0].algos_name}

    for pmd in data:
        for name,time,std in zip(pmd.algos_name,pmd.algos_time,pmd.algos_time_sd):
            algs_data[name].append((time,std,pmd.num_threads))

    plt.figure(figsize=(6, 6))

    for name,alg_data in algs_data.items():
        if "serial" in name.lower(): continue
        # print(name,alg_data)
        times,stds,nthreads = np.array(alg_data,dtype=np.float64).T
        speedup = times[0]/times
        plt.errorbar(nthreads,speedup,yerr=stds, label=name, ecolor='black', capsize=5)
    plt.plot([1,8],[1,8],":")
    plt.xlabel("Number of threads")
    plt.ylabel("Speedup")

    plt.title(f"Speedup plot for case {data[0].test_case}")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout()
    plt.savefig(f"plots/speedup/speedup_{data[0].test_case}.png", dpi=400)
    plt.close()


def main():

    if not os.path.isdir(Hyperparams.DATA_FILES_PATH):
        print(f"[WARNING]: {Hyperparams.DATA_FILES_PATH} does not exist. Autocreating it")
        os.mkdir(f"{Hyperparams.DATA_FILES_PATH}")
    if len(os.listdir(Hyperparams.DATA_FILES_PATH)) == 0:
        print(f"[ERROR]: No data files were found in {Hyperparams.DATA_FILES_PATH}. Edit data.sh to output to the same path as {Hyperparams.DATA_FILES_PATH} and run it. Exiting...")
        sys.exit(1)
    if not os.path.isdir(Hyperparams.PLOTS_DIR):
        print(f"[WARNING]: {Hyperparams.PLOTS_DIR} does not exist. Autocreating it")
        os.mkdir(f"{Hyperparams.PLOTS_DIR}")


    speedup_data = {}
    for filename in tqdm(os.listdir(Hyperparams.DATA_FILES_PATH)):
        filepath = Path(f"./plots/{filename.replace('.txt', '')}")
        if not filepath.exists(): filepath.mkdir()
        PMD_container = parse_data_file(f"{Hyperparams.DATA_FILES_PATH}", filename)
        for pmd in PMD_container:
            key = (pmd.title)
            if key not in speedup_data.keys(): speedup_data[key] = []
            speedup_data[key].append(pmd)
        for pmd in PMD_container:
            plot_data(pmd.algos_name, pmd.algos_time, pmd.algos_time_sd, f"Time | {pmd.title}", "Time [s]", pmd.plot_dst_dir)
            plot_data(pmd.algos_name, pmd.algos_perf, pmd.algos_perf_sd, f"Performance | {pmd.title}", "GFLOPS", pmd.plot_dst_dir)

    speeduppath = Path(f"./plots/speedup")
    if not speeduppath.exists(): speeduppath.mkdir()
    for list_pmd in speedup_data.values():
        speedup_plot(list_pmd)


if __name__ == "__main__":
    main()
