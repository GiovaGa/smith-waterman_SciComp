import numpy as np
import os
class Hyperparams:
    # editable params
    SEQ_FOLDER = r"test/input"
    PATH_TO_DATA = f"data/raw/G_DATA"
    PARALLEL_TYPE = 3
    MAX_N_THREADS = 4
    N_RUNS = 5
    MAX_TESTCASE_DIM = 6

    # internals, do not edit
    PATH_TO_MERGED = "merged.txt"
    REPORT_FOLDER = "reports"
    PATH_TO_REPORT = f"{REPORT_FOLDER}/Report_D1.txt"


def avg(array) -> float:
    return sum(array)/len(array)


def std_dev(array, avg):
    acc = 0
    for i in array:
        acc += (i - avg)**2
    return (acc/len(array))**0.5


def transpose_array(array):
    a = np.array(array)
    return list(np.transpose(a))


def get_seq_from_file(path: str) -> str:
    with open(path, "r") as file:
        seq = file.readlines()[1]
    return seq


def create_dummy_files():
    os.system(f"mkdir {Hyperparams.REPORT_FOLDER}")
    os.system(f"touch {Hyperparams.PATH_TO_REPORT}")
    os.system(f"touch {Hyperparams.PATH_TO_MERGED}")


def cleanup_dummy_files():
    os.system(f"rm -rf {Hyperparams.REPORT_FOLDER}")
    os.system(f"rm -f {Hyperparams.PATH_TO_MERGED}")


def run_impl(seqA_path: str, seqB_path: str, n_threads, c_exec_filename="./smith_waterman"):
    seq_a = get_seq_from_file(seqA_path)
    seq_b = get_seq_from_file(seqB_path)
    with open("merged.txt", "w") as ofile:
        ofile.write(f"Q:\t{seq_a}\nD:\t{seq_b}")
    cmd = ' '.join([c_exec_filename, f'-parallel {Hyperparams.PARALLEL_TYPE}', f'-threads {n_threads}', f'-path {Hyperparams.PATH_TO_MERGED}', f'-id D1 -match 1 -gap 0 -mismatch -1'])
    proc = os.popen(cmd).read().split("\n")
    os.system(f"rm -f {Hyperparams.PATH_TO_REPORT}")  # needed because of 'mkdir: cannot create directory ‘reports’: File exists' error
    return proc

    
def parse(out: list[str]) -> float:
    for line in out:
        if line.startswith("E)"):
            parsed = line.split(":")[-1].split("seconds")[0].strip()
            # todo: check if fail
            return float(parsed)


def dump_data_file(data_path: str, test_case_id_container: list[str], time_avg_container: list[float], time_stdev_container: list[float]):
    with open(data_path, "w") as ofile:
        for test_case_id, time_avg, time_stdev in zip(test_case_id_container, time_avg_container, time_stdev_container):
            ofile.write(f"TIME | CASE: {test_case_id}\t{time_avg}\t{time_stdev}\n")


def generate_g_data():
    create_dummy_files()
    os.system(f"mkdir {Hyperparams.PATH_TO_DATA}")
    file_container = os.listdir(Hyperparams.SEQ_FOLDER)    
    test_case_id_container = []
    time_avg_container = []
    time_stdev_container = []
    
    time_container = []

    for n_thread in range(1, Hyperparams.MAX_N_THREADS+1):
        print(f"doing thread {n_thread}")
        for file_idx in range(0, len(file_container), 2):   
            test_case_id = file_container[file_idx].split("_")[0]
            if int(test_case_id[-1]) > Hyperparams.MAX_TESTCASE_DIM:
                continue
            print(test_case_id)
            time_container.clear()
            for _ in range(Hyperparams.N_RUNS):
                impl_output = run_impl(
                    f"{Hyperparams.SEQ_FOLDER}/{file_container[file_idx]}", 
                    f"{Hyperparams.SEQ_FOLDER}/{file_container[file_idx+1]}",
                    n_thread)
                time_container.append(parse(impl_output))
            time_avg = avg(time_container)
            test_case_id_container.append(test_case_id)
            time_avg_container.append(time_avg)
            time_stdev_container.append(std_dev(time_container, time_avg))
        dump_data_file(
            f"{Hyperparams.PATH_TO_DATA}/data{n_thread}.txt", 
            test_case_id_container, 
            time_avg_container,
            time_stdev_container
        )
        test_case_id_container.clear()
        time_avg_container.clear()
        time_stdev_container.clear()

    cleanup_dummy_files()

    
def parse_data_file(data_path: str):
    test_id_container = []
    time_avg_container = []
    time_stdev_container = []

    with open(data_path, "r") as rfile:
        for line in rfile.readlines():
            test_id, time_avg, time_stdev = line.split("\t")
            test_id_container.append(test_id)
            time_avg_container.append(time_avg)
            time_stdev_container.append(time_stdev)
    return test_id_container, time_avg_container, time_stdev_container


def get_speedup_data():
    """converts 'per-test case' 2D array into 'per-threads' 2D array.
    FROM:
        [TEST CASE 11 with 1 THREAD , TEST CASE 12 with 1 THREAD , TEST CASE 13 with 1 THREAD , TEST CASE 14 with 1 THREAD , ...],
        [TEST CASE 21 with 2 THREADS, TEST CASE 22 with 2 THREADS, TEST CASE 23 with 2 THREADS, TEST CASE 24 with 2 THREADS, ...],
        ...
    TO: 
        [TEST CASE 11 with 1 THREAD , TEST CASE 11 with 2 THREAD , TEST CASE 11 with 3 THREAD , TEST CASE 11 with 4 THREAD , ...],
        [TEST CASE 12 with 1 THREAD , TEST CASE 12 with 2 THREAD , TEST CASE 12 with 3 THREAD , TEST CASE 12 with 4 THREAD , ...],
        ...
        """
    time_avg_container_2D = []
    time_stdev_container_2D = []
    for filename in os.listdir(Hyperparams.PATH_TO_DATA):
        _, time_avg_container, time_stdev_container = parse_data_file(f"{Hyperparams.PATH_TO_DATA}/{filename}")
        time_avg_container_2D.append(time_avg_container)
        time_stdev_container_2D.append(time_stdev_container)
    return transpose_array(time_avg_container_2D), transpose_array(time_stdev_container_2D)


def main():
    generate_g_data()
    

if __name__ == "__main__":
    main()