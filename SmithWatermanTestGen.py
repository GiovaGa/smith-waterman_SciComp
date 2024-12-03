import sys
import os
import subprocess
import random
from Bio import pairwise2

from random import randint


class Hyperparams:
    SCORE_MATCH = 4
    SCORE_MISMATCH = -2
    SCORE_GAP_OPEN = -3
    SCORE_GAP_EXT = -1

    # TEST_ID_<SEQ_LEN> (=that is, the second digit of TEST_ID) will be inferred from index.
    # Update TEST_SEQ_LEN to the desired dimensions
    TEST_SEQ_LEN = (50, 100, 500, 1_000, 5_000, 10_000)

    # BEFORE RUNNING YOU NEED TO CREATE THIS FOLDERS
    TEST_FOLDER = r"./test"
    INPUT_FOLDER = f"{TEST_FOLDER}/input/"
    RAW_OUTPUT_FOLDER = f"{TEST_FOLDER}/output/"

    SCORE_DUMP_FULL_PATH = rf"{TEST_FOLDER}/score_dump.txt"
    COMPRESSED_OUTPUT_PATH = f"./archive.tar.gz"


class Globals:
    NUCLEOTIDES = ["A", "C", "G", "T"]
    NUCLEOTIDES_TO_IDX = {n: i for i, n in enumerate(NUCLEOTIDES)}
    IDX_TO_NUCLEOTIDE = {i: n for i, n in enumerate(NUCLEOTIDES)}


def write_seq(relative_path: str, seq: str, seq_len: int) -> None:
    """create and write a file containing only the given DNA sequence"""
    with open(relative_path, "w") as file:
        file.write(f">chr{random.randint(1, 23)}\t0\t{seq_len}\n")
        file.write(seq)


def read_seq(path_to_file: str) -> str:
    """read DNA sequence from fasta file"""
    with open(path_to_file, "r") as file:
        return file.readlines()[1].replace("\n", "").strip()


def write_expected_output(path_to_seqA: str, path_to_seqB: str, output_path: str) -> None:
    seqA = read_seq(path_to_seqA)
    seqB = read_seq(path_to_seqB)
    print(f"aligning {path_to_seqA} -> {path_to_seqB}")
    alignments = pairwise2.align.localms(seqA, seqB, Hyperparams.SCORE_MATCH, Hyperparams.SCORE_MISMATCH,
                                         Hyperparams.SCORE_GAP_OPEN, Hyperparams.SCORE_GAP_EXT)
    with open(f"{output_path}", "w") as of:
        of.write(f"{alignments[0].score}\n")


def _write_expected_output(path_to_seqA: str, path_to_seqB: str, output_path: str) -> None:
    # deprecated
    """create and write a file containing a string with the expected output
    For now it will just be the alignment scores, each on a new line"""

    proc = subprocess.Popen(f"./ssw_test -m {Hyperparams.SCORE_MATCH} -x {abs(Hyperparams.SCORE_MISMATCH)}"
                            f" -o {abs(Hyperparams.SCORE_GAP_OPEN)} -e {abs(Hyperparams.SCORE_GAP_EXT)} {path_to_seqA}"
                            f" {path_to_seqB}", stdout=subprocess.PIPE,
                            shell=True)
    (out, _) = proc.communicate()
    out = str(out)
    score_identifier = "_score:"
    temp: str = out.split(score_identifier)[1].strip()
    i = 0
    while temp[i].isdigit():
        i += 1
    with open(f"{output_path}", "w") as of:
        of.write(f"{temp[:i]}\n")


def insert_at_idx(a: list, idx: int, cargo: list) -> list:
    # TODO: fix because it's a sloppy solution
    return a[:idx] + cargo + a[idx:]


def get_fixed_length_random_fragment(seq: str, frag_length: int) -> str:
    start_pos: int = randint(0, len(seq) - frag_length)
    return seq[start_pos:start_pos + frag_length]


def remove_N(a: str) -> str:
    """removes undefined nucleotide N from any given seq"""
    return a.replace("N", random.choice(Globals.NUCLEOTIDES))


def mismatch_n_nucleotides(s: str, n: int) -> str:
    """assure that exactly N nucleotides get substituted (no nucleotide can be substituted twice)"""
    seq_len: int = len(s)
    if n > seq_len:
        print(f"[ERROR]: The number of nucleotides to be replaced [given {n}] cannot exceed {len(s)}. Exiting...")
        sys.exit(1)
    seq: list[str] = list(s)
    idxs = random.sample(range(0, seq_len), n)
    for i in idxs:
        # TODO: refactor this crap
        if seq[i] == "A":
            seq[i] = random.choice(["C", "G", "T"])
        elif seq[i] == "C":
            seq[i] = random.choice(["A", "G", "T"])
        elif seq[i] == "G":
            seq[i] = random.choice(["A", "C", "T"])
        elif seq[i] == "T":
            seq[i] = random.choice(["A", "C", "G"])
    return "".join(seq)


def insert_rr_gaps(s: str, max_nof_gap: int = 0, gap_lengths=None) -> str:
    """inserts a RANDOM amount of gaps (bounded from 1 to len(s)-1, inclusive) each of RANDOM length.
    The 'ATCGA' seq has 4 gap opening idxs marked with | -> 'A|T|C|G|A'.
    The gap opening idx refers to the char BEFORE which the gap is opened. With seq='ATCGA' and gap_opening_idx=1 the gap will open at | -> AT|CGA
    A single gap size is bound to the 40% of the total seq length.
    max_nof_gap allows to
    Should be self opening proof (that is no gap can open from a previously opened gap)
    """
    seq_len: int = len(s)
    seq = list(s)
    if max_nof_gap > seq_len - 1:
        print(
            f"[ERROR]: The maximum number of gaps [given {max_nof_gap = }] exceeds the length of the seq [given {seq_len = }].. Exiting...")
        sys.exit(1)
    nof_gaps = max_nof_gap if max_nof_gap > 0 else random.randint(1, seq_len - 1)
    max_single_gap_length: int = round(seq_len * 0.4)

    idxs = sorted(random.sample(range(0, seq_len - 1), nof_gaps))
    if gap_lengths is None:
        gap_lengths = [random.randint(1, max_single_gap_length) for _ in range(nof_gaps)]

    acc_len_delta = 0
    for idx, gap_opening_idx in enumerate(idxs):
        # generate random gap seq
        gap_length = gap_lengths[idx]
        gap_nucleotides = [random.choice(Globals.NUCLEOTIDES) for _ in range(gap_length)]

        gap_opening_idx += acc_len_delta
        seq = insert_at_idx(seq, gap_opening_idx, gap_nucleotides)
        acc_len_delta += gap_length
    return "".join(seq)


def insert_nr_gaps(s: str, n: int) -> str:
    """insert n gaps of RANDOM length. See insert_rr_gaps for more details"""
    return insert_rr_gaps(s, n)


def insert_nm_gaps(s: str, n: int, m: list[int]) -> str:
    """insert n gaps of m length.
    calling insert_nm_gaps('AAAAAA', 2, [3, 2]) will create a new sequence with 2 (because n=2) new gaps (inserted at
     random positions) of length 3 (the first one) and 2 (the second one):
    See insert_rr_gaps for more details"""
    if n > len(m):
        print(f"[ERROR]: too few gap lengths provided ({n=} while {len(m)=}. Exiting...")
    return insert_rr_gaps(s, n, m)


def insert_custom_gap_at_idx(s: str, i: int, gap_nucleotides: list[str]) -> str:
    """custom_idx_insertion and multiple_custom_gap_insertion are easily implementable"""
    seq = list(s)
    return "".join(insert_at_idx(seq, i, gap_nucleotides))


class TestCase:
    TEST_ID_PERFECT_MATCH = "1"
    TEST_ID_MULTIPLE_MISMATCH = "2"
    TEST_ID_SINGLE_GAP = "3"
    TEST_ID_MULTIPLE_GAPS = "4"

    # TODO: refactor ALL_* to common builder
    @staticmethod
    def ALL_perfect_match(src_seq: str):
        for idx, seq_len in enumerate(Hyperparams.TEST_SEQ_LEN):
            s = get_fixed_length_random_fragment(src_seq, seq_len)

            path_to_seqA: str = f"{Hyperparams.INPUT_FOLDER}/{TestCase.TEST_ID_PERFECT_MATCH}{idx + 1}_seqA.fa"
            path_to_seqB: str = path_to_seqA.replace("seqA", "seqB")

            write_seq(path_to_seqA, s, seq_len)
            write_seq(path_to_seqB, s, seq_len)

            write_expected_output(path_to_seqA, path_to_seqB,
                                  f"{Hyperparams.RAW_OUTPUT_FOLDER}/{TestCase.TEST_ID_PERFECT_MATCH}{idx + 1}.txt")

    @staticmethod
    def ALL_multiple_mismatch(src_seq: str):
        """we will mismatch 30% of the original nucleotides"""
        mismatch_rate = 0.3
        for idx, seq_len in enumerate(Hyperparams.TEST_SEQ_LEN):
            seqA = get_fixed_length_random_fragment(src_seq, seq_len)
            path_to_seqA = f"{Hyperparams.INPUT_FOLDER}/{TestCase.TEST_ID_MULTIPLE_MISMATCH}{idx + 1}_seqA.fa"
            write_seq(path_to_seqA, seqA, seq_len)

            seqB = mismatch_n_nucleotides(seqA, round(seq_len * mismatch_rate))
            path_to_seqB: str = path_to_seqA.replace("seqA", "seqB")
            write_seq(path_to_seqB, seqB, seq_len)

            write_expected_output(path_to_seqA, path_to_seqB,
                                  f"{Hyperparams.RAW_OUTPUT_FOLDER}/{TestCase.TEST_ID_MULTIPLE_MISMATCH}{idx + 1}.txt")

    @staticmethod
    def ALL_single_gap(src_seq: str):
        """we will open a gap with a length of 20% of the original length"""
        gap_len_factor = 0.2
        for idx, seq_len in enumerate(Hyperparams.TEST_SEQ_LEN):
            seqA = get_fixed_length_random_fragment(src_seq, seq_len)
            path_to_seqA = f"{Hyperparams.INPUT_FOLDER}/{TestCase.TEST_ID_SINGLE_GAP}{idx + 1}_seqA.fa"
            write_seq(path_to_seqA, seqA, seq_len)

            path_to_seqB = path_to_seqA.replace("seqA", "seqB")
            seqB = insert_rr_gaps(seqA, 1, [round(seq_len * gap_len_factor)])
            write_seq(path_to_seqB, seqB, seq_len)

            write_expected_output(path_to_seqA, path_to_seqB,
                                  f"{Hyperparams.RAW_OUTPUT_FOLDER}/{TestCase.TEST_ID_SINGLE_GAP}{idx + 1}.txt")

    @staticmethod
    def ALL_multiple_gaps(src_seq: str):
        """we will open 1 gap per 50 nucleotides with a length of 50 nucleotides"""
        nof_gaps_factor = 0.02
        max_gap_len = 50
        for idx, seq_len in enumerate(Hyperparams.TEST_SEQ_LEN):
            seqA = get_fixed_length_random_fragment(src_seq, seq_len)
            path_to_seqA = f"{Hyperparams.INPUT_FOLDER}/{TestCase.TEST_ID_MULTIPLE_GAPS}{idx + 1}_seqA.fa"
            write_seq(path_to_seqA, seqA, seq_len)

            path_to_seqB = path_to_seqA.replace("seqA", "seqB")
            max_nof_gaps = round(seq_len * nof_gaps_factor)
            seqB = insert_rr_gaps(seqA, max_nof_gaps, [max_gap_len] * max_nof_gaps)
            write_seq(path_to_seqB, seqB, seq_len)

            write_expected_output(path_to_seqA, path_to_seqB,
                                  f"{Hyperparams.RAW_OUTPUT_FOLDER}/{TestCase.TEST_ID_MULTIPLE_GAPS}{idx + 1}.txt")


def generate_test_cases():
    db_seq = read_seq("./1M.fa")
    parsed_str = remove_N(db_seq)
    TestCase.ALL_perfect_match(parsed_str)
    print("DONE PERF MATCH")
    TestCase.ALL_multiple_mismatch(parsed_str)
    print("DONE MULTIPLE MISMATCH")
    TestCase.ALL_single_gap(parsed_str)
    print("DONE SINGLE GAP")
    TestCase.ALL_multiple_gaps(parsed_str)
    print("DONE MULTIPLE GAPS")
    subprocess.Popen(f"tar -zcvf {Hyperparams.COMPRESSED_OUTPUT_PATH} {Hyperparams.TEST_FOLDER}")
    print("DONE COMPRESSING")


def dump_score_params_to_file(file_full_path: str):
    with open(file_full_path, "w") as file:
        file.write(f"-m {Hyperparams.SCORE_MATCH} -x {abs(Hyperparams.SCORE_MISMATCH)} -o {abs(Hyperparams.SCORE_GAP_OPEN)} -e {abs(Hyperparams.SCORE_GAP_EXT)}")


def main():
    # TODO: maybe this folders (and all the sequence files) can be deleted after compressing
    if not os.path.isdir(Hyperparams.TEST_FOLDER):
        print(f"[ERROR]: {Hyperparams.TEST_FOLDER} does not exist. Exiting...")
        sys.exit(1)
    if not os.path.isdir(Hyperparams.INPUT_FOLDER):
        print(f"[ERROR]: {Hyperparams.INPUT_FOLDER} does not exist. Exiting...")
        sys.exit(1)
    if not os.path.isdir(Hyperparams.RAW_OUTPUT_FOLDER):
        print(f"[ERROR]: {Hyperparams.RAW_OUTPUT_FOLDER} does not exist. Exiting...")
        sys.exit(1)
    generate_test_cases()
    dump_score_params_to_file(Hyperparams.SCORE_DUMP_FULL_PATH)


if __name__ == '__main__':
    main()
