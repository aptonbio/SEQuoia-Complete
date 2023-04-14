# Reads a tab-delimited file with eight columns - 4 lines from the insert fastq
# followed by 4 lines from the UMI fastq - removes lines with invalid UMIs, adds
# the UMI to the insert record name, and writes out the insert reads as fastq.
from argparse import ArgumentParser
import csv
import sys
from typing import Tuple


def format_record(record, umi):
    return f"@{umi}_{record[0][1:]}\n{record[1]}\n{record[2]}\n{record[3]}\n"


def debarcode(
    input_tsv, output_fastq, excluded, expected_umi_length=8, max_umi_missing=1
) -> Tuple[int, int]:
    num_good = 0
    num_bad = 0
    for record in csv.reader(input_tsv, delimiter="\t"):
        umi = record[5].rstrip()
        umi_length = len(umi)
        if umi_length > expected_umi_length:
            raise Exception(
                f"UMI {umi} longer than expected length {expected_umi_length}"
            )
        # exclude any reads where the number of non-N bases is less than
        # expected_umi_length - max_umi_missing
        missing = expected_umi_length - umi_length
        i = 0
        while True:
            if missing > max_umi_missing:
                if excluded:
                    excluded.write(format_record(record, umi))
                num_bad += 1
                break
            elif i >= umi_length:
                output_fastq.write(format_record(record, umi))
                num_good += 1
                break
            else:
                if umi[i] == "N":
                    missing += 1
                i += 1

    return (num_good, num_bad)


def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input-file")
    parser.add_argument("-o", "--output-file", default=None)
    parser.add_argument("-x", "--excluded-file", default=None)
    parser.add_argument("-s", "--stats-file", default=None)
    parser.add_argument("--no-stats", action="store_true", default=False)
    parser.add_argument("--umi-length", type=int, default=8)
    parser.add_argument("--max-umi-missing", type=int, default=1)
    args = parser.parse_args()

    assert args.umi_length > 0
    assert args.max_umi_missing < args.umi_length

    close_output = args.output_file is not None
    output_fastq = open(args.output_file, "w") if close_output else sys.stdout
    excluded = open(args.excluded_file, "w") if args.excluded_file else None

    try:
        with open(args.input_file) as input_tsv:
            num_good, num_bad = debarcode(
                input_tsv, output_fastq, excluded, args.umi_length, args.max_umi_missing
            )
    finally:
        if close_output:
            output_fastq.close()
        if excluded:
            excluded.close()

    close_stats = args.stats_file is not None
    if not args.no_stats and (close_stats or close_output):
        stats = open(args.stats_file, "w") if close_stats else sys.stdout
        try:
            print(f"Total Reads: {num_good + num_bad}", file=stats)
            print(f"Good Reads: {num_good}", file=stats)
            print(f"Bad Reads: {num_bad}", file=stats)
        finally:
            if close_stats:
                stats.close()


if __name__ == "__main__":
    main()
