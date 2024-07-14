import sys
import logging
from Bio import SeqIO
import numpy as np

def setup_logger(log_file):
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger()
    return logger

def calculate_summary_statistics(fasta_file, output_file, log_file):
    logger = setup_logger(log_file)
    try:
        logger.info("Starting summary statistics calculation.")
        read_lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

        mean_read_length = np.mean(read_lengths)
        median_read_length = np.median(read_lengths)
        num_reads = len(read_lengths)
        read_length_n50 = calculate_n50(read_lengths)
        stdev_read_length = np.std(read_lengths)
        total_bases = np.sum(read_lengths)

        with open(output_file, 'w') as f:
            f.write(f"General summary:\n")
            f.write(f"Mean read length:\t{mean_read_length:.1f}\n")
            f.write(f"Median read length:\t{median_read_length:.1f}\n")
            f.write(f"Number of reads:\t{num_reads:.1f}\n")
            f.write(f"Read length N50:\t{read_length_n50:.1f}\n")
            f.write(f"STDEV read length:\t{stdev_read_length:.1f}\n")
            f.write(f"Total bases:\t\t{total_bases:.1f}\n")

        logger.info("Summary statistics calculation completed successfully.")
    except Exception as e:
        logger.error(f"Error occurred: {e}")
        raise

def calculate_n50(read_lengths):
    sorted_lengths = sorted(read_lengths, reverse=True)
    cumsum = np.cumsum(sorted_lengths)
    total_length = cumsum[-1]
    n50 = sorted_lengths[np.where(cumsum >= total_length * 0.5)[0][0]]
    return n50

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    output_file = sys.argv[2]
    log_file = sys.argv[3]
    calculate_summary_statistics(fasta_file, output_file, log_file)

