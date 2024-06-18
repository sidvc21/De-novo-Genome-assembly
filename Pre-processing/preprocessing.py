!pip install Biopython
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import Counter
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def calculate_gc_percentage(fastq_file):
    total_bases = 0
    gc_bases = 0

    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            total_bases += len(record.seq)
            gc_bases += record.seq.count("G") + record.seq.count("C")

    gc_percentage = (gc_bases / total_bases) * 100 if total_bases > 0 else 0
    return gc_percentage

def get_basic_statistics(fastq_file):
    # Initialize counters
    total_sequences = 0
    total_bases = 0
    poor_quality_sequences = 0
    sequence_lengths = Counter()
    gc_counts = Counter()

    quality_threshold = 28  # Minimum quality score for considering a base good quality
    n_percentage_threshold = 10  # Maximum allowed percentage of 'N' bases in a sequence

    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            total_sequences += 1
            total_bases += len(record.seq)

            # Check for poor quality sequences
            n_count = record.seq.count("N")
            if (n_count / len(record.seq)) * 100 > n_percentage_threshold or \
                    min(record.letter_annotations["phred_quality"]) < quality_threshold:
                poor_quality_sequences += 1
                continue  # Skip poor-quality sequences

            # Update sequence length and GC content counters
            sequence_lengths[len(record.seq)] += 1
            gc_counts[record.seq.count("G") + record.seq.count("C")] += 1

    return {
        "Filename": fastq_file,
        "File type": "FASTQ",
        "Encoding": "Illumina",  # Assuming Sanger encoding for simplicity
        "Total Sequences": total_sequences,
        "Total bases": total_bases,
        "Sequences flagged as poor quality": poor_quality_sequences,
        "Sequence length": dict(sequence_lengths)
    }

    print("Basic Statistics:")
    for key, value in basic_stats.items():
        print(f"{key}: {value}")

    print(f"The GC% in the FASTQ file is: {gc_percentage:.2f}%")

#PER BASE QUALITY
def calculate_per_base_quality(fastq_file):
    qualities = []
    max_sequence_length = 0

    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            sequence_length = len(record)
            max_sequence_length = max(max_sequence_length, sequence_length)
            quality_values = record.letter_annotations["phred_quality"]
            quality_values = quality_values[:sequence_length] + [0] * (max_sequence_length - len(quality_values))
            qualities.append(quality_values)

    qualities_array = np.array(qualities)
    per_base_quality = np.mean(qualities_array, axis=0)
    return per_base_quality


def plot_quality(per_base_quality):
    plt.plot(per_base_quality)
    plt.xlabel("Position in Read")
    plt.ylabel("Mean Quality Score")
    plt.title("Per Base Sequence Quality")
    plt.show()

#PER SEQUENCE QUALITY
def calculate_per_sequence_quality(fastq_file):
    sequence_qualities = []

    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            sequence_qualities.append(sum(record.letter_annotations["phred_quality"]) / len(record))

    return sequence_qualities

def plot_per_sequence_quality(sequence_qualities):
    sequence_indices = list(range(1, len(sequence_qualities) + 1))

    plt.bar(sequence_indices, sequence_qualities, color='blue')
    plt.xlabel("Sequence Index")
    plt.ylabel("Mean Sequence Quality Score")
    plt.title("Per Sequence Quality Scores")
    plt.show()

#PER BASE SEQUENCE CONTENT
def calculate_per_base_sequence_content(fastq_file):
    sequence_length = 0
    base_counts = {'A': [], 'T': [], 'G': [], 'C': []}

    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            sequence_length = len(record)
            if sequence_length == 0:
                continue

            for i, base in enumerate(record.seq):
                if base in base_counts:
                    base_counts[base].append(i)

    percentage_base_content = {base: np.histogram(counts, bins=sequence_length, range=(0, sequence_length))[0] / len(counts) * 100 for base, counts in base_counts.items()}
    return percentage_base_content

def plot_per_base_sequence_content(percentage_base_content):
    positions = range(1, len(percentage_base_content['A']) + 1)

    plt.plot(positions, percentage_base_content['A'], label='A', color='blue')
    plt.plot(positions, percentage_base_content['T'], label='T', color='green')
    plt.plot(positions, percentage_base_content['G'], label='G', color='orange')
    plt.plot(positions, percentage_base_content['C'], label='C', color='red')

    plt.xlabel("Position in Read")
    plt.ylabel("Percentage")
    plt.title("Per Base Sequence Content")
    plt.legend()
    plt.show()

#SEQUENCE LENGTH DISTRIBUTION
def plot_sequence_length_distribution(fastq_file):
    sequences = list(SeqIO.parse(fastq_file, "fastq"))
    lengths = [len(seq) for seq in sequences]

    plt.hist(lengths, bins=30, edgecolor="black")
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Number of Sequences")
    plt.title("Sequence Length Distribution")
    plt.show()

#STORE POOR QUALITY SEQUENCES
def get_poor_quality_sequences(fastq_file):
    poor_quality_sequences = []

    quality_threshold = 28  # Minimum quality score for considering a base good quality
    n_percentage_threshold = 10  # Maximum allowed percentage of 'N' bases in a sequence

    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            n_count = record.seq.count("N")
            if (n_count / len(record.seq)) * 100 > n_percentage_threshold or \
                    min(record.letter_annotations["phred_quality"]) < quality_threshold:
                poor_quality_sequences.append(record)

    return poor_quality_sequences

#ADAPTERS
def detect_adapters(fastq_file):
    adapter_sequences = {
        "Illumina Universal Adapter": "AGATCGGAAGAG",
        "Illumina Small RNA 3' Adapter": "TGGAATTCTCGG",
        "Illumina Small RNA 5' Adapter": "GATCGTCGGACT",
        "Nextera Transposase Sequence": "CTGTCTCTTATA",
        "PolyA": "A"*10,
        "PolyG": "G"*10
    }

    adapter_counts = {adapter: 0 for adapter in adapter_sequences}

    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            for adapter_name, adapter_sequence in adapter_sequences.items():
                if re.search(adapter_sequence, str(record.seq)):
                    adapter_counts[adapter_name] += 1

    return adapter_counts

def plot_adapter_counts(adapter_counts):
    adapters = list(adapter_counts.keys())
    counts = list(adapter_counts.values())

    plt.figure(figsize=(10, 6))
    plt.bar(adapters, counts, color=['blue', 'orange', 'green', 'red', 'purple', 'brown'])
    plt.title('Adapter Content QC')
    plt.xlabel('Adapter Sequences')
    plt.ylabel('Count')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

#REMOVE ADAPTERS
def remove_adapters(sequence, adapters):
    for adapter in adapters.values():
        if adapter in sequence:
            sequence = sequence.replace(adapter, '')
    return sequence

#HIGH QUALITY READS
def categorize_quality(quality_score):
    if quality_score >= 28:
        return "Best"
    elif 20 <= quality_score < 28:
        return "Medium"
    else:
        return "Poor"

def filter_and_write_high_quality_reads(input_fastq, output_high_quality):
    with open(input_fastq, "r") as input_handle, open(output_high_quality, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fastq"):
            min_quality = min(record.letter_annotations["phred_quality"])
            category = categorize_quality(min_quality)

            if category != "Poor":
                SeqIO.write(record, output_handle, "fastq")


#FINAL TRIM ORIGINAL FASTQ FILE TO TRIMMED FASTQ FILE
def trim_fastq(input_file, output_file, adapters):
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, 'fastq'):
            trimmed_sequence = remove_adapters(str(record.seq), adapters)
            trimmed_quality = record.letter_annotations["phred_quality"][:len(trimmed_sequence)]

            trimmed_record = SeqRecord(Seq(trimmed_sequence), id=record.id, name=record.name, description=record.description)
            trimmed_record.letter_annotations["phred_quality"] = trimmed_quality

            trimmed_record.description = f"length={len(trimmed_sequence)}"

            SeqIO.write(trimmed_record, out_handle, 'fastq')

#CONVERT TRIMMED FASTQ FILE TO FASTA FILE
def convert_fastq_to_fasta(fastq_file, fasta_file):
    with open(fasta_file, "w") as output_handle:
        for record in SeqIO.parse(fastq_file, "fastq"):
            SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    fastq_file = r'/content/sars_cov2_100.fastq'  # Replace with your FASTQ file path

    # 1. Basic Stats
    gc_percentage = calculate_gc_percentage(fastq_file)
    basic_stats = get_basic_statistics(fastq_file)

    print("Basic Statistics:")
    for key, value in basic_stats.items():
        print(f"{key}: {value}")

    print(f"The GC% in the FASTQ file is: {gc_percentage:.2f}%")

    # 2. Per Base Sequence Quality
    per_base_quality = calculate_per_base_quality(fastq_file)
    plot_quality(per_base_quality)

    # 3. Per Sequence Quality Scores
    per_sequence_qualities = calculate_per_sequence_quality(fastq_file)
    plot_per_sequence_quality(per_sequence_qualities)

    # 4. Per Base Sequence Content
    per_base_content = calculate_per_base_sequence_content(fastq_file)
    plot_per_base_sequence_content(per_base_content)

    # 5. Sequence Length Distribution
    plot_sequence_length_distribution(fastq_file)

    # 6. Poor Quality Sequences
    output_file = "poor_quality_sequences.txt"  # Output file to store poor quality sequences

    poor_quality_seqs = get_poor_quality_sequences(fastq_file)

    print(f"Total Poor Quality Sequences: {len(poor_quality_seqs)}")

    # Write poor quality sequences to a text file
    with open(output_file, "w") as out_handle:
        for idx, seq in enumerate(poor_quality_seqs):
            out_handle.write(f"Poor Quality Sequence {idx + 1}:\n")
            out_handle.write(f"ID: {seq.id}\n")
            out_handle.write(f"Sequence: {seq.seq}\n")
            out_handle.write(f"Quality Scores: {seq.letter_annotations['phred_quality']}\n\n")

    print(f"Poor quality sequences have been stored in '{output_file}'")

    # 7. Adapter Content
    adapter_counts = detect_adapters(fastq_file)
    print("Adapter Counts:")
    for adapter, count in adapter_counts.items():
        print(f"{adapter}: {count}")
    plot_adapter_counts(adapter_counts)

    #8 . High Quality Trimming

    input_fastq_file = r"/content/sars_cov2_100.fastq"  # Replace with your input FASTQ file path
    output_high_quality_file = r"/content/sars_cov2_100_high_quality_trimmed.fastq"

    filter_and_write_high_quality_reads(input_fastq_file, output_high_quality_file)


    # 8. Trimming
    input_fastq = r"/content/sars_cov2_100_high_quality_trimmed.fastq"
    output_fastq = r"/content/sars_cov2_100_final_trimmed.fastq"
    adapters = {
        "Illumina Universal Adapter": "AGATCGGAAGAG",
        "Illumina Small RNA 3' Adapter": "TGGAATTCTCGG",
        "Illumina Small RNA 5' Adapter": "GATCGTCGGACT",
        "Nextera Transposase Sequence": "CTGTCTCTTATA",
        "PolyA": "A"*10,
        "PolyG": "G"*10
    }
    trim_fastq(input_fastq, output_fastq, adapters)

    # 9. Fastq to Fasta Conversion
    fasta_file = r"/content/sars_cov2_100.fasta"
    convert_fastq_to_fasta(output_fastq, fasta_file)
    print(f"FASTQ file '{output_fastq}' converted to FASTA file '{fasta_file}'.")
