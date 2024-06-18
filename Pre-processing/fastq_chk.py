def validate_fastq(file_path):
    try:
        with open(file_path, 'r') as file:
            line_number = 0
            while True:
                line_number += 1
                identifier = file.readline().strip()
                if not identifier:  # End of file
                    break
                if not identifier.startswith('@'):
                    print(f"Error at line {line_number}: Sequence identifier does not start with '@'.")
                    return False

                line_number += 1
                sequence = file.readline().strip()
                if not sequence:
                    print(f"Error at line {line_number}: Sequence line is missing.")
                    return False

                line_number += 1
                separator = file.readline().strip()
                if not separator.startswith('+'):
                    print(f"Error at line {line_number}: Separator line does not start with '+'.")
                    return False

                line_number += 1
                quality_scores = file.readline().strip()
                if not quality_scores:
                    print(f"Error at line {line_number}: Quality score line is missing.")
                    return False

                if len(sequence) != len(quality_scores):
                    print(f"Error at line {line_number}: Length of sequence and quality score line do not match.")
                    return False

        print("The FASTQ file is valid.")
        return True

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False

# Fixed file path
file_path = r'/content/sars_cov2_100.fastq'
validate_fastq(file_path)
