import matplotlib.pyplot as plt
import numpy as np
import os


def extract_and_transform(input_file_path, output_file_path):
    # Initialize a list to hold all columns (each column is a list)
    columns = []

    with open(input_file_path, 'r') as file:
        current_column = []

        for line in file:
            # Strip leading/trailing whitespace
            stripped_line = line.strip()

            # Check if the line is empty (section divider)
            if stripped_line == '':
                if current_column:  # If the current column is not empty
                    columns.append(current_column)
                    current_column = []  # Reset for the next section
            else:
                # Split the line by whitespace and take the second value
                _, value = stripped_line.split()
                current_column.append(value)

        # Don't forget to add the last column if the file doesn't end with a blank line
        if current_column:
            columns.append(current_column)

    # Write the extracted data to the output file
    with open(output_file_path, 'w') as file:
        # Find the length of the longest column for iteration
        max_length = max(len(col) for col in columns)

        for i in range(max_length):
            # Extract the ith value from each column, if available
            line_values = [columns[j][i] if i < len(
                columns[j]) else '' for j in range(len(columns))]
            # Write the values as a comma-separated line
            file.write(', '.join(line_values) + '\n')


# Example usage
# Change this to your actual input file path
input_file_path = 'scripts/Data/QE_irrBZ/freq4x4x4.plot'
# Change this to your desired output file path
output_file_path = 'scripts/Data/QE_irrBZ/transformed_freq4x4x4.plot'

extract_and_transform(input_file_path, output_file_path)


# opens a matdyn.modes file (QE output) and retrieves the frequencies and eigenvectors
# Also performs energy conversion from THz to meV
def extract_data_grouped_by_q_vector(filename, skip_lines=0):

    # Conversion factor from THz to meV
    freq_to_energy = 4.136

    # Pattern to identify the line with q vector
    q_pattern = "q ="
    # Pattern to identify the line with frequency
    freq_pattern = "freq ("
    # Start with an empty list to collect the data
    data = []
    # Temporary storage for frequencies and vectors associated with the current q vector
    temp_frequencies_and_vectors = []

    with open(filename, 'r') as file:
        lines = file.readlines()

        for i, line in enumerate(lines):

            if i < skip_lines:
                continue
            
            if q_pattern in line:
                # Check if there's any collected data for the previous q vector
                if temp_frequencies_and_vectors:
                    # Append the collected data before resetting for the next q vector
                    data.append({
                        'q_vector': current_q_vector,
                        'entries': temp_frequencies_and_vectors
                    })
                    temp_frequencies_and_vectors = []  # Reset for the next q vector

                # Extract and update the current q vector
                q_vector_line = line.strip().split()
                current_q_vector = [float(q_vector_line[j]) for j in [2, 3, 4]] # indices anpassen nach datei!
                # convert units from 2pi/a to 1/a
                current_q_vector = [v * 2 * np.pi for v in current_q_vector]

            elif freq_pattern in line:
                # Corrected: Extract frequency in THz correctly by taking the value before [THz]
                freq_thz = float(line.split('=')[1].split('[')[0].strip())

                # The vector is on the next line, extract relevant components
                vector_line = lines[i+1].strip().replace('(',
                                                         '').replace(')', '').split()
                vector = [float(vector_line[j]) for j in [0, 2, 4]]

                vector = [v for v in vector]

                # Append frequency and vector to the temp storage
                temp_frequencies_and_vectors.append(
                    {'freq_thz': freq_thz * freq_to_energy, 'vector': vector})

        # Don't forget to add the last collected items after the loop
        if temp_frequencies_and_vectors:
            data.append({
                'q_vector': current_q_vector,
                'entries': temp_frequencies_and_vectors
            })

    return data



# opens a matdyn.modes file (QE output) and retrieves the frequencies and eigenvectors
# Also performs energy conversion from THz to meV
def extract_data_grouped_by_q_vector_dyn(filename):

    # Conversion factor from THz to meV
    freq_to_energy = 4.136

    # Pattern to identify the line with q vector
    q_pattern = "q ="
    # Pattern to identify the line with frequency
    freq_pattern = "freq ("
    # Start with an empty list to collect the data
    data = []
    # Temporary storage for frequencies and vectors associated with the current q vector
    frequencies_and_vectors = []

    q_vectors = []

    with open(filename, 'r') as file:
        lines = file.readlines()

        for i, line in enumerate(lines):
            
            if freq_pattern in line:
                # Corrected: Extract frequency in THz correctly by taking the value before [THz]
                freq_thz = float(line.split('=')[1].split('[')[0].strip())

                # The vector is on the next line, extract relevant components
                vector_line = lines[i+1].strip().replace('(',
                                                         '').replace(')', '').split()
                vector = [float(vector_line[j]) for j in [0, 2, 4]]

                vector = [v for v in vector]

                # Append frequency and vector to the temp storage
                frequencies_and_vectors.append(
                    {'freq_thz': freq_thz * freq_to_energy, 'vector': vector})

    with open(filename, 'r') as file:
        lines = file.readlines()

        for i, line in enumerate(lines):
            
            if q_pattern in line:
                
                # Extract and update the current q vector
                q_vector_line = line.strip().split()
                current_q_vector = [float(q_vector_line[j]) for j in [3, 4, 5]] # indices anpassen nach datei!
                # convert units from 2pi/a to 1/a
                current_q_vector = [v * 2 * np.pi for v in current_q_vector]
                q_vectors.append(current_q_vector)


    for q in q_vectors:
        data.append({'q_vector': q, 'entries': frequencies_and_vectors})

    return data

# Usage example:
# Replace 'data.txt' with the path to your actual file
#data = extract_data_grouped_by_q_vector('scripts/Data/QE_Pol_Disp/8x8x8_irrBZ.txt')

'''
for item in data:
    print(f"Q Vector: {item['q_vector']}")
    for entry in item['entries']:
        print(f"  Frequency: {entry['freq_thz']} THz, Vector: {entry['vector']}")
'''


def write_data_to_file(data, filename):
    header = "qx qy qz freq1 freq2 freq3 v1x v1y v1z v2x v2y v2z v3x v3y v3z\n"

    with open(filename, 'w') as file:
        # Write the header first
        file.write(header)

        for item in data:
            # Start with the q vector
            line_parts = [' '.join(map(str, item['q_vector']))]

            # Add the frequencies
            freqs = [str(entry['freq_thz']) for entry in item['entries']]
            line_parts.extend(freqs)

            # Add the vectors, concatenating components of each vector
            vectors = [' '.join(map(str, entry['vector']))
                       for entry in item['entries']]
            line_parts.extend(vectors)

            # Combine all parts into a single line
            line = ' '.join(line_parts)
            file.write(line + "\n")


# DIES HERE
#data = extract_data_grouped_by_q_vector('scripts/Data/dyn20x20x20')
#write_data_to_file(data, filename='Parameters/20x20x20.txt')

data = extract_data_grouped_by_q_vector('Outputs/full_path_new/bccfe.eig')
write_data_to_file(data, filename='Outputs/full_path_new/full_path_new_disp.txt')


def extract_data_from_dir(directory_path,outfile):
    # Initialize a list to hold the extracted data
    all_data = []

    of = open(outfile, 'w')

    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)
        if os.path.isfile(file_path):
            print(f"Processing {filename}...")
            data = extract_data_grouped_by_q_vector_dyn(file_path)

            print(data)

            for item in data:
                # Start with the q vector
                line_parts = [' '.join(map(str, item['q_vector']))]

                # Add the frequencies
                freqs = [str(entry['freq_thz']) for entry in item['entries']]
                line_parts.extend(freqs)

                # Add the vectors, concatenating components of each vector
                vectors = [' '.join(map(str, entry['vector']))
                        for entry in item['entries']]
                line_parts.extend(vectors)

                # Combine all parts into a single line
                line = ' '.join(line_parts)
                of.write(line + "\n")

            all_data.extend(data)
    




dir = 'scripts/Data/dyn20x20x20'
out_file = 'Parameters/20x20x20_disp.txt'
data = extract_data_from_dir(dir,out_file)


def plot_disp_from_file(filename):
    # Lists to store the q vectors and frequencies
    q_vectors = []
    frequencies = []

    with open(filename, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            parts = line.strip().split()
            # Extract q vector and frequencies
            q_vector = [float(parts[i]) for i in range(3)]
            freqs = [float(parts[i]) for i in range(3, 6)]

            q_vectors.append(q_vector)
            frequencies.append(freqs)

    # Convert to a format suitable for plotting
    # Transpose to have each frequency in its own list
    frequencies = list(zip(*frequencies))

    # Plotting
    plt.figure(figsize=(10, 6))
    for i, freqs in enumerate(frequencies):
        plt.plot(freqs, label=f'Frequency {i+1}')

    plt.xlabel('Index of q vector')
    plt.ylabel('disp (meV)')
    plt.legend()
    plt.show()


# Example usage:
# Replace 'formatted_data.txt' with the path to your actual formatted output file
#plot_disp_from_file('scripts/Data/QE_Pol_Disp/formatted_4x4x4_path_G_H.txt')
