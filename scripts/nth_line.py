# Specify the path to your input and output files
input_file_path = 'scripts/Data/gH_20000/path_formatted.txt'
output_file_path = 'scripts/Data/gH1000/path.txt'
# Specify the interval of lines you want to write to the new file (e.g., every 5th line)
n = 20

# Open the input file in read mode and the output file in write mode
with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
    # Always write the first line (column headers or information) to the output file
    output_file.write(next(input_file))
    # Iterate through each line in the input file, along with its index
    # Using enumerate starting from 2 because the first line (headers) is already written
    for i, line in enumerate(input_file, start=2):
        # Check if the line is the nth line, adjusting for the fact that the first line is always included
        if (i-1) % n == 0:
            # Write the current line to the output file
            output_file.write(line)

print(f'The header and every {n}th line from "{
      input_file_path}" starting from the second line has been written to "{output_file_path}".')
