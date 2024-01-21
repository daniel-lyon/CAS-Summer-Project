import pandas as pd

# Load the TSV file
df = pd.read_csv('params.tsv', delimiter='\t')

# Get the second column
second_column = df.iloc[:, 1]

# Print the second column
print(second_column)
