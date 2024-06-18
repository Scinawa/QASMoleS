import csv
from pysmiles import read_smiles
import networkx as nx


# Define the name of the CSV file
csv_filename = 'smiles_atomwise.txt'

# Open the CSV file for reading
with open(csv_filename, mode='r', newline='') as file:
    # Create a CSV reader object
    reader = csv.reader(file)
    
    # Iterate through each row in the CSV file
    for row in reader:
      #print(row)
      mol = read_smiles(row[0])
      adjmat = nx.to_numpy_array(mol)
      print(adjmat.shape)


