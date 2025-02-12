import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Process icSHAPE and RNAplFold data.')
parser.add_argument('--pname', required=True, help='Protein name for identification')
parser.add_argument('--vitroShapeFolder', required=True, help='Folder containing in vitro icSHAPE files')
parser.add_argument('--vivoShapeFolder', required=True, help='Folder containing in vivo icSHAPE files')
parser.add_argument('--plFolder', required=True, help='Folder containing RNAplFold files')
parser.add_argument('--plshuffle', required=True, help='Folder containing shuffle files')
args = parser.parse_args()

# Assign arguments to variables
pname = args.pname
vitroShapeFolder = args.vitroShapeFolder
vivoShapeFolder = args.vivoShapeFolder
plFolder = args.plFolder
plshuffle = args.plshuffle

# Function to load RNAplFold data
def load_plfold(shape_path):
    shape = pd.read_csv(shape_path, skiprows=2, header=None, delim_whitespace=True)
    return shape

# Function to load icSHAPE data
def load_shape_file(shape_path):
    shape = pd.read_csv(shape_path, sep="\t", header=None)
    return shape

# Function to load shuffle data
def load_shuffle_file(shape_path):
    shape = pd.read_csv(shape_path, skiprows=2, header=None, delim_whitespace=True)
    return shape

# List files in the specified folders
pl_files = [os.path.join(plFolder, f) for f in os.listdir(plFolder)]
vitro_files = [os.path.join(vitroShapeFolder, f) for f in os.listdir(vitroShapeFolder)]
vivo_files = [os.path.join(vivoShapeFolder, f) for f in os.listdir(vivoShapeFolder)]
shuffle_files = [os.path.join(plshuffle, f) for f in os.listdir(plshuffle)]

# Load data from files
pl_tabs = [load_plfold(f) for f in pl_files]
vitro_shape_tabs = [load_shape_file(f) for f in vitro_files]
vivo_shape_tabs = [load_shape_file(f) for f in vivo_files]
shuffle_tabs = [load_plfold(f) for f in shuffle_files]

# Function to combine columns from multiple dataframes
def combine_cols(shape_dfs):
    combined = pd.concat([df.iloc[:, 1] for df in shape_dfs], axis=1)
    return combined

# Combine data from multiple files
pl_combine = combine_cols(pl_tabs)
vitro_shape = combine_cols(vitro_shape_tabs)
vivo_shape = combine_cols(vivo_shape_tabs)
shuffle_combine = combine_cols(shuffle_tabs)

# Convert combined data to numpy arrays
pl_matrix = pl_combine.to_numpy()
vitro_shape_matrix = vitro_shape.to_numpy()
vivo_shape_matrix = vivo_shape.to_numpy()
shuffle_matrix = shuffle_combine.to_numpy()

# Calculate median values for each dataset
pl_col_avg = np.median(pl_matrix, axis=1)
vitro_col_avg = np.median(vitro_shape_matrix, axis=1)
vivo_col_avg = np.median(vivo_shape_matrix, axis=1)
shuffle_col_avg = np.median(shuffle_matrix, axis=1)

# Generate a random background sample
np.random.seed(1)
background = np.random.choice(pl_col_avg, size=len(pl_col_avg), replace=False)

# Plot and save the median icSHAPE profile
plt.figure(figsize=(8, 8))
plt.plot(pl_col_avg, color='gray', linewidth=3, label='RNAplFold')
plt.plot(vitro_col_avg, color='blue', linewidth=3, label='vitro')
plt.plot(vivo_col_avg, color='orange', linewidth=3, label='vivo')
plt.plot(shuffle_col_avg, color='green', linewidth=3, label='shuffle')
plt.xlabel('Position')
plt.ylabel('Unpaired Probability')
plt.title(pname)
plt.legend(loc='upper left')
plt.savefig(f"{pname}_median_plot.pdf")
plt.close()
