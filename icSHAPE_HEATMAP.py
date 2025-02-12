import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.mixture import GaussianMixture
from argparse import ArgumentParser

# Parse command-line arguments
parser = ArgumentParser(description='Generate heatmaps and clustering analysis for icSHAPE data.')
parser.add_argument('--vitroShapeFolder', required=True, help='Folder containing in vitro icSHAPE files')
parser.add_argument('--vivoShapeFolder', required=True, help='Folder containing in vivo icSHAPE files')
args = parser.parse_args()

vitroShapeFolder = args.vitroShapeFolder
vivoShapeFolder = args.vivoShapeFolder

# Function to load icSHAPE data
def load_shape_file(shape_path):
    shape = pd.read_csv(shape_path, sep="\t", header=None)
    return shape

# List files in the specified folders
vitro_files = [os.path.join(vitroShapeFolder, f) for f in os.listdir(vitroShapeFolder) if f.endswith('.SHAPE')]
vivo_files = [os.path.join(vivoShapeFolder, f) for f in os.listdir(vivoShapeFolder) if f.endswith('.SHAPE')]

# Load data from files
vitro_shape_tabs = [load_shape_file(f) for f in vitro_files]
vivo_shape_tabs = [load_shape_file(f) for f in vivo_files]

# Function to combine columns from multiple dataframes
def combine_cols(shape_dfs):
    combined = pd.concat([df.iloc[:, 1] for df in shape_dfs], axis=1)
    return combined

# Combine data from multiple files
vitro_shape = combine_cols(vitro_shape_tabs)
vivo_shape = combine_cols(vivo_shape_tabs)

# Set row and column names
vitro_shape.columns = range(1, vitro_shape.shape[1] + 1)
vivo_shape.columns = range(1, vivo_shape.shape[1] + 1)
vitro_shape.index = [os.path.basename(f).replace('.SHAPE', '') for f in vitro_files]
vivo_shape.index = [os.path.basename(f).replace('.SHAPE', '') for f in vivo_files]

# Calculate differences between in vitro and in vivo shapes
common_indices = vitro_shape.index.intersection(vivo_shape.index)
vitro_shape_common = vitro_shape.loc[common_indices]
vivo_shape_common = vivo_shape.loc[common_indices]
shape_diffs = vitro_shape_common - vivo_shape_common

# Generate heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(shape_diffs, cmap='coolwarm', center=0)
plt.title('Heatmap of icSHAPE Differences')
plt.savefig('icSHAPE_heatmap.pdf')
plt.close()

# Elbow method for optimal number of clusters
wss = []
for i in range(1, 16):
    kmeans = KMeans(n_clusters=i, random_state=13)
    kmeans.fit(shape_diffs)
    wss.append(kmeans.inertia_)
plt.figure()
plt.plot(range(1, 16), wss, marker='o')
plt.xlabel('Number of clusters')
plt.ylabel('Within-cluster Sum of Squares')
plt.title('Elbow Method for Optimal Number of Clusters')
plt.savefig('elbow_method.pdf')
plt.close()

# Silhouette analysis
silhouette_scores = []
for i in range(2, 16):
    kmeans = KMeans(n_clusters=i, random_state=13)
    kmeans.fit(shape_diffs)
    silhouette_scores.append(silhouette_score(shape_diffs, kmeans.labels_))
plt.figure()
plt.plot(range(2, 16), silhouette_scores, marker='o')
plt.xlabel('Number of clusters')
plt.ylabel('Silhouette Score')
plt.title('Silhouette Analysis for Optimal Number of Clusters')
plt.savefig('silhouette_analysis.pdf')
plt.close()

# Hierarchical clustering
linked = linkage(shape_diffs, 'single')
plt.figure(figsize=(10, 7))
dendrogram(linked, labels=shape_diffs.index, distance_sort='descending', show_leaf_counts=True)
plt.title('Hierarchical Clustering Dendrogram')
plt.savefig('hierarchical_clustering.pdf')
plt.close()

# Gaussian Mixture Model clustering
gmm = GaussianMixture(n_components=10, random_state=13)
gmm.fit(shape_diffs)
labels = gmm.predict(shape_diffs)
plt.figure()
sns.scatterplot(x=shape_diffs.iloc[:, 0], y=shape_diffs.iloc[:, 1], hue=labels, palette='viridis')
plt.title('Gaussian Mixture Model Clustering')
plt.savefig('gmm_clustering.pdf')
plt.close()
