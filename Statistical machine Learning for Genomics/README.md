# Single-Cell Gene Expression Analysis Project
This repository contains an end-to-end workflow for analyzing single-cell gene expression data, specifically focusing on ~6K PBMCs (peripheral blood mononuclear cells) from a healthy donor. The data includes the 500 most variable genes and has been pre-filtered to ensure quality. The project applies state-of-the-art computational methods for normalization, dimensionality reduction, clustering, and differential gene expression analysis.

# Key Features
1. Data Preprocessing
Normalization: Scaled data to the median library size and applied a log transformation.
Dimensionality Reduction: Performed PCA followed by t-SNE or UMAP to visualize data in 2D or 3D space using the top 20 principal components.
2. Clustering
K-means Clustering: Applied to identify distinct cell clusters. The number of clusters (K) was carefully optimized and visualized on the PCA-reduced embedding.
Graph-Based Clustering: Built a k-nearest neighbor (kNN) graph and implemented the Louvain algorithm for community detection. Results were compared with K-means for robustness.
3. Differential Expression Analysis
Identified differentially expressed genes (DEGs) for specific clusters using statistical tests.
Visualized DEGs within the embedding and characterized gene expression patterns across clusters.
Explored alternative differential expression methods tailored to the data distribution.
4. Graph Analysis
Computed a 30-NN adjacency graph and visualized it as a heatmap to uncover underlying relationships between cells.
Justified the selection of distance metrics for constructing the graph.
# Requirements
Python 3.8 or higher
Key dependencies: numpy, pandas, scikit-learn, scanpy, matplotlib, and seaborn
For detailed installation instructions, see requirements.txt.
