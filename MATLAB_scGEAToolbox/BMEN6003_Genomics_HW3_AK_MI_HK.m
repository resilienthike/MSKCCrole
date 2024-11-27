%% COMP MODELLING BMEN-6003_HW3
%HARRY KABODHA
%AYMAN KHALEQ
%MARCUS IBRAHIM
%% Q0
[X,genelist,barcodelist] = sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv',1);
Metadata_Table = readtable('metadata.csv');
%% Q1-1
% We just read the workspace section on the right
%% Q1-2 part 1
patient_samples = unique(Metadata_Table.sample);
% To find out the number of unique patient samples
number_of_patients = length(patient_samples);
% Display the result
disp(['Number of patients in the dataset: ', num2str(number_of_patients)]);
%% Q1-2 part 2
% Sum the transcripts across all genes for each cell
total_transcripts_per_cell = sum(X, 1);
% Now, create the violin plot for cells
figure;
violinplot(total_transcripts_per_cell, Metadata_Table.sample);
title('Violin Plot: Total Transcripts per Cell Across Patient Samples');
% Count the number of genes with non-zero expression per cell
genes_per_cell = sum(X > 0, 1);
% Create the violin plot for genes
figure;
violinplot(genes_per_cell, Metadata_Table.sample);
title('Violin Plot: Total Genes per Cell Across Patient Samples');
%% Q1-3 normalization
X = sc_norm(X,'type','deseq');
%% Q2-1 plot UMAP
% Identify highly variable genes
[T, ~, genelist_sorted] = sc_hvg(X, genelist);
% Select the top 2000 over-dispersed genes
top_genes = genelist_sorted(1:2000);
top_X = X(1:2000, :);
% Generate UMAP embedding using the top 2000 over-dispersed genes
sc_plotcells(top_X, top_genes, "Top Genes", 3, false, true, 0.4, 15);
%sc_plotcells(X,genelist,"TP53")
%% Q2-2 Plot UMAPs in loop
% Parameters that we vary
min_dist_values = [0.1, 0.5, 0.79];
n_neighbors_values = [5, 20, 100];
% Loop over all combinations of min_dist and n_neighbors
for i = 1:length(min_dist_values)
    for j = 1:length(n_neighbors_values)
        % Set the current values for min_dist and n_neighbors
        min_dist = min_dist_values(i);
        n_neighbors = n_neighbors_values(j);
        
        % Create a figure for each combination
        figure;
        sc_plotcells(X, genelist, "none", 2, false, true, min_dist, n_neighbors);
        title(sprintf('UMAP with min_dist=%.2f, n_neighbors=%d', min_dist, n_neighbors));
    end
end
%% Q2-3
% First, calculate the UMAP coordinates using the normalized data for all genes
X_umap = sc_umap(X);
% Now, plot the expression of each marker gene
marker_genes = {'ITGAM', 'CD5', 'MOG'};
for i = 1:length(marker_genes)
    figure;
    sc_plotcells(X, genelist, marker_genes{i}, 2, false, true, 0.4, 15);
    title(['Expression of ', marker_genes{i}, ' across cells']);
end
% Lastly, plot the cell_assignment using the UMAP coordinates
figure;
gscatter(X_umap(:,1), X_umap(:,2), Metadata_Table.cell_assignment);
title('Distribution of cell types across UMAP');
	
%% Q3.1 Extract only malignant cells from the metadata table
Metadata_Malignant = Metadata_Table(strcmp(Metadata_Table.cell_assignment, 'Malignant'), :);
% Subset the expression matrix to include only malignant cells
X_malignant = X(:, strcmp(Metadata_Table.cell_assignment, 'Malignant'));
% Determine which of the malignant cells are from pediatric and adult tumors
is_pediatric = strcmp(Metadata_Malignant.GBM_type, 'Pediatric');
is_adult = strcmp(Metadata_Malignant.GBM_type, 'Adult');
% Subset the expression matrix to pediatric and adult tumors
X_pediatric = X_malignant(:, is_pediatric);
X_adult = X_malignant(:, is_adult);
% Perform differential expression analysis between pediatric and adult tumors
deg_results = sc_deg(X_pediatric, X_adult, genelist);
% Remove genes with Inf or NaN values in avg_log2FC and p_val_adj greater than the threshold
deg_filtered = deg_results(~isinf(deg_results.avg_log2FC) & ...
                           ~isnan(deg_results.avg_log2FC) & ...
                           deg_results.p_val_adj < 0.001, :);
% Sort the results by log2-fold change to find top genes for pediatric tumors
[~, idx_pediatric] = sort(deg_filtered.avg_log2FC, 'descend');
top_pediatric_genes = deg_filtered.gene(idx_pediatric(1:5));
% Sort the results by log2-fold change to find top genes for adult tumors
% Note that for adult tumors we need the negative log2-fold change
[~, idx_adult] = sort(deg_filtered.avg_log2FC, 'ascend');
top_adult_genes = deg_filtered.gene(idx_adult(1:5));
% Display the top 5 genes defining each tumor type
fprintf('Top 5 genes defining pediatric tumors:\n');
disp(top_pediatric_genes);
fprintf('Top 5 genes defining adult tumors:\n');
disp(top_adult_genes);
%% %Question 3.2
[T, ~, X_pediatric] = sc_hvg(X, genelist);
% Select the top 2000 over-dispersed genes
top_genes = X_pediatric(1:2000);
top_X = X(1:2000, :);
% Generate UMAP embedding using the top 2000 over-dispersed genes
sc_plotcells(top_X, top_genes, "TP53", 3, false, true, 0.4, 15);
sc_plotcells(X,genelist,"TP53")
set(gca,'Color','k')