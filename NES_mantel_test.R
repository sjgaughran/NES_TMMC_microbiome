# Load required libraries
library(adegenet)
library(vcfR)
library(vegan)
library(dartR)

male_ids <- c("ES4721", "ES4716", "ES4676", "ES4709", "ES4731", "ES4681", "ES4682",
              "ES4707", "ES4685", "ES4689", "ES4730", "ES4711", "ES4691")

all_ids <- c("ES4721", "ES4716", "ES4676", "ES4709", "ES4731", "ES4681", "ES4682",
              "ES4707", "ES4685", "ES4689", "ES4730", "ES4711", "ES4691", "ES4708",
             "ES4695", "ES4675", "ES4679", "ES4672", "ES4686")


# Read VCF and create genlight object (calculate genetic distance once)
vcf <- read.vcfR("~/Downloads/seal_maf03_geno80_mind30_neutral_unlinked_chrEdited_renamed.filtered.recode.vcf")
gl <- vcfR2genlight(vcf)
indNames(gl)
gl <- gl[indNames(gl) %in% all_ids, ]
genetic_dist <- dist(as.matrix(gl))  # Calculate Euclidean genetic distance once

# Read Bray-Curtis distance matrix for T1
bray_curtis_matrix_t1 <- as.matrix(read.table("~/Downloads/bray-curtis-distance-matrix_2024_T1.tsv", header = TRUE, row.names = 1, sep = "\t"))
common_samples_t1 <- intersect(rownames(bray_curtis_matrix_t1), rownames(as.matrix(genetic_dist)))

# Subset genetic distance and Bray-Curtis matrices for common samples (T1)
genetic_dist_matrix <- as.matrix(genetic_dist)[common_samples_t1, common_samples_t1]
bray_curtis_matrix_subset_t1 <- bray_curtis_matrix_t1[common_samples_t1, common_samples_t1]

# Read Bray-Curtis distance matrix for T2
bray_curtis_matrix_t2 <- as.matrix(read.table("~/Downloads/bray-curtis-distance-matrix_2024_T2.tsv", header = TRUE, row.names = 1, sep = "\t"))
common_samples_t2 <- intersect(rownames(bray_curtis_matrix_t2), rownames(genetic_dist_matrix))

# Subset Bray-Curtis matrices for common samples (T2)
bray_curtis_matrix_subset_t2 <- bray_curtis_matrix_t2[common_samples_t2, common_samples_t2]

# Ensure the common samples are identical in order for both Bray-Curtis matrices and genetic distance
common_samples <- intersect(common_samples_t1, common_samples_t2)
genetic_dist_matrix <- genetic_dist_matrix[common_samples, common_samples]
bray_curtis_matrix_subset_t1 <- bray_curtis_matrix_subset_t1[common_samples, common_samples]
bray_curtis_matrix_subset_t2 <- bray_curtis_matrix_subset_t2[common_samples, common_samples]

# Calculate the change in Bray-Curtis distance between T1 and T2
bray_curtis_diff <- bray_curtis_matrix_subset_t2 - bray_curtis_matrix_subset_t1

# Convert genetic distance and Bray-Curtis difference matrices to vectors for plotting
genetic_dist_vector <- as.vector(genetic_dist_matrix[lower.tri(genetic_dist_matrix)])
bray_curtis_t1_vector <- as.vector(bray_curtis_matrix_subset_t1[lower.tri(bray_curtis_matrix_subset_t1)])
bray_curtis_t2_vector <- as.vector(bray_curtis_matrix_subset_t2[lower.tri(bray_curtis_matrix_subset_t2)])

bray_curtis_diff_vector <- as.vector(bray_curtis_diff[lower.tri(bray_curtis_diff)])

# Save the plot as a high-quality PNG
png("genetic_distance_vs_bray_curtis_t1.png", width = 6, height = 6, units = "in", res = 400)

# Plot genetic distance (X) vs. change in Bray-Curtis distance (Y)
plot(genetic_dist_vector, bray_curtis_t1_vector,
     xlab = "Genetic Distance (Euclidean)",
     ylab = "Bray-Curtis Dissimilarity",
     ylim = c(0.3, 1),  # Set Y-axis range
     pch = 19, col = "#FEA502")  # Change dot color to orange (#FEA502)
abline(lm(bray_curtis_t1_vector ~ genetic_dist_vector), col = "black", lty = 1)  # Add solid black trend line

# Finish saving the plot
dev.off()

# Optional: Calculate the correlation between genetic distance and change in Bray-Curtis
cor_result <- cor.test(genetic_dist_vector, bray_curtis_t1_vector, method = "spearman")
print(cor_result)


# Save the plot as a high-quality PNG
png("genetic_distance_vs_bray_curtis_t2.png", width = 6, height = 6, units = "in", res = 400)

# Plot genetic distance (X) vs. change in Bray-Curtis distance (Y)
plot(genetic_dist_vector, bray_curtis_t2_vector,
     xlab = "Genetic Distance (Euclidean)",
     ylab = "Bray-Curtis Dissimilarity",
     ylim = c(0.3, 1),  # Set Y-axis range
     pch = 19, col = "#002EFF")  # Change dot color to orange (#FEA502)
abline(lm(bray_curtis_t2_vector ~ genetic_dist_vector), col = "black", lty = 1)  # Add solid black trend line

# Finish saving the plot
dev.off()

# Optional: Calculate the correlation between genetic distance and change in Bray-Curtis
cor_result <- cor.test(genetic_dist_vector, bray_curtis_t2_vector, method = "spearman")
print(cor_result)


## ADMISSIONS ALL

# Read VCF and create genlight object (calculate genetic distance once)
vcf <- read.vcfR("~/Downloads/seal_maf03_geno80_mind30_neutral_unlinked_chrEdited_renamed.filtered.recode.vcf")
gl <- vcfR2genlight(vcf)
indNames(gl)
genetic_dist <- dist(as.matrix(gl))  # Calculate Euclidean genetic distance once

# Read Bray-Curtis distance matrix for T1
bray_curtis_matrix <- as.matrix(read.table("~/Downloads/bray-curtis-distance-matrix_2024_full.tsv", header = TRUE, row.names = 1, sep = "\t"))
common_samples <- intersect(rownames(bray_curtis_matrix), rownames(as.matrix(genetic_dist)))

# Subset genetic distance and Bray-Curtis matrices for common samples
genetic_dist_matrix <- as.matrix(genetic_dist)[common_samples, common_samples]
bray_curtis_matrix_subset <- bray_curtis_matrix[common_samples, common_samples]

# Convert genetic distance and Bray-Curtis difference matrices to vectors for plotting
genetic_dist_vector <- as.vector(genetic_dist_matrix[lower.tri(genetic_dist_matrix)])
bray_curtis_vector <- as.vector(bray_curtis_matrix_subset[lower.tri(bray_curtis_matrix_subset)])


# Plot genetic distance (X) vs. change in Bray-Curtis distance (Y)
png("genetic_distance_vs_bray_curtis_admitall.png", width = 6, height = 6, units = "in", res = 400)

plot(genetic_dist_vector, bray_curtis_vector,
     xlab = "Genetic Distance (Euclidean)",
     ylab = "Bray-Curtis Dissimilarity",
     ylim = c(0.3, 1),  # Set Y-axis range
     pch = 19, col = "#CD5301")  # Change dot color to orange (#FEA502)
abline(lm(bray_curtis_vector ~ genetic_dist_vector), col = "black", lty = 1)  # Add solid black trend line

# Finish saving the plot
dev.off()

# Optional: Calculate the correlation between genetic distance and change in Bray-Curtis
cor_result <- cor.test(genetic_dist_vector, bray_curtis_vector, method = "spearman")
print(cor_result)


## MALES

male_ids <- c("ES4663", "ES4682", "ES4689", "ES4691", "ES4705", "ES4707", "ES4711",
              "ES4719", "ES4728", "ES4730", "ES4731")
