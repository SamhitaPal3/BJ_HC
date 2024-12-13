# BJ_HC

Samhita Pal and Xinge Jessie Jeng

This contains the R codes for siumulation and real data analysis for our project 'Discovering Candidate Genes Regulated by GWAS Signals in Cis and Trans'. This study introduces a novel approach to discover candidate genes regulated by GWAS signals in both cis and trans. We utilize adaptive statistical metrics (like Berk Jones test statistic and the Higher Criticism test statistic) that can reflect both the strong, sparse effects of cis-eQTLs and the weak, dense effects of trans-eQTLs, thereby helping us uncover joint effects that include functional impact of GWAS loci residing in non-coding regions (Jeng et al., 2013). 

### **Example: Evaluating Statistical Methods for SNP-Gene Relationships**

#### **Goal**
The simulation framework demonstrates the evaluation of various statistical tests (e.g., Higher Criticism (HC), Berk-Jones (BJ), Mean-based, and Minimum-P-value) in identifying SNPs associated with specific genes. It provides precision-recall analysis to assess the performance of these methods when distinguishing between active SNPs (signals) and non-active SNPs (noise).

#### **Input Data**
1. **SNP Data**:
   - The simulation uses SNP genotype data (e.g., `SNP_data_HapMap_ch21_n90.RData`) with 90 subjects and 2000 SNPs selected for analysis.

2. **Parameters**:
   - Number of genes (100).
   - Number of active SNPs per gene (20).
   - Signal strength for **strong cis-effects** (\(A_s = 4.25\)) and **weak trans-effects** (\(A_w = 0.2\)).
   - Variance of signals (\(\sigma_s = \sigma_w = 1.05\)).
   - Correlation structure of SNPs derived from the input genotype data.

3. **Statistical Tests**:
   - Four statistical measures are computed for each gene: 
     - **Berk-Jones (BJ)** (Berk and Jones, 1979)
     - **Higher Criticism (HC)** (Donoho and Jin, 2004)
     - **Mean-based**
     - **Minimum p-value**

#### **Methodology**
1. **Simulate Gene-SNP Effects**:
   - Strong cis-effects are generated for two SNPs near each gene.
   - Weak trans-effects are simulated for additional sets of SNPs (500 SNPs per set).
   - Noise is added to the SNP effects for non-active SNPs.

2. **Calculate Test Statistics**:
   - Compute \(z\)-scores for each SNP across genes, accounting for the correlation structure among SNPs.
   - Calculate the p-values and statistical measures for each gene using the `SetTest` package from Zhang et al. (2020).

3. **Rank Genes by Statistical Test**:
   - Each gene is ranked by the statistical measure (e.g., BJ, HC).
   - Genes containing active SNPs are expected to rank higher.

4. **Evaluate Performance**:
   - Compare the methods using **precision-recall (PR) curves** to assess how well each method identifies active SNPs.

#### **Output**
1. **Precision-Recall Curves**:
   - The PR curves visualize the trade-off between precision and recall for each method across all iterations of the simulation.
   - Example plot: 
     - **BJ** (Berk-Jones) and **HC** (Higher Criticism) methods tend to outperform the Mean-based and Minimum-p-value methods when signals are sparse.
   - Example visualization:

   ![PR Curves](PRcurve.png) 

2. **Summary Metrics**:
   - Mean precision and recall values across iterations for each statistical method.

#### **Code to Run the Example**
Here’s an example snippet to simulate SNP-gene associations and generate the PR curve:
```R
# Load SNP genotype data
load("SNP_data_HapMap_ch21_n90.RData")

# Run the simulation code provided
source("simulation1.R")

# Visualize Precision-Recall Curves
plot(BJ_mean_recall, BJ_mean_prec, 
     type = "l", col = "darkslateblue", ylab = "Precision", 
     xlab = "Recall", lwd = 2, main = "Precision-Recall Analysis")
lines(mean_mean_recall, mean_mean_prec, type = "l", col = "darkolivegreen4", lwd = 2, lty = 2)
lines(min_mean_recall, min_mean_prec, type = "l", col = "darkgoldenrod2", lwd = 2, lty = 6)
lines(HC_mean_recall, HC_mean_prec, type = "l", col = "deeppink3", lwd = 2, lty = 3)

legend("bottomleft", 
       legend = c("Mean-based", "BJ", "HC", "MinimumP"), 
       col = c("darkolivegreen4", "darkslateblue", "deeppink3", "darkgoldenrod2"),
       lty = c(2, 1, 3, 6), lwd = 2, cex = 0.7, horiz = FALSE, xpd = NA)
```

#### **Key Insights**
- **Berk-Jones (BJ)** and **Higher Criticism (HC)** outperform other methods in terms of precision and recall when signals are sparse and weak.

#### **Real-World Relevance**
This framework can be applied to analyze real SNP-gene data to evaluate the relative performance of different statistical measures for detecting weak signals in genomics.

---
### Reference

Berk, R. H. and D. H. Jones (1979). Goodness-of-fit test statistics that dominate the
kolmogorov statistics. Zeitschrift f¨ur Wahrscheinlichkeitstheorie und verwandte Gebiete 47 (1), 47–59.

Donoho, D. and J. Jin (2004). Higher criticism for detecting sparse heterogeneous mixtures.
Annals of Statistics 32 (3), 962–994.

Jeng, X. J., T. T. Cai, and H. Li (2013). Simultaneous discovery of rare and common segment
variants. Biometrika 100 (1), 157–172. PMCID:PMC3696347.

Zhang, H., J. Jin, and Z. Wu (2020). Distributions and power of optimal signal-detection
statistics in finite case. IEEE Transactions on Signal Processing 68, 1021–1033.
