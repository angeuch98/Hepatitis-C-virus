# Hepatitis-C-virus

### Main goal
The main objective is to perform statistical inference and use of machine learning approaches in Hepatitis C virus data analysis. Machine learning techniques were used to build a binary classifier capable of classifying a new observation into a sick (Hepatitic C) or healthy group. The analysis was performed using R programming language.

### Dataset characteristic
The dataset was downloaded from the Kaggle platform. Donor is Institute of Clinical Chemistry; Medical University Hannover (MHH) in Germany. Dataset include 615 cases of laboratory values of blood donors and Hepatitis C patients and demographic values like age. The target attribute for classification is Category: Blood donors vs. Hepatitis C (including its progress: Hepatitis C, Fibrosis, Cirrhosis).

Attribute Information:
All attributes except Category and Sex are numerical (discrete and continous). The laboratory data are the attributes 5-14.
1) X (Patient ID)
2) Category (diagnosis) (values: '0=Blood Donor', '0s=suspect Blood Donor', '1=Hepatitis', '2=Fibrosis', '3=Cirrhosis')
3) Age (in years)
4) Sex (f,m)
5) ALB  - serum albumin level
6) ALP –  level of alkaline phosphatase level of alkaline phosphatase
7) ALT –  level of alkaline phosphatase
8) AST –  aspartate aminotransferase level; liver enzyme
9) BIL –  level of bilirubin in the blood (bile pigment)
10) CHE –  cholinesterase level
11) CHOL –  blood cholesterol level
12) CREA –  creatinine level
13) GGT –  level of gamma-glutamyltransferase (it is useful in the diagnosis of acute and chronic liver diseases)
14) PROT – total protein

### Methodology 
![image](https://github.com/angeuch98/Hepatitis-C-virus/assets/122879873/c0e9b8c1-256f-4724-9ff0-bc53a7767000)


Missing data were observed for 0.39% of all measurements. The largest number of missing values is for the  ALP and CHOL variable. Due to the fact that dataset is multidimensional performed missing data imputations using the MICE method (Multivariate Imputation by Chained Equations).
On the histograms it was possible to see that for ALT, BIL and GGT the distributions are very close to logarithmic ones, so it was decided to carry out the Cox box transformation. As a result, the distributions for these variables became close to the normal distribution. Because the distributions of variables are definitely not symmetrical, outlier detection was performed using the Hubert method for skewed distributions. **132 outliers** were obtained for this dataset. Performed kNN method to replace outliers.
In order to test the normality of variable distributions, Q-Q charts were drawn and statistical inference was performed based on the Shapiro-Wilk test. The null hypothesis for S-W assumes that the observations come from a normally distributed population. However, the alternative hypothesis is that the observations do not come from a normally distributed population. I will assume a significance level of 5%. If the value of the test probability determined on the basis of the test statistics is lower than the assumed level of significance, there is evidence that the null hypothesis has been rejected in favor of the alternative hypothesis.
Due to the fact that most of the variables have a non-normal distribution, the Levene's test was used. The null hypothesis assume that for two or more groups the observations come from a population with equal variances, in other words, the variances for the examined groups are the same. The alternative hypothesis, in turn, assumes that the variances are not equal.

For variables with:
- normal distribution
- independent samples
- enequal variances

performed Welch's test to chceck differences in mean between two groups (blood donor and hepatits).

For variables with:
- non-normal distribution
- independent samples
- unequal variances

performed Mann-Whitney U test to examine differences between two groups (blood donor and hepatits).
P-value correction was performed using the Bonferroni method, and then the effect size values were determined using the hedges coefficient (normally distributed features) and the rank biserial correlation coefficient (non-normally distributed features). The values of correlation coefficients between variables were also determined using the Pearson and Spearman methods. The classifier was built using logistic regression and the backward feature selection method.







