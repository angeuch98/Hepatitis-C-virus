# Hepatitis-C-virus

### Main goal
The main objective is to perform statistical inference and use of machine learning approaches in Hepatitis C virus data analysis. The analysis was performed using R programming language.

### Dataset characteristic
The dataset was downloaded from the Kaggle platform. Donor is Institute of Clinical Chemistry; Medical University Hannover (MHH); Hannover, Germany. Dataset include 615 cases of laboratory values of blood donors and Hepatitis C patients and demographic values like age. The target attribute for classification is Category: Blood donors vs. Hepatitis C (including its progress: Hepatitis C, Fibrosis, Cirrhosis).

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
![image info](./pictures/Methodology.png)

Missing data were observed for 0.39% of all measurements. The largest number of missing values is for the  ALP and CHOL variable. 
Due to the fact that dataset is multidimensional performed missing data imputations using the MICE method (Multivariate Imputation by Chained Equations) 
![image](https://github.com/angeuch98/Hepatitis-C-virus/assets/122879873/b6fbe717-912e-453f-9dd3-105cce8491ca)



