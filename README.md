# Bile-acid-Diabetes
BA DM Custom Code
*This file contains custom Stata (StataCorp, 2021, Release 17) and R (R-4.1.3 for Windows 64 bit) coding used for the analysis reported in the research article "Comprehensive clinical and genetic analyses of bile acids and their associations with diabetes and indices".
*Lines with asterisk(*) indicates comments and lines without * indicate Stata or R coding used.

*Manuscript Title: Comprehensive clinical and genetic analyses of bile acids and their associations with diabetes and indices

*Note: custom code was used to generate Tables 1 and S3 and Figures 1-4,S1 and S2
*Dependencies: tableone, rms, ggplot2, heatmap.2

*Author: Deepthi P. Mallela

*Contents:

*Part 1.  Data management
*Part 1.1 Importing variables
*Part 1.2 Label dataset
*Part 1.3 Adding labels to variable names
*Part 1.4 Defining outcome variables

*Part 2.  Missing data evaluation

*Part 3.  Descriptive statistics
*Part 3.1 Shapiro-Wilk normality test for variable distribution
*Part 3.2 Summary data for continuous variables - mean, standard deviation, median and interquartile range
*Part 3.3 k x 2 tables and Chi-square tests for categorical variables

*Part 4.  Tabulating mean, median, standard deviation, mininmum and maximum values for all bile acids by quartiles

*Part 5.  Calculation of Odds Ratios and 95% CI from logistic model 
*Part 5.1 Display of logit coefficients
*Part 5.2 Unadjusted and adjusted Odds Ratios and 95% CI among subjects not taking anti-diabetic medication 
*Part 5.3 Unadjusted and adjusted Odds Ratios and 95% CI for grouped plasma BA concentrations

*Part 6.  Fit null models, calculate LR test
*Part 7.  Hosmer-Lemeshow goodness-of-fit tests for regression models
*Part 8.  Generation of Table 1 (baseline characteristics of study participants)
*Part 9.  Forestplots for odds ratios generated for Q4 vs Q1 of bile acids with outcome variables.
*Part 10  Spearman correlation heatmap generated using heatmap.2 function in R.

*The below coding was performed using Stata (StataCorp, 2021, Release 17)

==================================================================================================
*Part 1.  Data management
==================================================================================================
*Part 1.1 Importing data to the statistical software (StataCorp, 2021, Release 17)
import delimited data.csv file 

*Part 1.2 Label dataset
label data "BA Diabetes"

*Part 1.3 Add labels to variable names
label variable age "Age"
label variable male "Sex, 1=male 0=female"
label variable smoking "Smoking, 1=yes 0=no"
label variable bpsystolic "Systolic Blood Pressure, 1=yes 0=no"
label variable bpdiastolic "Diastolic Blood Pressure, 1=yes 0=no"
label variable bmi "Body Mass Index"
label variable glucose "Glucose"
label variable insulin "Insulin"
label variable hba1c "HBA1C"
label variable hdlpriority "High Density Lipoprotein, 1=yes 0=no"
label variable ldlpriority "Low Density Lipoprotein, 1=yes 0=no"
label variable tgpriority "Triglyceride, 1=yes 0=no"
label variable aspirin "Aspirin, 1=yes 0=no"
label variable statin "Statin, 1=yes 0=no"
label variable crp "C-reactive protein"

*Part 1.4 Defining outcome variables (Homeostatic Model Assessment for Insulin Resistance (HOMA-IR), an index of insulin resistance using cutoff of ≥1.97), and risk for obesity (defined by body mass index (BMI) ≥30 kg/m2)
generate homair=1 if homa_ir>=1.97
replace homair=0 if homa_ir<1.97
generate obesity=1 if bmi>=30
replace obesity=0 if bmi<30

==================================================================================================
*Part 2. Missing data evaluation
==================================================================================================
*x=exposure or outcome variables
codebook x
summarize x, detail

==================================================================================================
*Part 3. Descriptive Statistics
==================================================================================================
*Bile acids with values less than LOD were replaced with 1/2 of minimum concentration of the analyte

*Part 3.1 Shapiro-Wilk normality test for variable distribution
*z=continuous variables
swilk z

*Part 3.2 Summary data for continuous variables - mean, standard deviation, median and interquartile range
*m=continuous variables (age, bpsystolic, bpdiastolic, bmi, hdlpriority, ldlpriority, tgpriority, glucose, insulin, hba1c)
*diabetics==0 (in subjects without diabetes)
*diabetics==1 (in subjects with diabetes)

summ m, detail
summ m if diabetics==0, detail 
summ m if diabetics==1, detail 

*Part 3.3 k x 2 tables for categorical variables
*n=categorical exposure variables(sex, smoking, aspirin, statin)

tabulate n diabetics, chi2 exact column
tabulate n diabetics, chi2 exact column

==================================================================================================================
*Part 4. Tabulating mean, median, standard deviation, mininmum and maximum values for all bile acids by quartiles
==================================================================================================================
*x=bile acid
*q=bile acid divided into 4 quartiles (Q1-Q4)
xtile q = x, nq(4)
tabstat x, stat(n mean median min max sd p50) by (q)

Example: xtile isolca = isolithocholicacid, nq(4)
tabstat isolithocholicacid, stat(n mean median min max sd p50) by (isolca)

=======================================================================================================
*Part 5. Calculation of Odds Ratios and 95% CI from logistic model - ANALYSES FOR FIGURES 2,3,4 AND S1
=======================================================================================================
*FIGURES 2 AND 3 ANALYSES

*y=diabetes, homair and obesity
*x=bile acid
*q1=bile acid divided into 4 quartiles (Q1-Q4) among total cohort (n=2,145)
xtile q1 = x, nq(4)
tabstat x, stat(n mean median min max sd p50) by (q1)

i.q=bile acid quartile (Q4 vs Q1)

*Unadjusted model
logistic y i.q1

*Adjusted model
logistic y i.q1 age male smoking bpsystolic hdlpriority ldlpriority tgpriority crp

*Part 5.1. Display of logit coefficients
logit y i.q1
logit y i.q1 age male smoking bpsystolic hdlpriority ldlpriority tgpriority crp

*Part 5.2. Unadjusted and adjusted Odds Ratios and 95% CI among subjects not taking anti-diabetic medication (DM meds)
*FIGURE S1 ANALYSES

*endm=1 if subjects are taking anti-diabetic medication
drop if endm==1

*y=diabetes, homair and obesity
*x=bile acid
*q2=bile acid divided into 4 quartiles (Q1-Q4) among subjects not on DM meds (n=1,885)
xtile q2 = x, nq(4)
tabstat x, stat(n mean median min max sd p50) by (q2)

*Unadjusted model
logistic y i.q2

*Adjusted model
logistic y i.q2 age male smoking bpsystolic hdlpriority ldlpriority tgpriority crp

*Part 5.3. Unadjusted and adjusted Odds Ratios and 95% CI for grouped plasma BA concentrations
*FIGURE S2 ANALYSES

*y=diabetes, homair and obesity
*x=bile acid groups (total,free,conj,conj_glycine,conj_taurine,total_primary,total_secondary,12-hydroxy,6-hydroxy,12-hydroxy/total,6-hydroxy/total)
*q3=bile acid divided into 4 quartiles (Q1-Q4) among total cohort (n=2,145)
xtile q3 = x, nq(4)
tabstat x, stat(n mean median min max sd p50) by (q3)

*q4=bile acid divided into 4 quartiles (Q1-Q4) among subjects not taking DM meds (n=1,885)
xtile q4 = x, nq(4)
tabstat x, stat(n mean median min max sd p50) by (q4)

*Unadjusted models
logistic y i.q3
logistic y i.q4

*Adjusted models
logistic y i.q3 age male smoking bpsystolic hdlpriority ldlpriority tgpriority crp
logistic y i.q4 age male smoking bpsystolic hdlpriority ldlpriority tgpriority crp

==================================================================================================
*Part 6.  Fit null models, calculate LR test
==================================================================================================
*Null model
*xi: logistic y i.x age male smoking bpsystolic hdlpriority ldlpriority tgpriority crp 
est store A

*z=adjustment variables (age male smoking bpsystolic hdlpriority ldlpriority tgpriority crp)
LR Test: z

For example,
*LR Test: age
quietly xi: logistic y i.x male smoking bpsystolic hdlpriority ldlpriority tgpriority crp 
lrtest A, stats

*LR Test: male
quietly xi: logistic y i.x age smoking bpsystolic hdlpriority ldlpriority tgpriority crp 
lrtest A, stats

==================================================================================================
*Part 7.  Hosmer-Lemeshow goodness-of-fit tests for regression models
==================================================================================================
estat gof
estat gof, group(10)
*note: obs collapsed on 10 quantiles of estimated probabilities.


*Below coding was performed using R software (R-4.1.3 for Windows 64 bit)

*Importing datafile.csv into f1 in R
rm(list = ls())
f1<-read.csv(file.choose(), header=TRUE)

==================================================================================================
*Part 8. Generation of Table 1 (baseline characteristics of study participants) - TABLE 1 ANALYSES
==================================================================================================
dim(f1)
install.packages("tableone")
library(tableone)
ls(f1)
myVars <- c("age", "male", "bpsystolic", "bpdiastolic", "smoking", "bmi", "ldlpriority", "hdlpriority", "tgpriority", "glucose", "insulin", "hba1c", "statin", "aspirin")
catVars <- c("diabetics")
tab1<-CreateTableOne(vars = myVars, strata ="diabetics",  data=f1, factorVars = catVars )
print(tab1, nonnormal =  myVars)
f2<- print(tab1, nonnormal =myVars, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.table(f2,"BA DM Table 1.csv",sep=",",col.names=NA)

*This same tableone R package was also used to generate Table S3 in the supplement by using the below additional code.
*BA includes all 36 bile acids measured in this study.
myVars <- c("BA")
catVars <- c("diabetics")
tab2<-CreateTableOne(vars = myVars, strata ="diabetics",  data=f1, factorVars = catVars )
print(tab2, nonnormal =  myVars)
f3<- print(tab2, nonnormal =myVars, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.table(f3,"BA DM Table S3.csv",sep=",",col.names=NA)

================================================================================================================================
*Part 9. Forestplots for odds ratios generated for Q4 vs Q1 of bile acids with outcome variables - FIGURES 2, 3 and S1 ANALYSES
================================================================================================================================
*outcome variables=diabetes, homair, and obesity.
*rms and ggplot2 R packages were used.
*Analyte=bile acid

library(rms)
library(ggplot2)
p <- ggplot(f1, aes(x=Analyte, y=HR, ymin=Lower, ymax=Upper)) +
geom_linerange(size=0.1, colour="black") +
geom_hline(aes( yintercept=1), lty=2) +
geom_point(size=4, shape=21, fill="Black", colour="Black", diabetes = 0.5) +
scale_x_discrete(name="Diabetes") +
scale_y_continuous(name="Odds ratio(95% CI)", limits = c(0, 3.0)) +
coord_flip() +
theme_classic()
p
library(forcats)
p + aes(x = fct_inorder(Analyte))
p + aes(x = fct_rev(Analyte))

=======================================================================================================
*Part 10. Spearman correlation heatmap generated using heatmap.2 function in R - FIGURE S2 ANALYSES
=======================================================================================================
rm(list = ls())
f1<-read.csv(file.choose(), header=TRUE)
f2<-f1[,-1]
rownames(f2)<-f1[,1]
f3 <- scale(f2)
heatmap(f3, scale = "none")
heatmap.2(f3, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")


=======================END OF CUSTOM CODE FILE==============================================








