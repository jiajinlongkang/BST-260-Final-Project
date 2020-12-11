# BST260 Final Project

### Project Overview  
#### **Project title**: Predictive modeling for Insulin resistance using the NHANES

#### **Team Members**: Jack Kang, Peilu Wang, Rui Song, Sophie Kang

#### **Background**
Insulin resistance is an important component of metabolic syndrome and often precedes the development of type II diabetes. Mounting evidence suggests that insulin resistance could promote carcinogenesis and cancer progression (Wilcox, G., 2005. Insulin and insulin resistance. Clinical biochemist reviews, 26(2), p.19). The diagnosis of insulin resistance relies on invasive tests and the symptoms of insulin resistance could be easily overlooked (Wallace, T.M. and Matthews, D.R., 2002. The assessment of insulin resistance in man. Diabetic Medicine, 19(7), pp.527-534). We proposed to develop prediction models to predict insulin resistance, measured by Homeostatic model assessment (HOMA).

#### **Motivation**
We are interested in learning and developing prediction models for anthropometric measures or chronic diseases. We chose insulin resistance as our outcome because it is mechanistically associated with various diseases, including cardiometabolic diseases, liver diseases, and cancer development. Although the condition is reversible to some extent, its diagnosis largely relies on invasive testing. A parsimonious model with easy-to-access variables could be useful to predict insulin resistance in large scale epidemiological studies.

#### **Objectives**
**Project goals**: We aimed to develop a prediction model for HOMA among general US population and examine if demographic factors, lifestyle factors, anthropometric measures, and laboratory measures could yield an accurate prediction.

**Skills/features**: We will learn to (1) manipulate data sets with complex sampling schemes, (2) develop prediction models for both continuous and binary response variable, and (3) evaluate the performance of the prediction models.

#### **Data**
We used National Health and Nutrition Examination Survey (NHANES) 2009-2010 cycle to develop prediction models for Insulin Resistance among general US population.

**Inclusion criteria**: Adults (age>=18 years old) who participated in the NHANES 2009-2010 cycle.

**Exclusion criteria**: (1) Pregnant during the questionnaire cycle; (2) Previous diagnosis of diabetes or taking anti-diabetic medications (pills or insulin); (3) Missing the response variable (i.e. HOMA-IR) and/or predictors.

**Variables**: We considered demographic and lifestyle factors from questionnaires, physical examination measures, and laboratory measures.

### **Analysis**:

Step 1. Data wrangling & Exploratory data analysis

Step 2. Machine learning  
  * Binary HOMA-IR   
  * Continuous HOMA-IR   

Step 3. Shiny app development

Step 4. Result summary and presentation

<br/>

### **How to navigate the repository?** 
Our final Rmarkdown and HTML files are **'final report.Rmd'** and **'final report.html'**, and the shiny app is stored in the **'shiny_app'** folder with the file name **'app.R'**. 

/data_set: Raw datasets  

/shiny_app: Shiny app for this project. Run 'app.R' to execute the app.   

/pictures: Plots used to explain LASSO and AUC in our report.    

clean.csv: The cleaned dataset after data wrangling, preprocessing and filtering based on preset criteria. This is the dataset used for predictive model building.    

final report.Rmd: R markdown file for the project report.     

final report.html: HTML file for the project report.   

Our project website can be found [here](https://sites.google.com/view/insulinresistance). The screencase can be found [here](https://www.youtube.com/watch?v=IjeUkCC3yUQ).
