### **Predictive Modeling For Insulin Resistance Using the NHANES Data**

**Team member (alphabetical order):** Jack Kang, Peilu Wang, Rui Song, Sophie Kang

**Assigned TA:** Santiago Romero-Brufau

#### **Background**
Insulin resistance is an important component of metabolic syndrome and often precedes the development of type II diabetes. Mounting evidence suggests that insulin resistance could promote carcinogenesis and cancer progression. The diagnosis of insulin resistance relies on invasive tests and the symptoms of insulin resistance could be easily overlooked. We proposed to develop prediction models to predict insulin resistance, measured by [**Homeostatic model assessment (HOMA)**](https://en.wikipedia.org/wiki/Homeostatic_model_assessment).

#### **Motivation**
We are interested in learning and developing prediction models for anthropometric measures or chronic diseases. We chose insulin resistance as our outcome because it is mechanistically associated with various diseases, including cardiometabolic diseases, liver diseases, and cancer development. Although the condition is reversible to some extent, its diagnosis largely relies on invasive testing. A parsimonious model with easy-to-access variables could be useful to predict insulin resistance in large scale epidemiological studies. 

#### **Objectives**
**Project goals:** We aimed to develop a prediction model for HOMA among general US population and examine if demographic factors, lifestyle factors, anthropometric measures, and laboratory measures could yield an accurate prediction. 

**Skills/features:** We will learn to (1) manipulate data sets with complex sampling schemes, (2) develop prediction models for both continuous and binary response variable, and (3) evaluate the performance of the prediction models.

#### **Data**
We used National Health and Nutrition Examination Survey (NHANES) 2009-2010 cycle to develop prediction models for Insulin Resistance among general US population.

**Inclusion criteria:** Adults (age>=18 years old) who participated in the NHANES 2009-2010 cycle.

**Exclusion criteria:** (1) Pregnant during the questionnaire cycle; (2) Previous diagnosis of diabetes or taking anti-diabetic medications (pills or insulin); (3) Missing the response variable (i.e. HOMA-IR).

**Variables:**
We considered demographic and lifestyle factors from questionnaires, physical examination measures, and laboratory measures.

#### **General Approach**
The dataset was randomly split into training set (80%) and test set (20%). 

An L1 penalized binary lasso regression was conducted on multiple potential predictors to create a sparse model in order to improve the direct interpretability as well as to avoid the collinearity and overfitting of the data. To tune the regularization parameter lambda, a ten-fold cross-validation was conducted, in which the training dataset was split into ten uniformly sized chunks. One of the ten chunks was used for validation and the rest was used for training for each time. We evaluated the performance of the binary lasso model in the test dataset through accuracy, sensitivity, specificity, and ROC. 

Then, a k-nearest neighbors (KNN) model was conducted for the binary outcome. The parameter k was selected through ten-fold validation. Same metrics were used to evaluate the performance of the binary KNN model in the test set. Moreover, lasso and KNN models were also conducted for continuous HOMA-IR. RMSE in test set was calculated for the evaluation of continuous models.

#### **Conclusion**
We developed several prediction models for the binary and continuous HOMA-IR in the NHANES dataset with 1121 participants and 37 predictors. For the classification of the binary HOMA-IR, the lasso model performed better than the KNN model, with higher accuracy and a more balanced sensitivity and specificity. It also achieved a higher AUC than KNN. Compared with the simple logistic model, the lasso model only showed slightly better performance. However, given that the lasso was more parsimonious and that a simpler model generalized better in the external dataset by avoiding overfitting, we prefer the lasso model in predicting binary HOMA-IR. For the continuous model, the RMSE was relatively large for both the lasso and the knn model. Since the classification of insulin resistance condition could be a more relevant clinical question that the general population would care about, we suggest implementing the binary lasso model for the prediction of insulin resistance.