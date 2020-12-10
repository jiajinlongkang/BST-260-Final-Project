library(shiny)
library(shinythemes)
library(DT) # display data table
# prediction/training
library(caret)
library(glmnet)
library(tidyverse)
library(e1071)
library(pROC)
library(PerformanceAnalytics)
library(corrplot)

# All data
merge_data <- readRDS("merge_data.rds")

# Data used for prediction
clean_data_all <- read.csv("clean.csv")
clean_data <- clean_data_all %>% 
    select(-c(X,SEQN,inAnalysis,met_score,DM_insulin,pregtest,FG_mmol,Insulin_pmol,Weight_ogtt,
              HOMA_b,Weight_MEC,Variance_Stratum,Variance_PSU,Weight_fasting))
#Check missingness of variables. Do not use variables with high missingness
missing <- apply(is.na(clean_data)/dim(clean_data)[1], 2, sum)
all_Var <- names(missing)[missing<0.2]

########## For exploratory analysis parts
# for exploratory analysis tabs
explore_data <- clean_data_all %>% 
    select(-c(X,SEQN,inAnalysis,met_score,DM_insulin,pregtest,FG_mmol,Insulin_pmol,Weight_ogtt,
              HOMA_b, Weight_MEC,Variance_Stratum,Variance_PSU,Weight_fasting)) %>% # include homa-b
    select(all_of(all_Var)) %>% 
    mutate(
        Age=as.numeric(Age),
        ApoB_mg=as.numeric(ApoB_mg),
        FVC = as.numeric(FVC),
        FEV1 = as.numeric(FEV1),
        preDM = as.factor(preDM),
        cancer = as.factor(cancer),
        DM_fhx = as.factor(DM_fhx),
        smoking = as.factor(smoking),
        CVD = as.factor(CVD),
        HOMA_bin=as.factor(ifelse(HOMA_IR<=quantile(HOMA_IR,2/3),0,1))
    ) %>% .[complete.cases(.),] 
num_var<-names(select_if(explore_data,is.numeric))
cat_var<-c("preDM","cancer","DM_fhx","CVD","smoking","Gender","Race","Education","Income_ratio","HOMA_bin")



######## For model building
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
#Read in merged data, select potential variables, dichotomize HOMA-IR to normal vs. resistant
# clean_data for binary analysis
clean_data_bin <- clean_data_all %>% drop_na(HOMA_IR) %>%
    select(c(all_Var)) %>% .[complete.cases(.),] %>% 
    mutate(HOMA_bin=as.factor(ifelse(HOMA_IR<=quantile(HOMA_IR,2/3),0,1))) %>%
    select(-HOMA_IR)

clean_data_cont <- clean_data_all %>% 
    select(-c(X,SEQN,inAnalysis,met_score,DM_insulin,pregtest,FG_mmol,Insulin_pmol,Weight_ogtt,
              HOMA_b,Weight_MEC,Variance_Stratum,Variance_PSU,Weight_fasting)) %>% 
    select(all_of(all_Var)) %>% .[complete.cases(.),]

# for correlation matrix
corr_name<-names(clean_data_cont %>% select(-c(Gender,Race, Education, Income_ratio)))

lasso_predict <- function(seed_val,prob){
    set.seed(seed_val)
    train_index <- createDataPartition(clean_data_bin$HOMA_bin, times = 1, p = 0.8, list = FALSE)
    train_set <- clean_data_bin[train_index, ]
    test_set <- clean_data_bin[-train_index, ]
    
    set.seed(seed_val)
    model.lasso <- train(
        HOMA_bin~., data = train_set, method = "glmnet",
        trControl = trainControl("cv", number = 10),
        tuneGrid = expand.grid(alpha=1,lambda=seq(0.001,0.5,by = 0.001)) #Range of lambda
    )
    
    #Display regression coefficients
    result<-coef(model.lasso$finalModel, model.lasso$bestTune$lambda)
    result<-data.frame(coef.name = dimnames(result)[[1]], coef.value = matrix(result)) %>% 
        filter(coef.value != 0 & coef.name !="(Intercept)")
    result$coef.value<-specify_decimal(result$coef.value,5)
    lambda <- model.lasso$bestTune$lambda
    
    #Make predictions on the test data
    pred_value<-predict(model.lasso ,newdata=test_set,type="prob")
    pred_status<-ifelse(pred_value[,2]>prob/100,1,0) 
    conf_mtrix<-confusionMatrix(data=as.factor(pred_status),reference=test_set$HOMA_bin,positive="1")
    #Calculate AUC and draw ROC curve
    roc<-roc(test_set$HOMA_bin,pred_value[,2])
    return(list(
        "model"=model.lasso,
        "variable"=result,
        "lambda"=lambda,
        "predict_val"=pred_value,
        "pred_status"=pred_status,
        "true_status"=test_set$HOMA_bin,
        "confusion"=conf_mtrix,
        "ROC"=roc
    ))
}

knn_predict<-function(seed_val,prob){
    set.seed(seed_val)
    train_index <- createDataPartition(clean_data_bin$HOMA_bin, times = 1, p = 0.8, list = FALSE)
    train_set <- clean_data_bin[train_index, ]
    test_set <- clean_data_bin[-train_index, ]
    
    set.seed(seed_val)
    model.knn <- train(
        HOMA_bin~.,
        data=train_set,method="knn",
        trControl = trainControl("cv",number=10), #10-fold cross validation
        tuneGrid=expand.grid(k=c(2,3,4,seq(5,50,by=2))), #Specify possible k values to evaluate
        metric="Accuracy",
        preProcess = c("center","scale")
    )
    fit_train<-model.knn
    pred_value<-predict(model.knn,newdata=test_set,type="prob")
    pred_status<- ifelse(pred_value[,2]>prob/100,1,0)
    #KNN model performance statistics
    conf_mtrix<-confusionMatrix(data=as.factor(pred_status),reference=test_set$HOMA_bin,positive="1")
    roc<- roc(test_set$HOMA_bin,pred_value[,2])
    return(list(
        "model"=fit_train,
        "k"=fit_train,
        "predict_val"=pred_value,
        "pred_status"=pred_status,
        "true_status"=test_set$HOMA_bin,
        "confusion"=conf_mtrix,
        "ROC"=roc
    ))
}

logit_predict<-function(seed_val,prob){
    set.seed(seed_val)
    train_index <- createDataPartition(clean_data_bin$HOMA_bin, times = 1, p = 0.8, list = FALSE)
    train_set <- clean_data_bin[train_index, ]
    test_set <- clean_data_bin[-train_index, ]
    
    model.logistic <- train(
        HOMA_bin ~., data = train_set, method = "glm",
        trControl = trainControl("cv", number = 10),
        family = "binomial"
    )
    #Make predictions on the test data
    pred_value <- predict(model.logistic,newdata=test_set,type="prob")
    pred_status <- ifelse(pred_value[,2]>0.5,1,0)
    #Logistic model performance statistics
    conf_mtrix<-confusionMatrix(data=as.factor(pred_status),reference=test_set$HOMA_bin,positive="1")
    #Calculate AUC and draw ROC curve
    roc <- roc(test_set$HOMA_bin,pred_value[,2])
    return(list(
        "model"=model.logistic,
        "predict_val"=pred_value,
        "pred_status"=pred_status,
        "true_status"=test_set$HOMA_bin,
        "confusion"=conf_mtrix,
        "ROC"=roc
    ))
}

lasso_predict_cont<-function(seed_val){
    set.seed(seed_val)
    training.samples <- clean_data_cont$HOMA_IR %>% 
        createDataPartition(p = 0.8, list = FALSE)
    train.data  <- clean_data_cont[training.samples, ]
    test.data <- clean_data_cont[-training.samples, ]
    
    # predictor variables
    x <- model.matrix(HOMA_IR~., train.data)[,-1]
    # outcome 
    y <- train.data$HOMA_IR
    
    cv.lasso <- cv.glmnet(x, y, alpha = 1, nfolds = 10)
    
    # Fit the final model on the training data
    
    
    # Make predictions on the test data
    x.test <- model.matrix(HOMA_IR ~., test.data)[,-1]
    y.test <- test.data$HOMA_IR
    model <- glmnet(x, y, alpha = 1, lambda = cv.lasso$lambda.min)  
        
    pred.lasso<-model %>% predict(newx = x.test)
    # Selected variables and coefficients
    coef.las <- data.frame(coef.name = dimnames(coef(model))[[1]], 
                           coef.value = matrix(coef(model)))
    x.select <- coef.las %>% filter(coef.value != 0 & coef.name !="(Intercept)") 
    x.select$coef.value<-specify_decimal(x.select$coef.value,5)
    # RMSE
    RMSE<-RMSE(pred.lasso, y.test)
    
    return(list(
        "model"=cv.lasso,
        "lambda"=cv.lasso$lambda.min,
        "variable"=x.select,
        "RMSE"=RMSE
    ))
}

knn_predict_cont <- function(seed_val){
    set.seed(1)
    training.samples <- clean_data_cont$HOMA_IR %>% 
        createDataPartition(p = 0.8, list = FALSE)
    train.data  <- clean_data_cont[training.samples, ]
    test.data <- clean_data_cont[-training.samples, ]
    model.knn <- train(
        HOMA_IR~.,
        data=train.data, method="knn",
        trControl = trainControl("cv", number=10), #10-fold cross validation
        tuneGrid=expand.grid(k=c(2,3,4,seq(5,50,by=2))), #Specify possible k values to evaluate
        metric="RMSE",
        preProcess = c("center","scale")
    )
    pred.knn <- predict(model.knn, newdata=test.data)
    #Plot model RMSE vs different values of k
    plot(model.knn)
    #Make predictions on the test data
    RMSE(pred.knn,test.data$HOMA_IR)
    return(list(
        "model"=model.knn,
        "k"=model.knn$bestTune,
        "RMSE"=RMSE(pred.knn, test.data$HOMA_IR)
    ))
}

########### define UI
ui <- navbarPage(
    "BST260 Final Project",
     theme = shinytheme("yeti"),
     collapsible = TRUE,
################# Basics #################
    tabPanel("Brief description",
        includeMarkdown("background.md")
             
     ), # close tabPanel Basics             
################# DATA #################
    tabPanel("Data",
             sidebarLayout(
                 # Sidebar panel for selection
                 sidebarPanel(
                    # Data     
                     conditionalPanel(
                         'input.tabset_name === "Data"',
                         selectInput(
                             "data_name","Data set:",
                             c("Original data set" = "merge",
                               "Prediction data set" = "clean")
                             ),
                         radioButtons(
                             "data_type","Summary type:",
                             c("Preview" = "preview",
                               "Structure" = "structure",
                                "Summary" = "summary")
                         ) # radiobuttons
                     ),# close conditionalP1
                     
                     # View
                     conditionalPanel(
                         'input.tabset_name === "View"',
                         selectInput(
                             "data_name2","Data set:",
                             c("Original data set" = "merge",
                               "Prediction data set" = "clean")
                         ),
                         conditionalPanel(
                             'input.data_name2=="merge"',
                             checkboxGroupInput(
                                 "show_vars1","Columns in data sets to show:", 
                                 names(merge_data), selected = names(merge_data))),
                         conditionalPanel(
                             'input.data_name2=="clean"',
                             checkboxGroupInput(
                                 "show_vars2","Columns in data sets to show:", 
                                 names(explore_data), selected = names(explore_data))),
                         helpText("Click the column header to sort a column.")
                     ), # close conditionalP2
                     
                     # Explore
                     conditionalPanel(
                         'input.tabset_name === "Explore"',
                         selectInput(
                             "plot_type", "Plot type:", 
                             choices = list(
                                 "Histogram" = 1,
                                 "Line plot" = 2,
                                 "Box plot" = 3,
                                 "Correlation matrix" = 4
                             ), selected = 4),
                         conditionalPanel(
                             'input.plot_type==1',
                             selectInput(
                                 "show_vars3", "Numeric variables:", 
                                 num_var, selected = NULL),
                             radioButtons(
                                 "x_scale0","Scale transformation (X):",
                                 c("Identity" = "identity",
                                   "Natural log" = "log",
                                   "Log 10" = "log10",
                                   "Square root" = "sqrt")
                             ),
                             helpText("Pick one variable from the list.")
                         ),
                         conditionalPanel(
                             'input.plot_type==2',
                             selectInput(
                                 "show_vars4", "Numeric variables (X):", 
                                 num_var, selected = NULL),
                             radioButtons(
                                 "x_scale","Scale transformation (X):",
                                 c("Identity" = "identity",
                                   "Natural log" = "log",
                                   "Log 10" = "log10",
                                   "Square root" = "sqrt")
                             ), # radiobuttons
                             selectInput(
                                 "show_vars5", "Numeric variables (Y):", 
                                 num_var, selected = NULL),
                             radioButtons(
                                 "y_scale","Scale transformation (Y):",
                                 c("Identity" = "identity",
                                   "Natural log" = "log",
                                   "Log 10" = "log10",
                                   "Square root" = "sqrt")
                             ), # radiobuttons
                             helpText("Pick two different variables from the list.")
                         ),
                         conditionalPanel(
                             'input.plot_type==3',
                             selectInput(
                                 "show_vars6", "Numeric variables:", 
                                 num_var, selected = NULL),
                             radioButtons(
                                 "y_scale0","Scale transformation (Y):",
                                 c("Identity" = "identity",
                                   "Natural log" = "log",
                                   "Log 10" = "log10",
                                   "Square root" = "sqrt")
                             ), # radiobuttons
                             selectInput(
                                 "show_vars7", "Categorical variables:", 
                                 cat_var, selected = NULL),
                             helpText("Pick two different variables from the list.")
                         ),
                         conditionalPanel(
                             'input.plot_type==4',
                             checkboxGroupInput(
                                 "show_vars8","Variables in prediction data set:", 
                                 corr_name, selected = corr_name)
                         )
                     ) # close conditionalP3
                 ), # close sidebarP
                 
                 
                 mainPanel(
                     tabsetPanel(type="tabs", id = 'tabset_name',
                        # tab 1
                        tabPanel("Data",
                          conditionalPanel(
                              'input.data_type == "preview"',
                          br(),
                          p("The original data set includes all the participants in NHANES 2009-2010 (n=10537)."),
                          br(),
                          p(strong("Exclusion criteria for the prediction analysis include:"),"younger than 18 years old, pregant, diagnosed with diabetes before, 
                            taking anti-diabetic treatments, or have not reported HOMA-IR."),
                          br(),
                          p("The final data set for model building includes 1121 participants."),
                          br(),
                          p("The first 10 observations of the selected data set:")
                          ), # conditionalP
                          verbatimTextOutput("data_structure")
                          ), # tabP1
                        
                        # tab 2
                        tabPanel("View",DT::dataTableOutput("data_view")),
                        
                        # tab 3
                        tabPanel("Explore", 
                           p("All the following plots are generated based on",
                             span(" the data set for model building",style = "color:blue"),
                             "."),
                            br(),
                            plotOutput("data_plot"))
                     ) # close tabsetP
                 ) # close mainP
                 
             )#close sidebarL
    ), # close tabPanel Data
    
################# Model #################
    navbarMenu("Model",
               
######## Binary #########
         tabPanel(
             "Binary outcome",
             sidebarLayout(
                 # Sidebar panel for selection
                 sidebarPanel(
                     # LASSO     
                     conditionalPanel(
                         'input.tabset_name2 === "Lasso"',
                         numericInput("seed1", "Choose a number for seed:", value = 1),
                         sliderInput("slider1", "Choose a cut-off probability:", min = 0, 
                                     max = 100, value = 50, post  = "%"),
                         radioButtons("radio1", "Select one output:",
                                      choices = list("Variable selection" = 1, 
                                                     "Optimal lambda" = 2, 
                                                     "ROC curve" = 3,
                                                     "Confusion matrix" = 4), 
                                      selected = 1),
                         helpText("It may take a few minutes to run.")
                     ),# close conditionalP1
                     
                     # KNN
                     conditionalPanel(
                         'input.tabset_name2 === "KNN"',
                         numericInput("seed2", "Choose a number for seed:", value = 1),
                         sliderInput("slider2", "Choose a cut-off probability:", min = 0, 
                                     max = 100, value = 50),
                         radioButtons("radio2", "Select one output:",
                                      choices = list("Optimal k" = 1, 
                                                     "ROC curve" = 2,
                                                     "Confusion matrix" = 3), 
                                      selected = 1),
                         helpText("It may take a few minutes to run.")
                     ), # close conditionalP2
                     
                     # Compare
                     conditionalPanel(
                         'input.tabset_name2 === "Comparison"',
                         numericInput("seed3", "Choose a number for seed:", value = 1),
                         sliderInput("slider3", "Choose a cut-off probability:", min = 0, 
                                     max = 100, value = 50),
                         radioButtons("radio3", "Select one output:",
                                      choices = list("Confusion matrix on test set" = 1, 
                                                     "Accuracy and Kappa comparison based on 10-fold cross-validation on train set" = 2), 
                                      selected = 1),
                         helpText("It may take a few minutes to run.")
                     ) # close conditionalP3
                 ), # close sidebarPanel
                 
                 
                 mainPanel(
                     tabsetPanel(type="tabs", id = 'tabset_name2',
                         # tab 1
                         tabPanel(
                             "Lasso",
                             conditionalPanel(
                                 'input.radio1==1',
                                 br(),
                                 strong("Variables selected from LASSO:"),
                                 br(),
                                 tableOutput("lasso_bin_var")
                             ),
                             conditionalPanel(
                                 'input.radio1==2',
                                 br(),
                                 strong("Optimal lambda based on cross-validation in training set:"),
                                 br(),
                                 verbatimTextOutput("lasso_bin_lambda")
                             ),
                             conditionalPanel(
                                 'input.radio1==3',
                                 br(),
                                 strong("ROC curve based on selected probability cutoff:"),
                                 br(),
                                 plotOutput("lasso_bin_roc"),
                                 p("AUC is:"),
                                 verbatimTextOutput("lasso_bin_AUC")
                             ),
                             conditionalPanel(
                                 'input.radio1==4',
                                 br(),
                                 strong("Confusion matrix based on selected probability cutoff:"),
                                 verbatimTextOutput("lasso_bin_confusion"),
                                 
                             )
                         ), # tabP1
                         
                         # tab 2
                         tabPanel(
                             "KNN",
                             conditionalPanel(
                                 'input.radio2==1',
                                 br(),
                                 strong("Optimal k based on cross-validation in training set:"),
                                 br(),
                                 verbatimTextOutput("knn_bin_k"),
                                 br(),
                                 plotOutput("knn_bin_cv")
                             ),
                             conditionalPanel(
                                 'input.radio2==2',
                                 br(),
                                 strong("ROC curve based on selected probability cutoff:"),
                                 br(),
                                 plotOutput("knn_bin_roc"),
                                 br(),
                                 p("AUC is:"),
                                 verbatimTextOutput("knn_bin_AUC")
                             ),
                             conditionalPanel(
                                 'input.radio2==3',
                                 br(),
                                 strong("Confusion matrix based on selected probability cutoff:"),
                                 br(),
                                 verbatimTextOutput("knn_bin_confusion")
                             )
                         ),
                         
                         # tab 3
                         tabPanel(
                             "Comparison",
                             conditionalPanel(
                                 'input.radio3==1',
                                 br(),
                                 strong("Confusion matrix based on selected probability cutoff:"),
                                 br(),
                                 strong("Logistic regression"),
                                 verbatimTextOutput("compare_bin_logit"),
                                 br(),
                                 strong("Lasso regression"),
                                 verbatimTextOutput("compare_bin_lasso"),
                                 br(),
                                 strong("KNN"),
                                 verbatimTextOutput("compare_bin_knn")
                             ),
                             conditionalPanel(
                                 'input.radio3==2',
                                 br(),
                                 strong("Comparison of three methods:"),
                                 br(),
                                 plotOutput("accuracy_bin")
                             )
                         )
                     ) # close tabsetP
                 ) # close mainP
                 
             ) #close sidebarLayout
        ), # close tabPanel for binary
        
        
######## Continuous #########        
        tabPanel(
            "Continuous outcome",
            sidebarLayout(
                # Sidebar panel for selection
                sidebarPanel(
                    # LASSO     
                    conditionalPanel(
                        'input.tabset_name3 === "Lasso"',
                        numericInput("seed4", "Choose a number for seed:", value = 1),
                        radioButtons("radio4", "Select one output:",
                                     choices = list("Variable selection" = 1, 
                                                    "Optimal lambda/RMSE" = 2), 
                                     selected = 1),
                        helpText("It may take a few minutes to run.")
                    ),# close conditionalP1
                    
                    # KNN
                    conditionalPanel(
                        'input.tabset_name3 === "KNN"',
                        numericInput("seed5", "Choose a number for seed:", value = 1),
                        radioButtons("radio5", "Select one output:",
                                     choices = list("Optimal k" = 1, 
                                                    "RMSE" = 2), 
                                     selected = 1),
                        helpText("It may take a few minutes to run.")
                    ) # close conditionalP2
                ), # close sidebarPanel
                
                
                mainPanel(
                    tabsetPanel(type="tabs", id = 'tabset_name3',
                        # tab 1
                        tabPanel(
                            "Lasso",
                            conditionalPanel(
                                'input.radio4==1',
                                br(),
                                strong("Variables selected from LASSO:"),
                                br(),
                                tableOutput("lasso_cont_var")
                            ),
                            conditionalPanel(
                                'input.radio4==2',
                                br(),
                                strong("Optimal lambda based on cross-validation in training set:"),
                                br(),
                                verbatimTextOutput("lasso_cont_lambda"),
                                br(),
                                plotOutput("lasso_cont_CV"),
                                br(),
                                strong("RMSE in test set:"),
                                verbatimTextOutput("lasso_cont_RMSE")
                            )
                        ), # tabP1
                        
                        # tab 2
                        tabPanel(
                            "KNN",
                            conditionalPanel(
                                'input.radio5==1',
                                br(),
                                strong("Optimal k selected from KNN:"),
                                br(),
                                verbatimTextOutput("knn_cont_k")
                            ),
                            conditionalPanel(
                                'input.radio5==2',
                                br(),
                                strong("RMSE:"),
                                br(),
                                verbatimTextOutput("knn_cont_RMSE"),
                                br(),
                                plotOutput("knn_cont_CV")
                            )
                        )
                    ) # close tabsetP
                ) # close mainP
                
            ) #close sidebarLayout
        ) # close tabPanel for binary     
             
             
    ) # close navbar Model
)# close UI






server <- function(input, output) {
################# DATA #################
    # Tab 1
    output$data_structure <- renderPrint({
        if (input$data_name=="merge"){
            if (input$data_type=="preview"){
                head(merge_data)
            } else if (input$data_type=="structure"){
                str(merge_data)
            } else {
                summary(merge_data)
            }
        } else if (input$data_name=="clean"){
            if (input$data_type=="preview"){
                head(explore_data)
            } else if (input$data_type=="structure"){
                str(explore_data)
            } else {
                summary(explore_data)
            }
        }
    })
    
    # Tab 2
    output$data_view <- DT::renderDataTable({
        if (input$data_name2=="merge"){
            DT::datatable(merge_data[, input$show_vars1, drop=FALSE])
        } else if (input$data_name2=="clean"){
            DT::datatable(explore_data[, input$show_vars2, drop=FALSE])
        }
    }) # renderDataTable
    
    
    # Tab 3
    output$data_plot<-renderPlot({
        if (input$plot_type==1){
            explore_data %>%  
                ggplot(aes_string(x=input$show_vars3)) + 
                geom_histogram(aes(y=..density..),colour="black", fill="white") +
                geom_density(alpha=.2, fill="light blue") + 
                scale_x_continuous(trans = input$x_scale0) +
                xlab(paste0(input$show_vars3," in ",input$x_scale0," scale")) +
                ggtitle(paste0("Distribution of ", input$show_vars3)) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5))
        } else if (input$plot_type==2){
            explore_data %>%  
                ggplot(aes_string(x=input$show_vars4,y=input$show_vars5)) +
                geom_point(color="grey") +
                geom_smooth(method = "lm") + 
                scale_x_continuous(trans = input$x_scale) +
                scale_y_continuous(trans = input$y_scale) +
                xlab(paste0(input$show_vars4," in ",input$x_scale," scale")) +
                ylab(paste0(input$show_vars5," in ",input$y_scale," scale")) +
                ggtitle(paste0(input$show_vars5," vs ",input$show_vars4)) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5))
        } else if (input$plot_type==3){
            explore_data %>%  
                ggplot(aes_string(y=input$show_vars6,x=input$show_vars7, fill=input$show_vars7)) +
                geom_boxplot() +
                scale_y_continuous(trans = input$y_scale0) +
                ylab(paste0(input$show_vars6," in ",input$y_scale0," scale")) +
                ggtitle(paste0(input$show_vars6," by ",input$show_vars7)) +
                scale_fill_brewer(palette="Set3") +
                theme_bw() +
                theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank()) +
                theme(plot.title = element_text(hjust = 0.5))
        } else if (input$plot_type==4){
            clean_data_cont %>% select(input$show_vars8) %>% 
                cor(., method = "spearman",use = "complete.obs") %>%  #, use = "complete.obs"
                corrplot(., order = "hclust", tl.cex=0.5, tl.col="black", 
                         addrect = 2)
            
        }
    })
    
################# Model #################
    ### Binary
    # Tab 1
    bin_lasso_1<-reactive(lasso_predict(seed_val = input$seed1,prob = input$slider1)) 
    
    output$lasso_bin_var<-renderTable(
        bin_lasso_1()$variable
    )
    
    output$lasso_bin_lambda<-renderText({bin_lasso_1()$lambda})
    
    output$lasso_bin_roc<-renderPlot({
        bin_lasso_1() %>% .$ROC %>% 
        ggroc(legacy.axes=T,col="red",lwd=2) +
            geom_abline(slope=1,intercept=0) + # add identity line
            theme(
                panel.background = element_blank(), 
                axis.title.x = element_text(size =18, face = 'bold'),
                axis.title.y = element_text(size =18, face = 'bold'),
                panel.border = element_rect(size = 2, fill = NA), 
                axis.text.x = element_text(size = 14, face ='bold'),
                axis.text.y = element_text(size = 14, face ='bold')) +
            xlab('Specificity %') +
            ylab('Sensitivity %') +
            scale_x_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25) * 100) + 
            scale_y_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25) * 100) +
         #   annotate("text",x=.25,y=.75,size=6,label=paste0("AUC=",round(auc(roc.lasso),4))) +
            ggtitle("ROC curve\n(lasso regression, binary outcome)") +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5))
    })
    output$lasso_bin_AUC<-renderText({auc(bin_lasso_1()$ROC)})
    output$lasso_bin_confusion<-renderPrint({
        #bin_lasso_1()$confusion
        bin_lasso_1()$confusion
    }
    )

    # Tab 2
    bin_knn_1<-reactive(knn_predict(seed_val = input$seed2,prob = input$slider2)) 
    
    output$knn_bin_k<-renderPrint({
       paste("k =",bin_knn_1()$k$k) 
    })
    output$knn_bin_cv<-renderPlot({
        plot(bin_knn_1()$model,xlab="Number of Neighbors",
             main="Cross-Validation in the training set")
    })
    output$knn_bin_roc<-renderPlot({
        bin_knn_1() %>% .$ROC %>% 
            ggroc(legacy.axes=T,col="red",lwd=2) +
            geom_abline(slope=1,intercept=0) + # add identity line
            theme(
                panel.background = element_blank(), 
                axis.title.x = element_text(size =18, face = 'bold'),
                axis.title.y = element_text(size =18, face = 'bold'),
                panel.border = element_rect(size = 2, fill = NA), 
                axis.text.x = element_text(size = 14, face ='bold'),
                axis.text.y = element_text(size = 14, face ='bold')) +
            xlab('Specificity %') +
            ylab('Sensitivity %') +
            scale_x_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25) * 100) + 
            scale_y_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25) * 100) +
            #   annotate("text",x=.25,y=.75,size=6,label=paste0("AUC=",round(auc(roc.lasso),4))) +
            ggtitle("ROC curve\n(KNN prediction, binary outcome)") +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5))
    })
    output$knn_bin_AUC<-renderText({auc(bin_knn_1()$ROC)})
    output$knn_bin_confusion<-renderPrint({
        bin_knn_1()$confusion
    })
    
    bin_lasso_2<-reactive(lasso_predict(seed_val = input$seed3,prob = input$slider3))
    bin_knn_2<-reactive(knn_predict(seed_val = input$seed3,prob = input$slider3))
    bin_logit<-reactive(logit_predict(seed_val = input$seed3,prob = input$slider3))
    
    # Tab 3
    output$compare_bin_logit<-renderPrint({bin_logit()$confusion})
    output$compare_bin_lasso<-renderPrint({bin_lasso_2()$confusion})
    output$compare_bin_knn<-renderPrint({bin_knn_2()$confusion})
    output$accuracy_bin<-renderPlot({
        resamples(list(
            "Logistic"=bin_logit()$model,
            "Lasso"=bin_lasso_2()$model,
            "knn"=bin_knn_2()$model)) %>% 
        dotplot(.,metric=c("Accuracy","Kappa"))
    })
    
    ### Continuous
    # Tab 1
    cont_lasso<-reactive(lasso_predict_cont(seed_val = input$seed4)) 
    
    output$lasso_cont_var<-renderTable(
        cont_lasso() %>% .$variable 
    )
    
    output$lasso_cont_lambda<-renderPrint({cont_lasso() %>% .$lambda})
    output$lasso_cont_CV<-renderPlot({
        plot(cont_lasso() %>% .$model)
    })
    output$lasso_cont_RMSE<-renderPrint({cont_lasso() %>% .$RMSE})
                        

    # Tab 2
    cont_knn<-reactive(knn_predict_cont(seed_val = input$seed5)) 
    
    output$knn_cont_k<-renderText({cont_knn()$k$k})
    output$knn_cont_RMSE<-renderText({cont_knn() %>% .$RMSE})
    output$knn_cont_CV<-renderPlot({
        plot(cont_knn() %>% .$model,main = "Cross-validation in the training set")
    })
    
    
} # close server

# Run the application 
shinyApp(ui = ui, server = server)
