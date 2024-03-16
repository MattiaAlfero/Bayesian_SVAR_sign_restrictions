# Bayesian SVAR with sign restrictions

The purpose of the code is to estimate a Bayesian SVAR with 2 different priors: flat and minnesota. Then the variance-covariance matrix is identified using sign restrictions. The code is very general, by changing the parameters in the section "Import the data and set main parameters", it can be readapted to datasets with any number of variables. Moreover it can be used to identify any type of shocks with sign restriction methodology by changing the variable "set.restrictions". Since we use the flat and minnesota priors, we don't have to set any prior parameter.

The main script is "Bayesian_SVAR_sign_restrictions.m", it is divided as follow:
1. Import the data and set the parameters, if you want to use the code, you need to readapt the import of the data and the manual parameters.
2. Estimation with flat and minnesota prior and identification with sign restriction: then it plots the IRFs, notice that the plot is made automatically.
3. Counterfactual path of the series made on the latest part of the sample using only the selected shock.
   
I created 2 main functions: 1. "minnesota_prior": very general function to create the minnesota prior; 2. "sign_restrictions": it draws a matrix exploiting the QR decomposition and it checks if it satisfies or not the sign restriction. 

In this specific case the code is used to identify 3 shocks: demand, energy and supply. The dataset is quarterly, based on U.S. data (downloaded from FRED), which starts in 1989:I and ends in 2023:IV. It includes n = 3 variables: (i) the logarithm of real Personal Consumption Expenditure (C); (ii) the logarithm of the chain-type price index of Personal Consumption Expenditure excluding Food and Energy (coreP); and (iii) the logarithm of the chain-type price index of Personal Consumption Expenditure in Energy Goods and Services (energyP).

Notice that the pdf file "Bayesian_SVAR_sign_restrictions.pdf" explains step by step the code and the formula used for the estimation. 

Let me know if you find this useful. Byee!!
