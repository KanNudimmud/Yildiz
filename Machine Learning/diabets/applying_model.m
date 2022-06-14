%% Machine Learning Project (Applying Model)

%% Load the data
data = readtable('diabetes.csv');

rng(1); % For reproducibility

%% Exploratory Data Analysis
headers = {'Pregnancies','Glucose','BloodPressure',...
    'SkinThickness','Insulin','BMI',...
    'DiabetesPedigreeFunction','Age'};
% Who have diabetes
diabet = data(find(data.Outcome == 1),:);

% Healthy 
healthy = data(find(data.Outcome == 0),:);

% Convert tables to double for computations
diabetD  = table2array(diabet);
healthyD = table2array(healthy);

% Bring together into data and convert
dataD = [diabetD; healthyD];

% Visualize summary statistics with box plot
figure, boxplot(diabetD(:,1:end-1))
xticklabels(headers)
title('Diabets Data Statistics with patients')

figure, boxplot(healthyD(:,1:end-1))
xticklabels(headers)
title('Diabets Data Statistics with healthy people')

% Compute correlation matrix
cormat = corr(dataD(:,1:end-1));
cormatL = tril(cormat);

% Make the heatmap
figure
h =heatmap(headers,headers,cormatL)
h.Title = 'Correlation Matrix';
colormap jet

%% Normalization

% Normalize the  diabetes data
% Determine max. values
max_diab = max(diabetD(:,1:end-1));

% Divide each column by their max. value
dataSD = bsxfun(@rdivide,diabetD(:,1:end-1),max_diab);

% Normalize the  healthy data
% Determine max. values
max_heal = max(healthyD(:,1:end-1));

% Divide each column by their max. value
dataSH = bsxfun(@rdivide,healthyD(:,1:end-1),max_heal);

%% Creating Model
% Seperate to train and test samples (%70 for training)
dataM = [dataSD;dataSH];

dataD = [diabetD;healthyD];

X = dataM(:,1:end);
Y = dataD(:,end);
cvpart = cvpartition(Y,'holdout',0.20);
Xtrain = X(training(cvpart),:);
Ytrain = Y(training(cvpart),:);
Xtest  = X(test(cvpart),:);
Ytest  = Y(test(cvpart),:);

% Fit a model
mdl4 = fitcensemble(Xtrain,Ytrain,'OptimizeHyperparameters','auto');
mdl = fitcsvm(Xtrain,Ytrain,'OptimizeHyperparameters','auto');
mdl2 = fitctree(Xtrain,Ytrain,'OptimizeHyperparameters','auto');
mdl3 = fitcnb(Xtrain,Ytrain,'OptimizeHyperparameters','auto');
mdl5 = fitcknn(Xtrain,Ytrain,'OptimizeHyperparameters','auto');

% Plot misclassification of the test data
figure
plot(loss(mdl4,Xtest,Ytest,'mode','cumulative'))
xlabel('Number of trees')
ylabel('Test classification error')

%% Evaluate the Model
% F-score
% Predict labels for the measurements in the second half of the data by using the trained classifier.
Yhat = predict(mdl,Xtest);
Yhat2 = predict(mdl2,Xtest);
Yhat3 = predict(mdl3,Xtest);
Yhat4 = predict(mdl4,Xtest);
Yhat5 = predict(mdl5,Xtest);

% Specify the group order and display the confusion matrix for the resulting classification.
C  = confusionmat(Ytest,Yhat);
C2 = confusionmat(Ytest,Yhat2);
C3 = confusionmat(Ytest,Yhat3);
C4 = confusionmat(Ytest,Yhat4);
C5 = confusionmat(Ytest,Yhat5);

% Plot the confusion matrix C
figure, 
subplot(231), confusionchart(C),title('SVM')
subplot(232), confusionchart(C2),title('Decision Tree')
subplot(233), confusionchart(C3),title('Naive Bayes')
subplot(234), confusionchart(C4),title('Ensemble')
subplot(235), confusionchart(C5),title('KNN')

% Use the function 
stats = statsOfMeasure(C, 1);
stats2 = statsOfMeasure(C2, 1);
stats3 = statsOfMeasure(C3, 1);
stats4 = statsOfMeasure(C4, 1);
stats5 = statsOfMeasure(C5, 1);

%% end.