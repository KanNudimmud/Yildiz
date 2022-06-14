%% Machine Learning Project (EDA Part)

%% Load the data
data = readtable('diabetes.csv');

rng(1); % For reproducibility

%% Exploratory Data Analysis
headers = {'Pregnancies','Glucose','BloodPressure',...
    'SkinThickness','Insulin','BMI',...
    'DiabetesPedigreeFunction','Age'};

% Data dimension
size(data)

% Dataframe contents
head(data)

% Check for missing values
sum(isempty(data))

% Who have diabetes
diabet = data(find(data.Outcome == 1),:);
head(diabet)

% Healthy 
healthy = data(find(data.Outcome == 0),:);
head(healthy)

% Find how many missing value for every feature
sum(ismissing(diabet))
sum(ismissing(healthy))

% Replace zeros with avarage values (do it seperately)
diabet.Glucose(diabet.Glucose == 0) = mean(diabet.Glucose);
diabet.BloodPressure(diabet.BloodPressure == 0) = mean(diabet.BloodPressure);
diabet.SkinThickness(diabet.SkinThickness == 0) = mean(diabet.SkinThickness);
diabet.Insulin(diabet.Insulin == 0) = mean(diabet.Insulin); % should ask the teacher
diabet.BMI(diabet.BMI == 0) = mean(diabet.BMI);

healthy.Glucose(healthy.Glucose == 0) = mean(healthy.Glucose);
healthy.BloodPressure(healthy.BloodPressure == 0) = mean(healthy.BloodPressure);
healthy.SkinThickness(healthy.SkinThickness == 0) = mean(healthy.SkinThickness);
healthy.Insulin(healthy.Insulin == 0) = mean(healthy.Insulin); % should ask the teacher
healthy.BMI(healthy.BMI == 0) = mean(healthy.BMI);

% Convert tables to double for computations
diabetD  = table2array(diabet);
healthyD = table2array(healthy);

% Bring together into data and convert
dataD = [diabetD; healthyD];

% Piechart to show either diabetes or non-diabetes
figure
pie([sum(dataD(:,end) == 1),sum(dataD(:,end) == 0)])
legend({'Diabetes','Healthy'})
title("Ratio of Patient's Conditions")

% histograms for every feature
% figure
% for i=1:size(dataD,2)
%     subplot(3,3,i)
%     histogram(dataD(:,i))
%     title(headers{i})
% end

% Descriptive statistics for who have diabetes
% compute mean
mdataDD = mean(diabetD);

% compute the standard deviation
sigmaDD = std(diabetD);

% compute the median
mediandatDD = median(diabetD);

% Loop over every feauture for who have diabetes
% Initialize matrices
Q1 = zeros(1,size(diabetD,2));
Q2 = zeros(1,size(diabetD,2));
Q3 = zeros(1,size(diabetD,2));
Noutliers = zeros(1,size(diabetD,2));

for i=1:size(diabetD,2)
    % rank-transform the data and scale to 1
    y = tiedrank(diabetD(:,i)) / size(diabetD,1);

    % find the values closest to 25% and 75% of the data
    q1 = dsearchn(y,.25);
    q3 = dsearchn(y,.75);

    % get two values in the data
    Q1(:,i) = diabetD(q1,i);
    Q3(:,i) = diabetD(q3,i);

    % compute Interquartile Range (IQR)
    IQR = Q3(:,i)-Q1(:,i);

    % Compute Semi Interquartile Deviation (SID)
    SID = IQR/2;

    % determine extreme Q1 outliers (e.g., x < Q1 - 3*IQR)
    ix = find(diabetD(:,i)<Q1(:,i)-3*IQR);
    if length(ix)>0
        outliersQ1 = diabetD(ix,i);
    else
        outliersQ1 = [];
    end

    % determine extreme Q3 outliers (e.g., x > Q1 + 3*IQR)
    iy = find(diabetD(:,i)>Q1(:,i)+3*IQR);
    if length(iy)>0
        outliersQ3 = diabetD(iy,i);
    else
        outliersQ3 = [];
    end

    % compute total number of outliers
    Noutliers(:,i) = length(outliersQ1)+length(outliersQ3);

    % replace outlier with mean data
    diabetD(ix,i) = mdataDD(i);
    diabetD(iy,i) = mdataDD(i);
end

% display results
disp(['Mean:                                ',num2str(mdataDD)]);
disp(['Standard Deviation:                  ',num2str(sigmaDD)]);
disp(['Median:                              ',num2str(mediandatDD)]);
disp(['Minimum Value:                       ',num2str(min(diabetD))]);
disp(['Maximum Value:                       ',num2str(max(diabetD))]);
disp(['25th Percentile:                     ',num2str(Q1)]);
disp(['50th Percentile:                     ',num2str(mediandatDD)]);
disp(['75th Percentile:                     ',num2str(Q3)]);
disp(['Semi Interquartile Deviation:        ',num2str(SID)]);
disp(['Number of outliers:                  ',num2str(Noutliers)]);

% Descriptive statistics for who have not diabetes
% compute mean
mdataDH = mean(healthyD);

% compute the standard deviation
sigmaDH = std(healthyD);

% compute the median
mediandatDH = median(healthyD);

% Loop over every feauture for who have diabetes
% Initialize matrices
Q1 = zeros(1,size(healthyD,2));
Q2 = zeros(1,size(healthyD,2));
Q3 = zeros(1,size(healthyD,2));
Noutliers = zeros(1,size(healthyD,2));

for i=1:size(healthyD,2)
    % rank-transform the data and scale to 1
    y = tiedrank(healthyD(:,i)) / size(healthyD,1);

    % find the values closest to 25% and 75% of the data
    q1 = dsearchn(y,.25);
    q3 = dsearchn(y,.75);

    % get two values in the data
    Q1(:,i) = healthyD(q1,i);
    Q3(:,i) = healthyD(q3,i);

    % compute Interquartile Range (IQR)
    IQR = Q3(:,i)-Q1(:,i);

    % Compute Semi Interquartile Deviation (SID)
    SID = IQR/2;

    % determine extreme Q1 outliers (e.g., x < Q1 - 3*IQR)
    ix = find(healthyD(:,i)<Q1(:,i)-3*IQR);
    if length(ix)>0
        outliersQ1 = healthyD(ix,i);
    else
        outliersQ1 = [];
    end

    % determine extreme Q3 outliers (e.g., x > Q1 + 3*IQR)
    iy = find(healthyD(:,i)>Q1(:,i)+3*IQR);
    if length(iy)>0
        outliersQ3 = healthyD(iy,i);
    else
        outliersQ3 = [];
    end

    % compute total number of outliers
    Noutliers(:,i) = length(outliersQ1)+length(outliersQ3);

    % replace outlier with mean data
    healthyD(ix,i) = mdataDH(i);
    healthyD(iy,i) = mdataDH(i);
end

% display results
disp(['Mean:                                ',num2str(mdataDH)]);
disp(['Standard Deviation:                  ',num2str(sigmaDH)]);
disp(['Median:                              ',num2str(mediandatDH)]);
disp(['Minimum Value:                       ',num2str(min(healthyD))]);
disp(['Maximum Value:                       ',num2str(max(healthyD))]);
disp(['25th Percentile:                     ',num2str(Q1)]);
disp(['50th Percentile:                     ',num2str(mediandatDH)]);
disp(['75th Percentile:                     ',num2str(Q3)]);
disp(['Semi Interquartile Deviation:        ',num2str(SID)]);
disp(['Number of outliers:                  ',num2str(Noutliers)]);

% Visualize summary statistics with box plot
figure, boxplot(diabetD(:,1:end-1))
xticklabels(headers)
title('Diabets Data Statistics with patients')

figure, boxplot(healthyD(:,1:end-1))
xticklabels(headers)
title('Diabets Data Statistics with healthy people')

dataD = [diabetD;healthyD];

% Compute correlation matrix
cormat = corr(dataD(:,1:end-1));
cormatL = tril(cormat);

% Make the heatmap
figure
h =heatmap(headers,headers,cormatL)
h.Title = 'Correlation Matrix';
colormap jet

% Scatter Plots
figure
scatter(diabetD(:,8),diabetD(:,1))
hold on,
scatter(healthyD(:,8),healthyD(:,1))
legend({'Diabetes','Healthy'})
xlabel('Age'),ylabel('Pregnancies')
title('Age vs. Pregnancies')

figure
scatter(diabetD(:,2),diabetD(:,5))
hold on,
scatter(healthyD(:,2),healthyD(:,5))
legend({'Diabetes','Healthy'})
xlabel('Glucose'),ylabel('Insulin')
title('Glucose vs. Insulin')

figure
scatter(diabetD(:,6),diabetD(:,4))
hold on,
scatter(healthyD(:,6),healthyD(:,4))
legend({'Diabetes','Healthy'})
xlabel('BMI'),ylabel('SkinThickness')
title('BMI vs. SkinThickness')

%% Data Normalization and outliers
% z-score for who have disease
datazD = (diabetD-mdataDD) / sigmaDD; 

mdatazD = mean(datazD);
sdatazD = std(datazD);

figure
plot(datazD,'s','markersize',8,'markerfacecolor','r')
xlabel('Data index'), ylabel('Data value')
title([ 'Mean = ' num2str(round(mdatazD,2)) '; std = ' num2str(round(sdatazD,2)) ])

% min-max scaling for who have disease
% get min and max
dataMinDD = min(diabetD);
dataMaxDD = max(diabetD);

% now min-max scale
dataSD = zeros(268,8);
for i=1:8
    dataSD(:,i) = (diabetD(:,i)-dataMinDD(i)) / (dataMaxDD(i)-dataMinDD(i));
end

% min-max scaling for who have not disease
% get min and max
dataMinDH = min(healthyD);
dataMaxDH = max(healthyD);

% now min-max scale
dataSH = zeros(500,8);
for i=1:8
    dataSH(:,i) =  (healthyD(:,i)-dataMinDH(i)) / (dataMaxDH(i)-dataMinDH(i));
end

%% Normalization (different normalization methods)

% Normalize the  diabetes data
% Determine max. values
max_diab = max(diabetD);

% Divide each column by their max. value
dataSD = bsxfun(@rdivide,diabetD,max_diab);

% Normalize the  healthy data
% Determine max. values
max_heal = max(healthyD);

% Divide each column by their max. value
dataSH = bsxfun(@rdivide,healthyD,max_heal);

%% Feature Extraction
% Concatanete data
dataM = [dataSD;dataSH];

% Glucose-Insulin
F1 = dataM(:,2) .* dataM(:,5);
F1 = F1/max(F1);

% BMI-SkinThickness
F2 = dataM(:,6) .* dataM(:,4);
F2 = F2 / max(F2);

% Add features to the data
dataF = [dataM(:,1:end-1),F1,F2];

% Show data in a table format
dataT = array2table(dataF,'VariableNames',[headers {'Glucose/Insulin',...
    'BMI/SkinThickness'}]);
head(dataT)

%% end.