%% Probability Density Function (PDF)
%% Create a Dataset
mu = 0;
sigma = 1;
data = randn(6000,1)*sigma+mu;

% Display the data
figure(1)
histogram(data,80)

%% Find mean,variance and standart deviation
% Compute mean 
n     = size(data,1);
meanDat = sum(data) / n;

% Compute variance
varDat = 0;

for i=1:length(data)
    varDat = varDat + (data(i)-meanDat).^2;
end

varDat = varDat / (n-1);

% Compute standard deviation
dataM = data-meanDat;

stdDat = 0;
for i=1:length(data)
    stdDat = stdDat + dataM(i).^2;
end

stdDat = sqrt(stdDat / (n-1));

%% Compute the covariance
%% Method 1 - Loop
% Initialize covariance matrix
covmatL = zeros(n);

% Nested loop and compute dot product scaled by N-1
for i=1:n
    for j=1:n
        
        % Mean-centered data
        centi = data(i,:) - mean(data(i,:));
        centj = data(j,:) - mean(data(j,:));
        
        % Compute covariance
        covmatL(i,j) = sum(centi.*centj) / (n -1);
    end 
end 

%% Method 2 - Matrix multiplication
% Mean-center over time
dataM = bsxfun(@minus,data,mean(data,2));

% Pairwise dot product 
covmatM = dataM*dataM' / (n-1);

%% Display PDF
% Calculate PDF
time  = -4:.1:4;
fx    = (1/sqrt(2*pi*(sigma^2))) * exp(((time-mu).^2) ./(-2*(sigma^2)))*600;
% Note : Multiplaying with 600 is for normalization respect to histogram
% Plot top op histogram
figure(1), hold on
plot(time,fx,'r','linew',3)

%% end