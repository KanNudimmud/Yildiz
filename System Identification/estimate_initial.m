%% Estimation of the System with Initial Model Structure Selection
%% Pre-processing the Data
load dry2 % hair-dryer data

% Form a data set for estimation of the first half, and a reference set for validation purposes of the second half
ze = dry2(1:500);
zr = dry2(501:1000);

% Detrend each of the sets
ze = detrend(ze);
zr = detrend(zr);

% Let us look at a portion of the estimation data
figure,plot(ze(200:350))

%% Estimating Input Delay
% Using delayest
delay = delayest(ze) % na = nb = 2 is used, by default

%  To gain insight into how delayest works, let us evaluate the loss function for various choices of delays explicitly. We select a second order model (na=nb=2), which is the default for delayest, and try out every time delay between 1 and 10. The loss function for the different models are computed using the validation data set
V = arxstruc(ze,zr,struc(2,2,1:10));

% We now select that delay that gives the best fit for the validation data
[nn,Vm] = selstruc(V,0); % nn is given as [na nb nk]
% which show the best model has a delay of nn(3) = 3.
% The choice of 3 delays is thus rather clear, since the corresponding loss is minimum.

% Using impulse
FIRModel = impulseest(ze);

%  plot this response with a confidence interval represented by 3 standard deviations
clf
h = impulseplot(FIRModel);
showConfidence(h,3)
% The filled light-blue region shows the confidence interval for the insignificant response in this estimation. There is a clear indication that the impulse response "takes off" (leaves the uncertainty region) after 3 samples. This points to a delay of three intervals.

% Using n4sid based state-space evaluation
% We may also estimate a family of parametric models to find the delay corresponding to the "best" model. In case of state-space models, a range of orders may be evaluated simultaneously and the best order picked from a Hankel Singular Value plot. Execute the following command to invoke n4sid in an interactive mode
m = n4sid(ze,1:15); % All orders between 1 and 15.

% The plot indicates an order of 3 as the best value. For this choice, let us compute the impulse response of the model m
m = n4sid(ze, 3);
showConfidence(impulseplot(m),3)
% As with non-parametric impulse response, there is a clear indication that the delay from input to output is of three samples.

%% Determining Model Order
%  Functions arxstruc and selstruc may be used for choosing the best order for ARX models. For our example, let us check the fit for all 100 combinations of up to 10 b-parameters and up to 10 a-parameters, all with a delay value of 3
V = arxstruc(ze,zr,struc(1:10,1:10,3));

% The best fit for the validation data set is obtained fo
nn = selstruc(V,0)

nns = selstruc(V) %invoke selstruc in an interactive mode

%The best fit is thus obtained for nn = [4 4 3], while we see that the improved fit compared to nn = [2 2 3] is rather marginal.
%We may also approach this problem from the direction of reducing a higher order model. If the order is higher than necessary, then the extra parameters are basically used to "model" the measurement noise. These "extra" poles are estimated with a lower level of accuracy (large confidence interval). If their are cancelled by a zero located nearby, then it is an indication that this pole-zero pair may not be required to capture the essential dynamics of the system.
%For our example, let us compute a 2th order model:
%The confidence intervals for the two complex-conjugate poles and zeros overlap, indicating they are likely to cancel each other. Hence, a second order model might be adequate. Based on this evidence, let us compute a 2nd order ARX model:
th2 = arx(ze,[2 2 3]);

%The plot indicates that there was no significant loss of accuracy in reducing the order from 4 to 2. We can also check the residuals ("leftovers") of this model, i.e., what is left unexplained by the model.
e = resid(ze,th2);
plot(e(:,1,[])), title('The residuals')

%We see that the residuals are quite small compared to the signal level of the output, that they are reasonably well (although not perfectly) uncorrelated with the input and among themselves. We can thus be (provisionally) satisfied with the model th2.
%Let us now check if we can determine the model order for a state-space structure. As before, we know the delay is 3 samples. We can try all orders from 1 to 15 with a total lag of 3 samples in n4sid. Execute the following command to try various orders and choose one interactively.
ms = n4sid(ze,[1:15],'InputDelay',2); %n4sid estimation with variable orders

%The "InputDelay" was set to 2 because by default n4sid estimates a model with no feedthrough (which accounts for one sample lag between input and output). The default order, indicated in the figure above, is 3, that is in good agreement with our earlier findings. Finally, we compare how the state-space model ms and the ARX model th2 compare in reproducing the measured output of the validation data:
ms = n4sid(ze,3,'InputDelay',2);
compare(zr,ms,th2)
%The comparison plot indicates that the two models are practically identical.

%% end.