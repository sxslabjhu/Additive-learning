clc;clear;

node=[20,20];
Normalization=0;
isSpring=0;

load('condition1.mat');
X_train=Log_All_data(1:20900,[1,4,8]);%setdiff(1:12,7));
Y_train=Log_All_data(1:20900,7);
[~,l]=size(X_train);
dataset_label='NucSize';
if Normalization
    X_prcs=(X_train-mean(X_train))/sqrtm(cov(X_train));
    Norm_label='NewNorm';
else
    X_prcs=(X_train-mean(X_train))./std(X_train);
    Norm_label='OldNorm';
end
Y_prcs=(Y_train-mean(Y_train))./std(Y_train);

X = X_prcs;
Y = Y_prcs;

rng("default") % For reproducibility of the partition
c = cvpartition(length(Y),"Holdout",0.20);
trainingIdx = training(c); % Indices for the training set
XTrain = X(trainingIdx,:);
YTrain = Y(trainingIdx);
testIdx = test(c); % Indices for the test set
XTest = X(:,:);%X(testIdx,:);
YTest = Y(:,:);%Y(testIdx);

Mdl = fitrnet(XTrain,YTrain,"Standardize",true, ...
    "LayerSizes",node);

testPredictions = predict(Mdl,XTest);

%testPredictions = testPredictions*std(All_data(:,4))+mean(All_data(:,4));

% plot(YTest*std(All_data(:,4))+mean(All_data(:,4)),testPredictions,".")
% hold on
% plot(YTest*std(All_data(:,4))+mean(All_data(:,4)),YTest*std(All_data(:,4))+mean(All_data(:,4)),'linewidth',2)
% hold off

if isSpring
    X_ori=X_train;
    Y_ori=Y_train;
    Y_predict_ori=testPredictions;
else
    X_ori=exp(X_train);
    Y_ori=exp(Y_train);
    Y_predict_ori=testPredictions*std(Y_train)+mean(Y_train);
    Y_predict_ori=exp(Y_predict_ori);
end

scatter(Y_ori,Y_predict_ori,".")
hold on
plot(Y_ori,Y_ori,'linewidth',2)
hold off

xlabel("True")
ylabel("Predicted")
set(gca,'FontSize',18,'Fontname', 'Arial','linewidth',1)
xlim([0,3e5]);
ylim([0,3e5]);

% if ~isSpring
%     ylim([0,800]);
% end

err=abs(Y_predict_ori-Y_ori)./Y_ori;
mean(err)