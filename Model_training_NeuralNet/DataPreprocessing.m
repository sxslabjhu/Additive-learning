%% Data preprocessing and save the processed data
clc;clear;
Datatype=3;
RatioTraining=0.8; 
if Datatype==1
    %Generate_data_SpringNetwork;
    %load('Spring_exp_data.mat');
    load('XXX.mat');
    randidx=randperm(length(dX_training));
    X_all=dF_training(randidx,:);
    Y_all=dX_training(randidx,:);
    X_train_prcs=X_all(1:round(length(Y_all)*RatioTraining),:);
    Y_train_prcs=Y_all(1:round(length(Y_all)*RatioTraining),:);
    X_test_prcs=X_all((round(length(Y_all)*RatioTraining)+1):end,:);
    Y_test_prcs=Y_all((round(length(Y_all)*RatioTraining)+1):end,:);
    save('Spring_prcs.mat',...
        'X_train_prcs','Y_train_prcs',...
        'Y_test_prcs','X_test_prcs',...
        'X_all','Y_all','randidx');
elseif Datatype==2
    load('AcH3_Complete_data_v2.mat');
    randidx=randperm(length(Log_All_data));
    X_all=Log_All_data(randidx,:);
    Y_all=Log_All_data(randidx,:);
    X_train=X_all(1:round(length(Y_all)*RatioTraining),:);
    Y_train=Y_all(1:round(length(Y_all)*RatioTraining),:);
    X_test=X_all((round(length(Y_all)*RatioTraining)+1):end,:);
    Y_test=Y_all((round(length(Y_all)*RatioTraining)+1):end,:);
    X_train_prcs=(X_train-mean(X_train))/sqrtm(cov(X_train));
    Y_train_prcs=(Y_train-mean(Y_train))./std(Y_train);
    X_test_prcs=(X_test-mean(X_train))/sqrtm(cov(X_train));
    Y_test_prcs=(Y_test-mean(Y_train))./std(Y_train);
    save('AcH3_prcs.mat',...
        'X_train_prcs','Y_train_prcs',...
        'Y_test_prcs','X_test_prcs',...
        'X_all','Y_all','randidx');
elseif Datatype==3
    load('SenescenceMxIFData_230224.mat');
    cell_condition=1; % 1: control; 2. quiescent 3. 50 um Bleo 4. 250nm Doxo
    inputnum=[2,4:6,8:11];
    outputnum=7;
    wellid_temp = (find(wellidinfo==cell_condition));
    temp=ismember(wellidout,wellid_temp);
    X_raw=Mintout(temp,inputnum);
    Y_raw=Mintout(temp,outputnum);
    All_raw=[X_raw,Y_raw];
    All_rmoutlier=rmoutliers(All_raw,'gesd');
    randidx=randperm(length(All_rmoutlier));

    X_all=All_rmoutlier(randidx,1:(end-1));
    Y_all=All_rmoutlier(randidx,end);
    X_train=X_all(1:round(length(Y_all)*RatioTraining),:);
    Y_train=Y_all(1:round(length(Y_all)*RatioTraining),:);
    X_test=X_all((round(length(Y_all)*RatioTraining)+1):end,:);
    Y_test=Y_all((round(length(Y_all)*RatioTraining)+1):end,:);
    X_train_prcs=(X_train-mean(X_train))/sqrtm(cov(X_train));
    Y_train_prcs=(Y_train-mean(Y_train))./std(Y_train);
    X_test_prcs=(X_test-mean(X_train))/sqrtm(cov(X_train));
    Y_test_prcs=(Y_test-mean(Y_train))./std(Y_train);
    save(strcat('Senescence_cellcondition',...
        num2str(cell_condition),'_prcs.mat'));
end