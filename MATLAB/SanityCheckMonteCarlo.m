%load data from csv as table GTEx version 8

T= readtable('Male+Female.csv');
%T should be (the number of samples by how ever many input and output features)
%T = (n_samples x n_features)

data = T(:,:);
%Extract ID feature
ID = table2array(data(:, {'ID'}));

%Set input parameters
x_names = {'HDAC1','HDAC2','HDAC3','HDAC4','HDAC5','HDAC6','HDAC7','HDAC8','HDAC9','HDAC10','HDAC11',...
    'SIRT1','SIRT2','SIRT3','SIRT4','SIRT5','SIRT6','SIRT7',...
    'KAT2A', 'KAT2B', 'HAT1', 'ATF2', 'KAT5', 'KAT6A', 'KAT6B', 'KAT7', 'EP300', 'CREBBP', 'NCOA1', 'NCOA3', 'TAF1', 'GTF3C1', 'CLOCK'};
%'FOXO1', 'FOXO3', 'GATA4', 'GATA6', 'HIF1A', 'TRIM28', 'KLF4', 'KLF5', 'MEF2A', 'NFAT5', 'NFKB1', 'NKX25', 'NOTCH1', 'RUNX1', 'SHMT2', 'SOD1', 'TBX5', 'TGFB1', 'YY1'};

%     {'HDAC1','HDAC2','HDAC3','HDAC4','HDAC5','HDAC6','HDAC7','HDAC8','HDAC9','HDAC10','HDAC11',...
%     'SIRT1','SIRT2','SIRT3','SIRT4','SIRT5','SIRT6','SIRT7',...
%     'KAT2A', 'KAT2B', 'HAT1', 'ATF2', 'KAT5', 'KAT6A', 'KAT6B', 'KAT7', 'EP300', 'CREBBP', 'NCOA1', 'NCOA3', 'TAF1', 'GTF3C1', 'CLOCK'};
%'FOXO1', 'FOXO3', 'GATA4', 'GATA6', 'HIF1A', 'TRIM28', 'KLF4', 'KLF5', 'MEF2A', 'NFAT5', 'NFKB1', 'NKX25', 'NOTCH1', 'RUNX1', 'SHMT2', 'SOD1', 'TBX5', 'TGFB1', 'YY1'};
X_data = table2array(data(:,x_names));
%X_data = (n_samples x n_inputs)

%
y_names = {'SCN5A', 'CACNA1C', 'KCNH2', 'KCNQ1', 'KCNJ2', 'ATP1A1', 'SLC8A1', 'ATP2A2', 'RYR2', 'GJA1'};
%'FOXO1', 'FOXO3', 'GATA4', 'GATA6', 'HIF1A', 'TRIM28', 'KLF4', 'KLF5', 'MEF2A', 'NFAT5', 'NFKB1', 'NKX25', 'NOTCH1', 'RUNX1', 'SHMT2', 'SOD1', 'TBX5', 'TGFB1', 'YY1'};

%extract features to be predicted
y_data = table2array(data(:,y_names));

%y_data = (n_samples x n_outputs)

numComponents = 4; %number of desired components to use in plsregression
color1 = '#2F5597'; %Male Female Colors
color2 = '#C55A11';
%color1 = '#DA70D6';% Random Mix Colors
%color2 = '#9370DB';

%Extract Sex data from csv file
Sex = table2array(data(:,{'Sex'}));

%seperate into male and female data
X_data_male = X_data(Sex == 0,:);
y_data_male = y_data(Sex == 0,:);

X_data_female = X_data(Sex == 1,:);
y_data_female = y_data(Sex == 1,:);

%% Call Function
n = 1;
All_error_male = [];
All_error_female = [];
parfor i = 1:n
    %% Grab 84 random male samples
    cv = cvpartition(size(X_data_male,1), 'Holdout', 0.535);
    idx_84 = cv.test;
    
    X_data_84 = X_data_male(idx_84,:);
    y_data_84 = y_data_male(idx_84,:);
    
    %split samples in half and mix
    cv = cvpartition(size(X_data_84,1), 'Holdout', 0.5);
    idx_m = cv.test;
    
    cv = cvpartition(size(X_data_female,1), 'Holdout', 0.5);
    idx_f = cv.test;
    
    X_data_1 = [X_data_84(~idx_m,:); X_data_female(~idx_f,:)];
    X_data_2 = [X_data_84(idx_m,:); X_data_female(idx_f,:)];
    y_data_1 = [y_data_84(~idx_m,:); y_data_female(~idx_f,:)];
    y_data_2 = [y_data_84(idx_m,:); y_data_female(idx_f,:)];
    
    
    tic
    [BETA_pls1, All_error_male_pls1, All_error_female_pls1, absErrArray_male, absErrArray_female, PLSDataS] = SanityCheck(X_data_84,y_data_84, X_data_female, y_data_female, x_names, y_names, numComponents, 1000);
    toc
    
    All_error_male = [All_error_male; All_error_male_pls1];
    All_error_female = [All_error_female; All_error_female_pls1];
end
%% Store information 
T_female = array2table(All_error_female, 'VariableNames', y_names);
female = string.empty(1000*n,0);
female(1:1000*n) = 'B';
T_female_sex = array2table(female', 'VariableNames', {'Sex'});
T_All_error_female = [T_female T_female_sex];

T_male = array2table(All_error_male, 'VariableNames', y_names);
male = string.empty(1000*n,0);
male(1:1000*n) = 'A';
T_male_sex = array2table(male', 'VariableNames', {'Sex'});
T_All_error_male = [T_male T_male_sex];

T_All_error = [T_All_error_male; T_All_error_female];
%%
figure(1);
sTable = stack(T_All_error, 1:10,'NewDataVariableName', 'values', 'IndexVariableName','Channels');
t = boxchart(sTable.Channels, sTable.values, 'GroupByColor', sTable.Sex);
t(2).BoxFaceColor = color2;
t(2).MarkerColor = color2;
t(1).BoxFaceColor = color1;
t(1).MarkerColor = color1;
legend('Male', 'Female');
hold on
yline(0,'k--', 'LineWidth', 2)
hold off

%%
figure(2);
for i = 1:length(y_names)
    subplot(5,2,i);
    
    histogram(All_error_female(:,i),'BinWidth',.01,'FaceColor',color2)
    hold on;
    histogram(All_error_male(:,i),'BinWidth',.01,'FaceColor',color1)
    hold off;
    title(y_names(i));
    xlim([-.15, .15001])
    
end
%%
function[BETA1, All_error_male, All_error_female, absErrArray_male, absErrArray_female,PLSDataS] = SanityCheck(X_data_male, y_data_male, X_data_female, y_data_female, x_names, y_names, numComponents, nRuns)

z = floor(length(X_data_male)*.1);
k = round(z/length(X_data_female), 2);

%allocate space for all arrarys
BETA1 = zeros(length(x_names)+1, length(y_names), nRuns);
absErrArray_male = zeros(nRuns,length(y_names));
absErrArray_female = zeros(nRuns,length(y_names));
All_error_male = zeros(nRuns,length(y_names));
All_error_female = zeros(nRuns,length(y_names));

%Object that contains most of the PLS data. Allocating all space.
PLSDataS.XL = zeros(length(x_names),numComponents, nRuns);
PLSDataS.YL = zeros(length(y_names), numComponents, nRuns);
PLSDataS.XS = zeros(ceil(length(X_data_male)*.9),numComponents, nRuns);
PLSDataS.YS = zeros(ceil(length(X_data_male)*.9), numComponents, nRuns);
PLSDataS.PCTVAR = zeros(2, numComponents, nRuns);
PLSDataS.MSE = zeros(2, numComponents+1,nRuns);
PLSDataS.R2 = zeros(nRuns,1);
PLSDataS.RMSE_male = zeros(nRuns,1);
PLSDataS.RMSE_female = zeros(nRuns,1);
PLSDataS.y_pred_female = zeros(floor(length(y_data_female)*k), length(y_names), nRuns);
PLSDataS.y_pred_male = zeros(floor(length(y_data_male)*0.1), length(y_names), nRuns);
PLSDataS.y_test_female = zeros(floor(length(y_data_female)*k), length(y_names), nRuns);
PLSDataS.RPDp_male = zeros(0,nRuns);
PLSDataS.RPDp_female = zeros(0,nRuns);

%for loop to run nRun tiems
for i = 1:nRuns
    
    %partition train and test data
    cv = cvpartition(size(X_data_male,1),'HoldOut', 0.1); %0.1 is the percent holdout.
    idx_male = cv.test;
    
    cv = cvpartition(size(X_data_female,1),'HoldOut', k);
    idx_female = cv.test;
    
    X_test_female = X_data_female(idx_female,:);
    y_train_female = y_data_female(~idx_female,:);
    y_test_female = y_data_female(idx_female,:);
    
    
    PLSDataS.y_test_female(:,:,i) = y_test_female;
    
    X_train_male = X_data_male(~idx_male,:);
    X_test_male = X_data_male(idx_male,:);
    y_train_male = y_data_male(~idx_male,:);
    y_test_male = y_data_male(idx_male,:);
    
    %PLS
    [XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(X_train_male, y_train_male, numComponents);
    
    %Store data for individual run
    PLSDataS.XL(:,:,i) = XL;
    PLSDataS.YL(:,:,i) = YL;
    PLSDataS.XS(:,:,i) = XS;
    PLSDataS.YS(:,:,i) = YS;
    PLSDataS.PCTVAR(:,:,i) = PCTVAR;
    PLSDataS.MSE(:,:,i) = MSE;
    
    
    BETA1(:,:,i) = BETA;
    
    %Prediction Calculation
    y_pred_male = [ones(size(X_test_male,1),1) X_test_male]*BETA;
    y_pred_female = [ones(size(X_test_female,1),1) X_test_female]*BETA;
    
    PLSDataS.y_pred_female(:,:,i) = y_pred_female;
    PLSDataS.y_pred_male(:,:,i) = y_pred_male;
    
    
    %calculate male error
    absErrArray_male = mean(abs(y_pred_male - y_test_male));
    error_male = y_pred_male - y_test_male;
    
    All_error_male(i,:) = mean(error_male);
    
    %calculate female error
    absErrArray_female(i,:) = mean(abs(y_pred_female - y_test_female));
    error_female = y_pred_female - y_test_female;
    
    %RPDp
    PLSDataS.RPDp_female(i) = std(y_test_female)/sum((error_female).^2);
    PLSDataS.RPDp_male(i) = std(y_test_male)/sum((error_male).^2);
    
    %store all female error
    All_error_female(i,:) = mean(error_female);
    
    %Mean squared error calculation
    mse_male = immse(y_test_male, y_pred_male);
    RMSE_male = sqrt(mse_male);
    
    mse_female = immse(y_test_female, y_pred_female);
    RMSE_female = sqrt(mse_female);
    
    TSS = sum((y_train_male-mean(y_train_male)).^2); %Total Sum of Squares
    RSS_PLS = sum(stats.Yresiduals.^2); %Total sum of squares of residuals
    R2 = 1 - RSS_PLS/TSS; %How well are model explains the output
    
    PLSDataS.RMSE_male(i) = RMSE_male;
    PLSDataS.RMSE_female(i) = RMSE_female;
    PLSDataS.R2(i) = R2;
    
end
end
