% multiRun(X_data_male, y_data_male, X_data_female, y_data_female, x_names, y_names, numComponents, nRuns)
%load data from csv as table GTEx version 8

T= readtable('Male+Female.csv'); %Male+Female.csv contains both Male and Female samples with rows(samples) and columns(genes).
%T should be (the number of samples by how ever many input and output features)
%T = (n_samples x n_features)
%%
sex = 'male';
path = ''; %path for saving files
vdir = 1; %Param to flip vertically
hdir = -1; %Param to flip horizontally
inlim = [-3 3];
outlim = [-4 4];

data = T(:,:);

%Set input parameters
x_names = {'HDAC1','HDAC2','HDAC3','HDAC4','HDAC5','HDAC6','HDAC7','HDAC8','HDAC9','HDAC10','HDAC11',...
    'SIRT1','SIRT2','SIRT3','SIRT4','SIRT5','SIRT6','SIRT7',...
    'KAT2A', 'KAT2B', 'HAT1', 'ATF2', 'KAT5', 'KAT6A', 'KAT6B', 'KAT7', 'EP300', 'CREBBP', 'NCOA1', 'NCOA3', 'TAF1', 'GTF3C1', 'CLOCK'};
%'FOXO1', 'FOXO3', 'GATA4', 'GATA6', 'HIF1A', 'TRIM28', 'KLF4', 'KLF5', 'MEF2A', 'NFAT5', 'NFKB1', 'NKX25', 'NOTCH1', 'RUNX1', 'SHMT2', 'SOD1', 'TBX5', 'TGFB1', 'YY1'};
X_data = table2array(data(:,x_names));
%X_data = (n_samples x n_inputs)

%
y_names = {'SCN5A', 'CACNA1C', 'KCNH2', 'KCNQ1', 'KCNJ2', 'SLC8A1', 'ATP2A2', 'RYR2', 'GJA1'};
%'FOXO1', 'FOXO3', 'GATA4', 'GATA6', 'HIF1A', 'TRIM28', 'KLF4', 'KLF5', 'MEF2A', 'NFAT5', 'NFKB1', 'NKX25', 'NOTCH1', 'RUNX1', 'SHMT2', 'SOD1', 'TBX5', 'TGFB1', 'YY1'};

%extract features to be predicted
y_data = table2array(data(:,y_names));

%y_data = (n_samples x n_outputs)

numComponents = 4; %number of desired components to use in plsregression

%Extract Sex data from csv file
Sex = table2array(data(:,{'Sex'}));

X_data_male = X_data(Sex == 0,:);
y_data_male = y_data(Sex == 0,:);

X_data_female = X_data(Sex == 1,:);
y_data_female = y_data(Sex == 1,:);

if strcmp(sex,'male')
    X_data_n = X_data_male;
    y_data_n = y_data_male;
    color = [0.2, 0.36, 0.506];
    fprintf('male has been selected\n');
else
    X_data_n = X_data_female;
    y_data_n = y_data_female;
    color = [197/255, 90/255, 17/255];
    fprintf('female has been selected\n');
end

%% Correlation Matrix Generation
corrMat = corrcoef([X_data, y_data]);

%% Call Function
[BETAarr,avgBETA, All_error, Sum_errArray, Sum_errArray_abs, PLSData] = multiRun(X_data_n, y_data_n, x_names, y_names, numComponents, 1);

%%
% plot and save average

figure;

%create biplots
biplot_1([hdir*PLSData.avgXL(:,1),vdir*PLSData.avgXL(:,2)], x_names,color);
%title('Average Inputs Loading and Scores ')
daspect([1 1 1]);
xlabel(['PC1 ',num2str(round(PLSData.avgPCTVAR(1,1)*100),2),'%']);
ylabel(['PC2 ',num2str(round(PLSData.avgPCTVAR(1,2)*100),2),'%']);
xlim(inlim);
ylim(inlim);
line([inlim NaN 0 0],[0 0 NaN inlim], 'Color','black');
savefig([path,'avgInputBiplot_',sex,'.fig']);

figure;
ax2 = subplot(1,1,1);
biplot_1([hdir*PLSData.avgYL(:,1),vdir*PLSData.avgYL(:,2)],y_names, color);
%title(['Average Outputs Loading and Scores '])
daspect([1 1 1]);
xlabel(['PC1 ',num2str(round(PLSData.avgPCTVAR(2,1)*100),2),'%']);
ylabel(['PC2 ',num2str(round(PLSData.avgPCTVAR(2,2)*100),2),'%']);
xlim(outlim);
ylim(outlim);
line([outlim NaN 0 0],[0 0 NaN outlim], 'Color', 'black');
savefig([path,'avgOutputBiplot_',sex,'.fig']);

%creates scores plots
samples = string(1:length(PLSData.avgXS(:,:)));
figure;
sgtitle('Average X Scores vs Y Scores ');
subplot(2, 2, 1);
plot(PLSData.avgXS(:,1), PLSData.avgYS(:,1), 'b.');
text(PLSData.avgXS(:,1), PLSData.avgYS(:,1),samples);
xlabel(['X Component 1 (Percent Explained ',num2str(round(PLSData.avgPCTVAR(1,1)*100),2),'%)']);
ylabel(['Y Component 1 (Percent Explained ',num2str(round(PLSData.avgPCTVAR(2,1)*100),2),'%)']);

subplot(2, 2, 2);
plot(PLSData.avgXS(:,2), PLSData.avgYS(:,2), 'b.');
text(PLSData.avgXS(:,2), PLSData.avgYS(:,2),samples);
xlabel(['X Component 2 (Percent Explained ',num2str(round(PLSData.avgPCTVAR(1,2)*100),2),'%)']);
ylabel(['Y Component 2 (Percent Explained ',num2str(round(PLSData.avgPCTVAR(2,2)*100),2),'%)']);

subplot(2, 2, 3);
plot(PLSData.avgXS(:,3), PLSData.avgYS(:,3), 'b.');
text(PLSData.avgXS(:,3), PLSData.avgYS(:,3),samples);
xlabel(['X Component 3 (Percent Explained ',num2str(round(PLSData.avgPCTVAR(1,3)*100),2),'%)']);
ylabel(['Y Component 3 (Percent Explained ',num2str(round(PLSData.avgPCTVAR(2,3)*100),2),'%)']);

subplot(2, 2, 4);
plot(PLSData.avgXS(:,4), PLSData.avgYS(:,4), 'b.');
text(PLSData.avgXS(:,4), PLSData.avgYS(:,4),samples);
xlabel(['X Component 4 (Percent Explained ',num2str(round(PLSData.avgPCTVAR(1,4)*100),2),'%)']);
ylabel(['Y Component 4 (Percent Explained ',num2str(round(PLSData.avgPCTVAR(2,4)*100),2),'%)']);
savefig([path,'avgXYScores_',sex,'.fig']);
%

%%
% save the last 6 runs to compare to average
for i = 1:6
    figure;
    
    %create biplots
    biplot_1(PLSData.XL(:,1:2,i), PLSData.XS(:,1:2,i),x_names, color);
    title(['Inputs Loading and Scores ', num2str(i)])
    xlabel(['PC1 ',num2str(round(PLSData.PCTVAR(1,1,i)*100),2),'%)']);
    ylabel(['PC2 ',num2str(round(PLSData.PCTVAR(1,2,i)*100),2),'%)']);
    savefig([path,'InputBiplot_',sex,'_', num2str(i),'.fig']);
    
    figure;
    biplot_1(PLSData.YL(:,1:2,i), PLSData.YS(:,1:2,i),y_names, color);
    title(['Outputs Loading and Scores ', num2str(i)])
    xlabel(['PC1 ',num2str(round(PLSData.PCTVAR(2,1,i)*100),2),'%)']);
    ylabel(['PC2 ',num2str(round(PLSData.PCTVAR(2,2,i)*100),2),'%)']);
    savefig([path,'OutputBiplot_',sex,'_', num2str(i),'.fig']);
    
    %creates scores plots
    samples = string(1:length(PLSData.XS(:,:,i)));
    figure;
    sgtitle(['X Scores vs Y Scores ', num2str(i)]);
    subplot(2, 2, 1);
    plot(PLSData.XS(:,1,i), PLSData.YS(:,1,i), 'b.');
    text(PLSData.XS(:,1,i), PLSData.YS(:,1,i),samples);
    xlabel(['X Component 1 (Percent Explained ',num2str(round(PLSData.PCTVAR(1,1,i)*100),2),'%)']);
    ylabel(['Y Component 1 (Percent Explained ',num2str(round(PLSData.PCTVAR(2,1,i)*100),2),'%)']);
    
    subplot(2, 2, 2);
    plot(PLSData.XS(:,2,i), PLSData.YS(:,2,i), 'b.');
    text(PLSData.XS(:,2,i), PLSData.YS(:,2,i),samples);
    xlabel(['X Component 2 (Percent Explained ',num2str(round(PLSData.PCTVAR(1,2,i)*100),2),'%)']);
    ylabel(['Y Component 2 (Percent Explained ',num2str(round(PLSData.PCTVAR(2,2,i)*100),2),'%)']);
    
    subplot(2, 2, 3);
    plot(PLSData.XS(:,3,i), PLSData.YS(:,3,i), 'b.');
    text(PLSData.XS(:,3,i), PLSData.YS(:,3,i),samples);
    xlabel(['X Component 3 (Percent Explained ',num2str(round(PLSData.PCTVAR(1,3,i)*100),2),'%)']);
    ylabel(['Y Component 3 (Percent Explained ',num2str(round(PLSData.PCTVAR(2,3,i)*100),2),'%)']);
    
    subplot(2, 2, 4);
    plot(PLSData.XS(:,4,i), PLSData.YS(:,4,i), 'b.');
    text(PLSData.XS(:,4,i), PLSData.YS(:,4,i),samples);
    xlabel(['X Component 4 (Percent Explained ',num2str(round(PLSData.PCTVAR(1,4,i)*100),2),'%)']);
    ylabel(['Y Component 4 (Percent Explained ',num2str(round(PLSData.PCTVAR(2,4,i)*100),2),'%)']);
    savefig([path,'XYScores_',sex,'_',num2str(i),'.fig']);
    
end

%%
figure;

plot([1:numComponents], cumsum(PLSData.avgPCTVAR(1,:)), 'bo-');
title('Percent Variance Explained X');
xlabel('Num Components');
ylabel('Percent Variance Explained X');
xticks([1:numComponents]);
ylim([0,1]);

figure;

plot([1:numComponents], cumsum(PLSData.avgPCTVAR(2,:)), 'bo-');
title('Percent Variance Explained Y');
xlabel('Num Components');
ylabel('Percent Variance Explained Y');
xticks([1:numComponents]);
ylim([0,1]);
%%
function[BETAarr,avgBETA, All_error, Sum_errArray, Sum_errArray_abs, PLSData] = multiRun(X_data, y_data, x_names, y_names, numComponents, nRuns)

%allocate space for all arrarys
BETAarr = zeros(length(x_names)+1, length(y_names), nRuns);
All_error = zeros(floor(length(X_data)*.1),length(y_names),nRuns);
Sum_errArray = zeros(nRuns,length(y_names));
Sum_errArray_abs = zeros(nRuns,length(y_names));

%Object that contains most of the PLS data. Allocating all space.
PLSData.XL = zeros(length(x_names),numComponents, nRuns);
PLSData.YL = zeros(length(y_names), numComponents, nRuns);
PLSData.XS = zeros(ceil(length(X_data)*.9),numComponents, nRuns);
PLSData.YS = zeros(ceil(length(X_data)*.9), numComponents, nRuns);
PLSData.PCTVAR = zeros(2, numComponents, nRuns);
PLSData.MSE = zeros(2, numComponents+1,nRuns);
PLSData.R2 = zeros(nRuns,1);
PLSData.RMSE = zeros(nRuns,1);

cv = cvpartition(size(X_data,1),'HoldOut',0.1); %0.1 is the percent holdout.
idx = cv.test;
TestSize = cv.TestSize;

PLSData.y_pred = zeros(TestSize, length(y_names));
PLSData.y_test = zeros(TestSize, length(y_names));

%for loop to run nRun tiems
for i = 1:nRuns
    
    %partition train and test data
    cv = cvpartition(size(X_data,1),'HoldOut',0.1); %0.1 is the percent holdout.
    idx = cv.test;
    
    X_train = X_data(~idx,:);
    X_test = X_data(idx,:);
    y_train = y_data(~idx,:);
    y_test = y_data(idx,:);
    
    %PLS
    [XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(X_train, y_train, numComponents);
    
    %Store data for individual run
    PLSData.XL(:,:,i) = XL;
    PLSData.YL(:,:,i) = YL;
    PLSData.XS(:,:,i) = XS;
    PLSData.YS(:,:,i) = YS;
    PLSData.PCTVAR(:,:,i) = PCTVAR;
    PLSData.MSE(:,:,i) = MSE;
    
    BETAarr(:,:,i) = BETA;
    
    %Prediction Calculation
    y_pred = [ones(size(X_test,1),1) X_test]*BETA;
    
    PLSData.y_pred(i*TestSize-(TestSize-1):i*TestSize,:) = y_pred;
    PLSData.y_test(i*TestSize-(TestSize-1):i*TestSize,:) = y_test;
    
    
    %calculate male error
    error_abs = abs(y_pred - y_test);
    error = y_pred - y_test;
    Sum_error_abs = sum(error_abs, 1);
    Sum_error = sum(error, 1);
    
    %store all the error
    All_error(:,:,i) = error;
    
    %store summed error array
    Sum_errArray(i,:) = Sum_error;
    Sum_errArray_abs(i,:) = Sum_error_abs;
    
    %Mean squared error calculation
    mse = immse(y_test, y_pred);
    RMSE = sqrt(mse);
    
    
    TSS = sum((y_train-mean(y_train)).^2); %Total Sum of Squares
    RSS_PLS = sum(stats.Yresiduals.^2); %Total sum of squares of residuals
    R2 = 1 - RSS_PLS/TSS; %How well are model explains the output
    
    PLSData.RMSE(i) = RMSE;
    PLSData.R2(i) = R2;
    
end
avgBETA = mean(BETAarr,3);
PLSData.avgXS = mean(XS,3);
PLSData.avgYS = mean(YS,3);
PLSData.avgXL = mean(XL,3);
PLSData.avgYL = mean(YL,3);
PLSData.avgMSE = mean(PLSData.MSE,3);
PLSData.avgPCTVAR = mean(PLSData.PCTVAR,3);
PLSData.avgR2 = mean(PLSData.R2);
PLSData.avgRMSE = mean(PLSData.RMSE);
end

