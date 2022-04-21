%load normalized data from csv as table GTEx version 8
S = 1;
%choose either male or female
if S == 1
    T= readtable('norm-gtexm.csv','ReadRowNames',true, 'VariableNamingRule', 'preserve');    
elseif S==2
    T= readtable('norm-gtexf.csv','ReadRowNames',true, 'VariableNamingRule', 'preserve');
end
%%
genes = T.Properties.RowNames;
%%
inputs = {'HDAC1','HDAC2','HDAC3','HDAC4','HDAC5','HDAC6','HDAC7','HDAC8','HDAC9','HDAC10','HDAC11',...
    'SIRT1','SIRT2','SIRT3','SIRT4','SIRT5','SIRT6','SIRT7',...
    'KAT2A', 'KAT2B', 'HAT1', 'ATF2', 'KAT5', 'KAT6A', 'KAT6B', 'KAT7', 'EP300', 'CREBBP', 'NCOA1', 'NCOA3', 'TAF1', 'GTF3C1', 'CLOCK'};
%'FOXO1', 'FOXO3', 'GATA4', 'GATA6', 'HIF1A', 'TRIM28', 'KLF4', 'KLF5', 'MEF2A', 'NFAT5', 'NFKB1', 'NKX25', 'NOTCH1', 'RUNX1', 'SHMT2', 'SOD1', 'TBX5', 'TGFB1', 'YY1'};

%
outputs = {'SCN5A', 'CACNA1C', 'KCNH2', 'KCNQ1', 'KCNJ2', 'ATP1A1', 'SLC8A1', 'ATP2A2', 'RYR2', 'GJA1'};
%'FOXO1', 'FOXO3', 'GATA4', 'GATA6', 'HIF1A', 'TRIM28', 'KLF4', 'KLF5', 'MEF2A', 'NFAT5', 'NFKB1', 'NKX25', 'NOTCH1', 'RUNX1', 'SHMT2', 'SOD1', 'TBX5', 'TGFB1', 'YY1'};%%
%%
%%
dataT = T;
%%
dataT = rows2vars(dataT);
%%
X_data = table2array(dataT(:, inputs));
y_data = table2array(dataT(:, outputs));
%%
samples = string(1:length(X_data));

remove = [187, 161, 19];
X_data_n = X_data;
X_data_n(remove,:) = [];
y_data_n = y_data;
y_data_n(remove,:) = [];
samples(remove) = [];
samples_removedX = X_data(remove,:);
samples_removedY = y_data(remove,:);

%%

[XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(X_data_n, y_data_n, 4);

%%


%create biplots

%creates scores plots
figure;
daspect([1 1 1]);
%sgtitle(['X Scores vs Y Scores ']);
subplot(2, 2, 1);
plot(XS(:,1), YS(:,1), 'b.');
text(XS(:,1), YS(:,1),samples);
xlabel(['XS PC1 ',num2str(round(PCTVAR(1,1)*100),2),'%']);
ylabel(['YS PC1 ',num2str(round(PCTVAR(2,1)*100),2),'%']);

subplot(2, 2, 2);
plot(XS(:,2), YS(:,2), 'b.');
text(XS(:,2), YS(:,2),samples);
xlabel(['XS PC2 ',num2str(round(PCTVAR(1,2)*100),2),'%']);
ylabel(['YS PC2 ',num2str(round(PCTVAR(2,2)*100),2),'%']);

subplot(2, 2, 3);
plot(XS(:,3), YS(:,3), 'b.');
text(XS(:,3), YS(:,3),samples);
xlabel(['XS PC3 ',num2str(round(PCTVAR(1,3)*100),2),'%']);
ylabel(['YS PC3 ',num2str(round(PCTVAR(2,3)*100),2),'%']);

subplot(2, 2, 4);
plot(XS(:,4), YS(:,4), 'b.');
text(XS(:,4), YS(:,4),samples);
xlabel(['XS PC4 ',num2str(round(PCTVAR(1,4)*100),2),'%']);
ylabel(['YS PC4 ',num2str(round(PCTVAR(2,4)*100),2),'%']);