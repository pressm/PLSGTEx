% D = Ashamn Coefficients
% S =  Binomial Seperration Coefficient
% BF = Bimodality Coefficient
% um = mean of male
% uf = mean of female
% sm = standard deviation male
% sf = standard deviation female
%
data= readtable('Male+Female.csv');

y_names = {'SCN5A', 'CACNA1C', 'KCNH2', 'KCNQ1', 'KCNJ2', 'ATP1A1', 'SLC8A1', 'ATP2A2', 'RYR2', 'GJA1'};

Sex = data.Sex;
%%
file = load('MaleFemaleErrorRun1.mat');
MF_Error = file.T_All_error;

D = zeros(width(MF_Error)-1,1);
S = zeros(width(MF_Error)-1,1);
beta = zeros(width(MF_Error)-1,1);
BF = zeros(width(MF_Error)-1,2);
uf = zeros(width(MF_Error)-1,1);
um = zeros(width(MF_Error)-1,1);
sf = zeros(width(MF_Error)-1,1);
sm = zeros(width(MF_Error)-1,1);
for i = 1:width(MF_Error)-1
    [pdca, gn, ~] = fitdist(MF_Error.(i),'Normal', 'By',MF_Error.Sex);
    um(i) = pdca{1}.mean;
    uf(i) = pdca{2}.mean;
    sm(i) = pdca{1}.sigma;
    sf(i) = pdca{2}.sigma;
    D(i) = ashman(uf(i),um(i), sf(i), sm(i));
    S(i) = bisep(uf(i), um(i), sf(i), sm(i));
    beta(i) = bicoeff(MF_Error.(i));
    [BF(i,1),BF(i,2)] = bimodalitycoeff(MF_Error.(i));
end
MF_T = table(uf,um,sf,sm);

%%
file = load('Random1Runs1.mat');
RandomError = file.T_All_error;

D_R = zeros(width(RandomError)-1,1);
S_R = zeros(width(RandomError)-1,1);
beta_R = zeros(width(RandomError)-1,1);
BF_R = zeros(width(RandomError)-1,2);
u1 = zeros(width(RandomError)-1,1);
u2 = zeros(width(RandomError)-1,1);
s1 = zeros(width(RandomError)-1,1);
s2 = zeros(width(RandomError)-1,1);
for i = 1:width(RandomError)-1
    [pdca, gr, ~] = fitdist(RandomError.(i),'Normal', 'By',RandomError.Sex);
    u1(i) = pdca{1}.mean;
    u2(i) = pdca{2}.mean;
    s1(i) = pdca{1}.sigma;
    s2(i) = pdca{2}.sigma;
    D_R(i) = ashman(u1(i),u2(i), s1(i), s2(i));
    S_R(i) = bisep(u1(i), u2(i), s1(i), s2(i));
    beta_R(i)= bicoeff(RandomError.(i));
    [BF_R(i,1),BF_R(i,2)] = bimodalitycoeff(RandomError.(i));
end
Random_T = table(u1,u2,s1,s2);
%%
function D = ashman(u1, u2, s1, s2)
D = (2^(1/2))*abs(u1-u2)/sqrt(s1^2 + s2^2);
end

function S = bisep(u1, u2, s1, s2)
S = abs(u1-u2)/(2*(s1+s2));
end

function beta = bicoeff(X)
skew = skewness(X);
kurt = kurtosis(X);
beta = (skew^2 + 1)/kurt;
end
