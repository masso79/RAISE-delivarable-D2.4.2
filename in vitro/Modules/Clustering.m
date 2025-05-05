%% Clusetring IFR distribution
clear all
clc

path = uigetdir(pwd, 'Select the IFR folder');

dirInfo= dir(path);
files=~[dirInfo.isdir];
fileNames={dirInfo(files).name};

IFR = [];
numSamples = 0;
numExp = 0;

for k = 1:length(fileNames)
    load(fullfile(path,fileNames{k}));
    numSamples = max (size(eval(fileNames{k}(1:end-4)),2), numSamples);
    numExp = size(eval(fileNames{k}(1:end-4)),1) + numExp;
end

IFR = zeros(numExp,numSamples);
start = 1;

for k = 1:length(fileNames)
    IFR(start:size(eval(fileNames{k}(1:end-4)),1)+start-1, 1:size(eval(fileNames{k}(1:end-4)),2)) = eval(fileNames{k}(1:end-4));
    start = size(eval(fileNames{k}(1:end-4)),1)+start;
end
p = [];
for j= 1:size(IFR,1)
    for k = j+1:size(IFR,1)
        [h,s] = kstest2(IFR(j,:), IFR(k,:));
        t(j,k) = h;
        t(k,j) = h;
        p(j,k) = s;
        p(k,j) = s;
    end
end
IFR2 = log(IFR);
IFR2./max(IFR2')';

tree = linkage(p,'average');
figure
[value, val, idx] = dendrogram(tree,0,'ColorThreshold', 0.65*max(tree(:,3)));
figure
imagesc(p)
Preshaped = p(idx,idx);

figure
imagesc(Preshaped)


Cortex = 1:26;
Striatum = 27:55;
Thalamus = 56:84;
