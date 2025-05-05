clear all
clc

%% IFR Analysis
% Select the IFR files
[path_IFR] = uigetdir (pwd,'Select the IFR folder');

% Select the Network Burst file
[nb_file, path_nb] = uigetfile('.mat','Select the NetworkBurstDetectionFiles');

% Load the IFR files
cd(path_IFR);
dir_IFR = dir;
count = 0;
for i = 3:length(dir_IFR)
    if strcmp(dir_IFR(i).name(end-17:end),'saveParameters.mat') == 0
        count =count+1;
        load(dir_IFR(i).name);
        name = split(dir_IFR(i).name,'.');
        name = split(name{1},'_');
        name_clusters(count) =  string(name{end});
        IFR_clusters (count,:) = cumIFR;
    end
end

% Load the nb file
load(fullfile(path_nb,nb_file));
%%
% Delete IFR of the electrodes covered by the rood
index_gray = find(strcmp(name_clusters, 'Gray'));
IFR_clusters(index_gray,:) = [];
name_clusters(index_gray) = [];

fs = 10000;
% Inter network burst interval
diffSample = diff(netBursts(:,1));
diff_ms = diffSample/fs*1000;

% For each nb, we identify the activation sequence

for i = 1:length(netBursts)
    tmp_clusters = name_clusters;
    chunck = full(IFR_clusters(:,netBursts(i,1):netBursts(i,1)+netBursts(i,4)));

    % to find the maximum values we need to smooth the curves
    % We consider a peak if it has an heigth of at least 70% of the maximum
    % value
    count = 1;
    for j = 1:size(chunck,1)
        smooth_chunck(j,:) = smooth(chunck(j,:),2000);
        [values{j}, index{j}] = findpeaks(smooth_chunck(j,:),'MinPeakHeight', 0.05*max(smooth_chunck(j,:)), 'MinPeakDistance', netBursts(i,4)-10);
        % [values{j}, index{j}] = findpeaks(smooth_chunck(j,:),'MinPeakHeight', 0.7*max(smooth_chunck(j,:)),'MinPeakDistance',300,'MinPeakProminence',0.1*max(smooth_chunck(j,:)));
        for k = 1:length(values{j})
            clust_act{count} = tmp_clusters{j};
            count = count +1;
        end
    end
    samples = cell2mat(index);
    [r,c] = sort(samples);
    sequences{i} = clust_act(c);
    clear smooth_chunck
    clear clust_act
end

%% Count the percentage of occurrences
n = max(cellfun('size',sequences,2));
temporal_sequence = string(zeros(length(sequences),n));
for i = 1:length(sequences)
    temporal_sequence(i,1:length(sequences{i})) = convertCharsToStrings(sequences{i});
end

temporal_sequence(strcmp(temporal_sequence,'0'))="";
temporal_sequence(ismissing(temporal_sequence,'0'))="";
seq = sortrows(temporal_sequence);
count = 1;
for i = 1:length(seq)
    if i == 1
        name = "";
        for k = 1: size(seq,2)
            name = append(name, seq(i,k));
        end
        ripet(count) = 1;
        seq_uni(count) = name;
        index_uni(count) = i;
    else
        toCon = "";
        for k = 1: size(seq,2)
            toCon = append(toCon, seq(i,k));
        end
        if ~strcmp(toCon, seq_uni(count))
            count = count +1;
            seq_uni(count) = toCon;
            index_uni (count)= i;
            ripet(count) = 1;
        else
            ripet(count) = ripet(count)+1;
        end    
    end
end

percentage = ripet./sum(ripet)*100;
sequence_ripet = seq(index_uni,:);


%% Save results
cd ..
phase = split(path_IFR,'_');
phase = phase{end};
xlswrite(strcat('IFRAnalysis_', phase),temporal_sequence,'Temporal_Activation_Sequences')
path = pwd;
newExcel = actxserver('excel.application');
excelWB = newExcel.Workbooks.Open(strcat(pwd,strcat('\IFRAnalysis_', phase)),0,false);
newExcel.Sheets.Item(1).Delete;
excelWB.Save();
excelWB.Close();
newExcel.Quit();
delete(newExcel);

% xlswrite(strcat('IFRAnalysis_', phase),delay,'Delay_Activation_ms')
xlswrite(strcat('IFRAnalysis_', phase),sequence_ripet,'Activation_Sequences','A1')
xlswrite(strcat('IFRAnalysis_', phase),percentage','Activation_Sequences','I1')
xlswrite(strcat('IFRAnalysis_', phase),diff_ms,'InterNetworkBurstInterval_ms')



