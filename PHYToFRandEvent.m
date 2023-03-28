%From SPKC to Event and FR
% DIO有3个，第一个是从start50到water50，第二第三个是选项2和选项3的start50到end
%% read NPY from phy
FileAddress=uigetdir('D:\TrodesData\','Select a PHYgui folder');
spike_clusters = readNPY([FileAddress,'\spike_clusters.npy'])+1;%所以分类指标都+1
spike_times = readNPY([FileAddress,'\spike_times.npy']);
fileID = fopen([FileAddress,'\cluster_group.tsv']);
X1 = textscan(fileID,'%s %s');
fclose(fileID);
for i = 1:length(X1{1})-1
% Nclust(i) = str2num(X1{1}{i+1});
Nclustgood(i) = (str2num(X1{1}{i+1})+1)*strcmp(X1{2}{i+1},'good');%所以分类指标都+1
end
Nclustgood(Nclustgood==0) = [];
SpikeTimeStamps = {};
for i = 1:length(Nclustgood)
    SpikeTimeStamps{end+1} = double(spike_times(spike_clusters==Nclustgood(i))); %ts直接为单位
end
%% trials
[FileName,FileAddress]=uigetfile('*.dat','Select a DIO file','E:\TrodesData\');
data1 = readTrodesExtractedDataFile([FileAddress,'\',FileName]);
TrialS = data1.fields(1).data;
OnTrialS = double(TrialS(2:2:end)-TrialS(1))/30; %ms为单位


for i = 1:length(OnTrialS)
OnTrialSmsStart = OnTrialS(i)-2500;
OnTrialSmsEnd = OnTrialS(i)+2500;
for j = 1:length(SpikeTimeStamps)
TST{i,j} = SpikeTimeStamps{1,j}(find((SpikeTimeStamps{1,j}>OnTrialSmsStart) & (SpikeTimeStamps{1,j}<OnTrialSmsEnd)))-OnTrialS(i);
end
end

for i = 1:size(TST,2)
FRN = zeros(size(TST,1),51);
for j = 1:size(TST,1)
FRN(j,:) = hist(TST{j,i},-2500:100:2500);
end
FR100{i} = FRN;
FR100m(:,i) = mean(FRN);
end


% clearvars -except SpikeTimeStamps OnTrialS TrialS TST FR100 FR100m

