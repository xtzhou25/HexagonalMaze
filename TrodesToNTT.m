clear
% 这玩意负责把spikegadgets的spikes转换成ntt文件交由MClust处理。权宜之计。
[FileName, FileAddress] = uigetfile('*.dat', 'Select a dat file', 'MultiSelect', 'on');
cd(FileAddress);
if isa(FileName, 'cell')
    nn = length(FileName);
else
    nn = 1;
    tempfile{1} = FileName;
    FileName = tempfile;
end

for i = 1:nn
    raw = readTrodesExtractedDataFile(FileName{i});
    
    ts = double(raw.fields(1).data)./30000;
    n=length(ts);
    
    wave1 = double(raw.fields(2).data);
    wave2 = double(raw.fields(3).data);
    wave3 = double(raw.fields(4).data);
    wave4 = double(raw.fields(5).data);
    
%     SampleWave =zeros(32,4,n);
%     SampleWave(:, 1, :) = resample(-wave1', 32, size(wave1, 2)); %老办法 降采样
%     可能引入误差
%     SampleWave(:, 2, :) = resample(-wave2', 32, size(wave1, 2));
%     SampleWave(:, 3, :) = resample(-wave3', 32, size(wave1, 2));
%     SampleWave(:, 4, :) = resample(-wave4', 32, size(wave1, 2));


    SampleWave =zeros(32,4,n);
    SampleWave(:, 1, :) = -wave1(:,3:34)';  % 掐头去尾保留原始数据点，干净的原始数据
    SampleWave(:, 2, :) = -wave2(:,3:34)';
    SampleWave(:, 3, :) = -wave3(:,3:34)';
    SampleWave(:, 4, :) = -wave4(:,3:34)';

    Mat2NlxSpike([FileName{i}(24:end),'.ntt'], 0, 1, [], [1 0 1 0 1 0], ts'*1000000,zeros(1,length(ts)), SampleWave);
end