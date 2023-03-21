% 检测有没有outlier，然后去csv的对应位置修改
clear
[FileName, FileAddress] = uigetfile('E:\HexagonData\Videos\20230223\example.csv');
csvraw = readmatrix([FileAddress, FileName]);
RawData(:, 1) = csvraw(:, 1) + 1;
RawData(:, 2:3) = csvraw(:, 2:3);
RawData(1, 4) = sqrt((RawData(1, 2) - RawData(2, 2)).^2+(RawData(1, 3) - RawData(2, 3)).^2);
for i = 2:length(RawData) - 1
    RawData(i, 4) = 0.5 * sqrt((RawData(i+1, 2) - RawData(i-1, 2)).^2+(RawData(i+1, 3) - RawData(i-1, 3)).^2);
end
RawData(end, 4) = sqrt((RawData(end, 2) - RawData(end-1, 2)).^2+(RawData(end, 3) - RawData(end, 3)).^2);
plot(RawData(:,4),'k.');
% figure;
% plot(RawData(:, 2),RawData(:, 3),'k-');
