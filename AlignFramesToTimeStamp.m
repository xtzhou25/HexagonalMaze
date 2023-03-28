% get framestamp and timestamp, align

%% Load Data
clear
datadate = '20230325';
MazeCenter = [362, 367];
% specify the date of recordings

% parameters
% misstrial = 0;
% % if the trialdivision not paired to DIO
% trialNEED = [1 42];
% % which trial do we need to analyse

[trialdivisionName, trialdivisionPath] = uigetfile('*.mat', 'Select a trialdivision file', ['E:\HexagonData\Videos\', datadate, '\']);
load([trialdivisionPath, trialdivisionName]);
%get trialdivision result:trialnum
LEDChange = [trialnum(:, 1); trialnum(end, 1)] - [trialnum(1, 1); trialnum(:, 1)];
trialFrame = find(LEDChange);
trialN = length(trialFrame);
% first frame of a trial


[dioName1, dioPath] = uigetfile('*.dat', 'Select first dio file', ['D:\TrodesData\', datadate, '\']);
dioName2 = [dioName1(1:41), '2.dat'];
dioraw1 = readTrodesExtractedDataFile([dioPath, dioName1]);
dioraw2 = readTrodesExtractedDataFile([dioPath, dioName2]);
%get dio raw

dioFirstTS = dioraw1.first_timestamp;
% all dio ts need to minus first_timestamp to align to spikeTS

%% liner model from Frame to TS

trialTS = double(dioraw1.fields(1).data(3:trialN+2)) - dioFirstTS;
% change it when the trials are not

% plot(trialFrame,trialTS,'rp');
linermdl = fitlm(trialFrame, trialTS);
% linermdl2 = fitlm(trialTS,trialFrame);
hold on
plot(linermdl);
hold off

% for this liner model, you can predict TS for each Frame

%% csv file

[csvName, csvPath] = uigetfile('*.csv', 'Select a track file', trialdivisionPath);
csvraw = readmatrix([csvPath, csvName]);

% RawData = zeros(length(csvraw), 6+length(neuron));

%% movement information


head(:, 1) = csvraw(:, 2) - MazeCenter(1);
head(:, 2) = -csvraw(:, 3) + MazeCenter(2);
headsmooth = smoothdata(head, 'gaussian', 10);
headsmooth = headsmooth ./ 5;
% 换算成厘米

speed = zeros(length(head), 1);
speed(1) = sqrt((headsmooth(1, 1) - headsmooth(2, 1)).^2+(headsmooth(1, 2) - headsmooth(2, 2)).^2);
for i = 2:length(head) - 1
    speed(i) = 0.5 * sqrt((headsmooth(i+1, 1) - headsmooth(i-1, 1)).^2+(headsmooth(i+1, 2) - headsmooth(i-1, 2)).^2);
end
speed(end) = sqrt((headsmooth(end, 1) - headsmooth(end-1, 1)).^2+(headsmooth(end, 2) - headsmooth(end-1, 2)).^2);

speed = speed * 30;
% 厘米每秒

speedsmooth = smoothdata(speed, 'gaussian', 10);
% totest if it is running at this time

%% time spent in each spatial bin
XbinEdge = -80:2:80;
YbinEdge = -80:2:80;
Tspatial = zeros(80, 80);
% 2*2 cm bin

% select the time scale
% trialFrame(1):trialFrame(end);
for xx = 1:80
    for yy = 1:80
        Tspatial(xx, yy) = length(find(headsmooth(trialFrame(1):trialFrame(end), 1) >= XbinEdge(xx) & headsmooth(trialFrame(1):trialFrame(end), 1) <= XbinEdge(xx+1) & headsmooth(trialFrame(1):trialFrame(end), 2) >= XbinEdge(yy) & headsmooth(trialFrame(1):trialFrame(end), 2) <= XbinEdge(yy+1) & speedsmooth(trialFrame(1):trialFrame(end)) > 10));

    end

end
% imshow(Tspatial,[],'Colormap',turbo,'Border','tight');
Tspatial = Tspatial ./ 30;
% 单位是秒

%% Load Spikes


allTS = feval(linermdl, 1:length(speedsmooth));
% load('D:\TrodesData\20230325\2-1test2\spikes.mat');
load('spikes.mat');
neuronNum = length(SpikeTimeStamps);

%% Rate maps

%trialFrame(1):trialFrame(end)
drawratemap = 0;





if drawratemap
    lxin = 23*cos(0:pi/3:2*pi);
    lyin = 23*sin(0:pi/3:2*pi);
    lxout = 39.5*cos(0:pi/3:2*pi);
    lyout = 39.5*sin(0:pi/3:2*pi);
    lxin = lxin +40.5;
    lyin = lyin +40.5;
    lxout = lxout +40.5;
    lyout = lyout +40.5;

    for neuronN = 1:neuronNum

        Thisneuron = SpikeTimeStamps{neuronN};

        Thisneuron = Thisneuron(Thisneuron > trialTS(1) & Thisneuron < trialTS(end));


        %第二三列是位置，第四列是速度 插值
        Thisneuron(:, 2) = interp1(allTS, headsmooth(:, 1), Thisneuron(:, 1));
        Thisneuron(:, 3) = interp1(allTS, headsmooth(:, 2), Thisneuron(:, 1));
        Thisneuron(:, 4) = interp1(allTS, speedsmooth, Thisneuron(:, 1));


        Fspatial = zeros(80, 80);

        for xx = 1:80
            for yy = 1:80
                Fspatial(xx, yy) = length(find(Thisneuron(:, 2) >= XbinEdge(xx) & Thisneuron(:, 2) <= XbinEdge(xx+1) & Thisneuron(:, 3) >= XbinEdge(yy) & Thisneuron(:, 3) <= XbinEdge(yy+1) & Thisneuron(:, 4) > 10));

            end

        end
        % 先把没经过的地方置零再卷积，但是画图的时候要把它画成NaN

        %
        Ratemapraw = Fspatial ./ Tspatial;
        Ratemap = Ratemapraw;
        Ratemap(isnan(Ratemapraw)) = 0;
        Ratemap(isinf(Ratemapraw)) = 0;
        gaussian_kernel = fspecial('gaussian', 5, 1);
        SmoothedRM = conv2(Ratemap, gaussian_kernel, "same");
        SmoothedRMdraw = SmoothedRM;
        SmoothedRMdraw(isnan(Ratemapraw)) = NaN;
        SmoothedRMdraw(isinf(Ratemapraw)) = NaN;

        % imshow(SmoothedRM,[],'Colormap',turbo,'Border','tight');
        figure;
        imagesc(SmoothedRM');
        ax = gca;
        ax.YDir = 'normal';
        colormap("turbo");
        colorbar;
        axis equal
        xlim([1, 80]);
        ylim([1, 80]);
        axis off
        line(lxin,lyin,'color','w','LineWidth',2);
        line(lxout,lyout,'color','w','LineWidth',2);


        filename = sprintf('place_cell_%d.png', neuronN);
        exportgraphics(gcf, filename);

    end
    % figure;
    % plot(Thisneuron(Thisneuron(:, 4) > 5, 2), Thisneuron(Thisneuron(:, 4) > 5, 3), 'k.');
    % axis equal;xlim([-80 80]); ylim([-80 80]);
end
close all

%只要这段时间

%% collect data for SVM

timewindow = 15000;
% 500ms
step = 1500;
% 50ms
RawData = zeros(length(trialTS(1):step:trialTS(end)), neuronNum+7);
%前六列分别是TS X Y theta Xn Yn Speed, 这里的n是normalized
%剩下右边的是神经元在这个时间窗的放电次数。


RawData(:, 1) = trialTS(1):step:trialTS(end);
RawData(:, 1) = RawData(:, 1) + timewindow / 2;
RawData(:, 2) = interp1(allTS, headsmooth(:, 1), RawData(:, 1));
RawData(:, 3) = interp1(allTS, headsmooth(:, 2), RawData(:, 1));
RawData(:, 4) = atan2(RawData(:, 3), RawData(:, 2));
RawData(:, 5) = cos(RawData(:, 4));
RawData(:, 6) = sin(RawData(:, 4));
RawData(:, 7) = interp1(allTS, speedsmooth, RawData(:, 1));

j = 1;

for i = trialTS(1):step:trialTS(end)

    for k = 1:neuronNum
        RawData(j, 7+k) = length(find(SpikeTimeStamps{k} > RawData(j, 1) & SpikeTimeStamps{k} <= RawData(j, 1)+timewindow));
    end
    j = j + 1;
end

ssfr = smoothdata(RawData(:, 8:end),'gaussian',10);

%% train SVM
% SVMXn = fitrsvm(RawData(RawData(:,7)>10, 8:end), RawData(RawData(:,7)>10, 5));
% disp('ok')
% 
% SVMYn = fitrsvm(RawData(RawData(:,7)>10, 8:end), RawData(RawData(:,7)>10, 6));

% SVMXn = fitrsvm(RawData(RawData(:,7)>10, 8:end), RawData(RawData(:,7)>10, 2));
% disp('ok')
% 
% SVMYn = fitrsvm(RawData(RawData(:,7)>10, 8:end), RawData(RawData(:,7)>10, 3));



%%
Xp = predict(SVMXn, RawData(:, 8:end));

Yp = predict(SVMYn, RawData(:, 8:end));

thetap = atan2(Yp, Xp);
theta = RawData(:, 4);

figure;
hold on

axis off
% set(gca, 'color', 'k')


v1 = VideoWriter('D:\TrodesData\20230325\2-1test1\demo1.mp4', 'MPEG-4');
v1.Quality = 15;
v1.FrameRate = 30;
open(v1);

for jj = fix(j/3):fix(j/2)
    %     if (RawData(i, 12) == 1)
    %         plot(cos(thetap(i)),sin(thetap(i)),'ro','LineWidth',2);
    %         plot(cos(theta(i)),sin(theta(i)),'r+','LineWidth',2);
    % %         plot(Xp(i),Yp(i),'ro','LineWidth',2);
    % %         plot(Xn(i),Yn(i),'r+','LineWidth',2);
    %     else
    % %         plot(Xp(i),Yp(i),'wo','LineWidth',2);
    % %         plot(Xn,Yn,'w+','LineWidth',2);
    %         plot(cos(thetap(i)),sin(thetap(i)),'wo','LineWidth',2);
    %         plot(cos(theta(i)),sin(theta(i)),'w+','LineWidth',2);
    %     end

    %     bar(0,cc(i),0.1,'blue');

    if RawData(jj,7)>10

        plot(cos(thetap(jj)), sin(thetap(jj)), 'go', 'LineWidth', 2, "MarkerSize", 10);

    end

    plot(cos(theta(jj)), sin(theta(jj)), 'r+', 'LineWidth', 2, "MarkerSize", 10);
    %     axis off


    set(gcf, 'color', 'k');
    xlim([-2, 2])
    ylim([-2, 2])


    drawnow;
    whatineed = getframe;
    writeVideo(v1, whatineed);
    cla;
end
close(v1);
close all
