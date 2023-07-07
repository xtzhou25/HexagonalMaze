% This script read mp4 videos to get the brightness of two LED indicators,
% thus we can divide each video into several trials.

clear
[file, path] = uigetfile('E:\HexagonData\Videos\20230626\example.mp4', 'MultiSelect', 'on');

range = 2; %calculate the pixels around the positions.

% 就用它�?�，一帧一帧读


if isa(file, 'cell')
    nn = length(file); % 有多少个文件�?分�?
else
    nn = 1;
    tempfile{1} = file;
    file = tempfile;
end

%%

tic

for n = 1:nn
    %     tic
    v0 = VideoReader([path, file{n}]); %#ok<TNMLP>
    l = v0.NumFrames;
    %         l = 10000;
    if n == 1
        [redlight, greenlight] = getlightposition(v0);
    end

    rvalue = zeros(l, 1);
    gvalue = zeros(l, 1);


    for i = 1:l
        whatineed = read(v0, i);
        %     imshow(whatineed);
        %     hold on
        %     plot(redlight(1),redlight(2),'b.','Markersize',13);
        %     plot(greenlight(1),greenlight(2),'b.','Markersize',13);
        rhsv = rgb2hsv(whatineed(redlight(2)-range:redlight(2)+range, redlight(1)-range:redlight(1)+range, :));
        ghsv = rgb2hsv(whatineed(greenlight(2)-range:greenlight(2)+range, greenlight(1)-range:greenlight(1)+range, :));
        rvalue(i) = sum(rhsv(:, :, 3), 'all');
        gvalue(i) = sum(ghsv(:, :, 3), 'all');

    end
    %     figure;
    %     hold on
    %     plot(1:l,rvalue(1:l),'r.-',1:l,gvalue(1:l),'g.-');
    %     drawnow;

    [idxr, Cr] = kmeans(rvalue, 2); %用�?�类算出�?��?�亮还是�?�
    [idxg, Cg] = kmeans(gvalue, 2);
    %     rgvalue = [rvalue,gvalue];
    changethreshold = [min(Cr) + (max(Cr) - min(Cr)) / 3, min(Cg) + (max(Cg) - min(Cg)) / 3]; %三等分点作为判断阈值
    trialnum = zeros(l, 1);
    trialnum(1) = 1; % 1是red，2是green

    i = 2;
    while i < l+1
        if mod(trialnum(i-1), 2) %�?一帧还是红
            if gvalue(i) > changethreshold(2)
                trialnum(i:min(i+29, l)) = trialnum(i-1) + 1;
                i = min(i+29, l); %
            else
                trialnum(i) = trialnum(i-1);
                i = i + 1;
            end
        else
            if rvalue(i) > changethreshold(1)
                trialnum(i:min(i+29, l)) = trialnum(i-1) + 1;
                i = min(i+29, l); %
            else
                trialnum(i) = trialnum(i-1);
                i = i + 1;
            end
        end
    end
    %     yyaxis right
    %     plot(1:l,trialnum(1:l),'k.-');
    %     hold off
    save([path, file{n}(1:end - 3), 'mat'], 'trialnum');
    %     toc
end
toc


function [redlight, greenlight] = getlightposition(v0)
figure;
whatineed = read(v0, 20);
imshow(whatineed);
set(gcf, 'outerposition', get(0, 'screensize'));
hold on;
redlight = [0, 0]; % 红�?�
[redlight(1), redlight(2)] = ginput(1);
redlight = round(redlight);
plot(redlight(1), redlight(2), 'b.', 'Markersize', 13);
greenlight = [0, 0]; % 绿�?�
[greenlight(1), greenlight(2)] = ginput(1);
greenlight = round(greenlight);
plot(greenlight(1), greenlight(2), 'b.', 'Markersize', 13);
hold off
cla
close all
end
