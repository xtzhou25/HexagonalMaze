% clear
load('D:\TrodesData\230307\Organized\train1-2-20230307test02_RawData.mat');
smoothfr = smoothdata(RawData(:, 13:end), 1, 'gaussian', 70);

MazeCenter = [362, 367];
pace = 0;
if pace == 1
    data = RawData((RawData(:, 12) == 1), :);
    ssfr = smoothfr((RawData(:, 12) == 1), :); %selected ssfr
else
    data = RawData;
    ssfr = smoothfr; %selected ssfr
end

X = data(:, 9) - MazeCenter(1);
Y = -data(:, 10) + MazeCenter(2);
theta = atan2(Y, X);
Xn = cos(theta);
Yn = sin(theta);

mymdl1 = fitrsvm(ssfr, Xn);
mymdl2 = fitrsvm(ssfr, Yn);
% mymdl3 = fitrsvm(ssfr,theta);

%%
Xp = predict(mymdl1, ssfr);
Yp = predict(mymdl2, ssfr);

% ssfr2 = smoothdata(RawData(:, 13:end), 1, 'gaussian', 30);
% Xp = predict(mymdl1, ssfr2);
% Yp = predict(mymdl2, ssfr2);

% thetap = predict(mymdl3,ssfr);
thetap = atan2(Yp, Xp);
confidence = sqrt(Xp.^2+Yp.^2);
cc=normalize(confidence,'range');

%% 可视化
% plot(X,Y);hold on;plot(Xp,Yp);

% figure;
% plot(Y);
% hold on;
% plot(Yp);
% hold off

% figure;plot(theta);hold on;plot(thetap);


figure;
hold on

axis off
% set(gca, 'color', 'k')


v1 = VideoWriter('D:\TrodesData\demo22.mp4','MPEG-4');
v1.Quality = 15;
v1.FrameRate = 30;
open(v1);


for i = 1:18390
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
    plot(cos(thetap(i)), sin(thetap(i)), 'go', 'LineWidth', 2,"MarkerSize",10);
    plot(cos(theta(i)), sin(theta(i)), 'r+', 'LineWidth', 2,"MarkerSize",10);
%     axis off

    
    set(gcf, 'color', 'k');
    xlim([-2, 2])
    ylim([-2, 2])


    drawnow;
        whatineed = getframe;
        writeVideo(v1,whatineed);
    cla;
end
close(v1);
%


% cr = length(find(label == (data(:,4)>367)))./size(data,1);
% scr = length(find(label == (data(randperm(size(data,1)),3)>362)))./size(data,1);

%% 定义loss
Xs = Xn(randperm(size(Xn,1)));
Ys = Yn(randperm(size(Yn,1)));
% totaltloss = 0;
% totalsloss = 0;


for k = 1:size(ssfr,1)
%     thetaloss(k) = acosd(dot([Xn(k) Yn(k)],[Xp(k),Yp(k)])/(norm([Xn(k) Yn(k)])*norm([Xp(k) Yp(k)])));
%     shuffleloss(k) = acosd(dot([Xs(k) Ys(k)],[Xp(k),Yp(k)])/(norm([Xs(k) Ys(k)])*norm([Xp(k) Yp(k)])));
    thetaloss(k) = dot([Xn(k) Yn(k)],[Xp(k),Yp(k)])/(norm([Xn(k) Yn(k)])*norm([Xp(k) Yp(k)]));
    shuffleloss(k) = dot([Xs(k) Ys(k)],[Xp(k),Yp(k)])/(norm([Xs(k) Ys(k)])*norm([Xp(k) Yp(k)]));
end

%%
j=1;
clear totaltloss
clear totalsloss
tw = 600;
for k = 1:tw:size(ssfr,1)-tw
    
    totaltloss(j) =  sum(thetaloss(k:k+tw))./tw;
    totalsloss(j) =  sum(shuffleloss(k:k+tw))./tw;


j=j+1;

end
plot(totaltloss,'r','LineWidth',2);
hold on
plot(totalsloss,'k','LineWidth',2);
legend('Decoder','Shuffle');

% j=0;
% for i = 1:90:size(ssfr,1)-90
%     j=j+1;
%     for k = i:i+80
%     thetaloss = acosd(dot([Xn(k) Yn(k)],[Xp(k),Yp(k)])/(norm([Xn(k) Yn(k)])*norm([Xp(k) Yp(k)])));
%     shuffleloss = acosd(dot([Xs(k) Ys(k)],[Xp(k),Yp(k)])/(norm([Xs(k) Ys(k)])*norm([Xp(k) Yp(k)])));
%     totaltloss(j) =  totaltloss(j) + thetaloss;
%     totalsloss(j) =  totalsloss(j) + shuffleloss;
%     end
%     totaltloss(j) = totaltloss(j)./90;
%     totalsloss(j) = totalsloss(j)./90;
% 
% end

%%
% totaltloss = totaltloss/size(ssfr,1);
% totalsloss = totalsloss/size(ssfr,1);
