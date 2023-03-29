%% 贝叶斯解码
function [p_x_n] = BayesianDecoder(spkraster,placecell,timewindow,timestep,sampFreq)
%   输入：
%   spkraster: 待解码：细胞数目*时间
%   placecell: 已有信息：细胞数目*位置 (36个位置）
%   timewindow: 0.05s
%   timetimestep:   0.01s
%   sampFreq:采样频率 0.05(一秒采样20,000次)
%   输出：
%   p_x_n：位置*时间

if mod(sampFreq,2) == 0 
    nsample = timewindow*sampFreq+1;
else
    nsample = timewindow*sampFreq;
end
%   记录时间长度
%   1.将数据转换成raster（采样率是20,000）
%   2.数据预处理：将与视频帧不对齐的头和尾去掉
%   3.将静止的帧数 剔除(velocity < 5 cm/s for > 3 s）

%   初始化；
p_x_n = NaN(size(placecell,2),ceil(size(spkraster,2)/timestep/sampFreq));
ind2pxn = 0;

%   计算后验概率：
for pp = floor(nsample/2)+1:timestep*sampFreq:size(spkraster,2)-floor(nsample/2)%构建时间窗中心点
    ind2pxn = ind2pxn+1;%建立时间窗索引
    range = pp-floor(nsample/2):pp+floor(nsample/2)-1;%时间窗区间
    % 如果至少存在一个细胞发放
    if sum(sum(spkraster(:,range),2)~=0) >= 1 
        n = repmat(sum(spkraster(:,range),2),[1 size(placecell,2)]);
        p_x_n(:,ind2pxn) = (prod(placecell.^n,1))'.*exp(-timewindow*sum(placecell,1))'; %后验概率
    end    
    p_x_n(:,ind2pxn) = p_x_n(:,ind2pxn)/nansum(p_x_n(:,ind2pxn));        
end
% 结束条件
if size(spkraster,2) > max(range)
    if sum(sum(spkraster(:,max(range)+1:end))) > 0 
        ind2pxn = ind2pxn+1;
        n = repmat(sum(spkraster(:,max(range)+1:end),2),[1 size(placecell,2)]);
        p_x_n(:,ind2pxn) = (prod(placecell.^n,1))'.*exp(-timewindow*sum(placecell,1))';        
    end
    p_x_n(:,ind2pxn) = p_x_n(:,ind2pxn)/nansum(p_x_n(:,ind2pxn));   
    if ind2pxn<size(p_x_n,2)
        p_x_n(:,ind2pxn+1:end) = repmat(p_x_n(:,ind2pxn),1,size(p_x_n,2)-ind2pxn);
    end   
end