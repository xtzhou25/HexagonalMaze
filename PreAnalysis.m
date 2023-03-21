%调用trialseq函数

clear

%% 读取文件
tic
path = 'D:\TrodesData\230308\'; %！！【文件所在路径】
cd(path) %！！【文件所在路径】cd F:\20220710
% 读取视频文件
% 注意：视频从第0帧开始，第一列是帧数，第2、3列是Head,第5、6列是Neck
VideoFiles = dir('*.csv'); %返回一个struct矩阵，包含当前路径下的所有指定文件类型的文件名、大小等信息，详见files变量
Fileslen = size(VideoFiles, 1); %所分析session数目
% 读取Trial划分文件
MatFiles = dir('*.mat'); %读取划分trail的.mat文件
DinPath = [path,'din\'];
cd(DinPath)
DinFiles = dir('*.mat');
% %读取电生理时间戳


%%

%Fileslen

for o = 1:Fileslen

    csvfile = VideoFiles(o).name; %文件名（循环内依次操作）
    VfilePath = [path, csvfile]; %路径
    csvraw = readmatrix(VfilePath); %读取当前所指.csv文件
    Spikefile = MatFiles(o).name;
    Trialfile = MatFiles(o+1).name;
    Dinfile = DinFiles(o).name;
    SpikesfilePath = [path, Spikefile];
    TDfilePath = [path, Trialfile];
    DinfilePath = [DinPath, Dinfile];
    load(TDfilePath); % 读取TrialDivision的结果；变量名为trialnum
    load(SpikesfilePath); % 读取TrialDivision的结果；变量名为neuron
    load(DinfilePath)
    %合并变量

    %先用第一个neuron的时间戳对齐

    %% 对齐时间戳
    LEDChange = [trialnum(:, 1); trialnum(end, 1)] - [trialnum(1, 1); trialnum(:, 1)]; %标记trial分界处
    LEDtimestamp = find(LEDChange); %trial的第一帧数

    Spiketimestamp = din.in1(1).data;
    if length(din.in1(1).data) == 45
        x = double(Spiketimestamp(3:end-1)); %电生理时间戳,是第二到倒一
    elseif length(din.in1(1).data) == 44
        x = double(Spiketimestamp(2:end-1));
    end
    y = LEDtimestamp; %视频时间戳 tv1
    plot(x, y, 'rp');
    %用最小二乘法拟合,
    Lxx = sum((x - mean(x)).^2);
    Lxy = sum((x - mean(x)).*(y - mean(y)));
    b1 = Lxy / Lxx;
    b0 = mean(y) - b1 .* mean(x);
    Y = b1 .* x + b0; %拟合曲线方程
    hold on
    plot(x, Y); %将电生理 换算成视频帧！
    a = 30000 * b1;
    % 试验一次的结果：
    % b1 =9.987091369092669e-04;
    % b0 =-4.324489410061124e+02;
    %继而每一个电生理发放时间都可通过拟合公式转换成视频帧数
    %事实：电生理数据会长于视频【所以，我们需要去头】

    %%
    %第1列为：视频帧数：
    RawData = zeros(length(csvraw), 6+length(neuron));
    RawData(:, 1) = csvraw(:, 1) + 1;
    %第2列为：trial划分
    RawData(:, 2) = trialnum;
    %第3-6列为：头部、颈部坐标
    RawData(:, 3:4) = csvraw(:, 2:3);
    RawData(:, 5:6) = csvraw(:, 5:6);
    %第7列是原始速度（前后两点距离的均值）
    %速度转换系数：*1.533 换算成cm/s 1cm=5pixels
    rv = a / 5;
    RawData(1, 7) = sqrt((RawData(1, 3) - RawData(2, 3)).^2+(RawData(1, 4) - RawData(2, 4)).^2) * rv;
    for i = 2:length(RawData) - 1
        RawData(i, 7) = 0.5 * sqrt((RawData(i+1, 3) - RawData(i-1, 3)).^2+(RawData(i+1, 4) - RawData(i-1, 4)).^2) * rv;
    end
    RawData(end, 7) = sqrt((RawData(end, 3) - RawData(end-1, 3)).^2+(RawData(end, 4) - RawData(end, 4)).^2) * rv;
    %第8列是原始速度阈值判断 (velocity < 5 cm/s for > 3 s) excluded.
    m = find(RawData(:, 7) < 5); %找速度小于5cm/s
    mm = -[m(1); m] + [m; m(end)];
    mm2 = find(mm > 1);
    mm3 = -[1; mm2] + [mm2; mm2(end)];
    mm4 = find(mm3 > 89); %找到持续时间长于3s的序列
    RawData(:, 8) = ones(length(RawData), 1);
    for i = 1:length(mm4)
        RawData(m(mm2(mm4(i))-mm3(mm4(i))):m(mm2(mm4(i))-1), 8) = 0;
    end
    %第9、10列是平滑后头部坐标
    RawData(:, 9:10) = smoothdata(RawData(:, 3:4), 'gaussian', 10);
    %第11列是平滑后速度（前后两点距离的均值）
    RawData(1, 11) = sqrt((RawData(1, 9) - RawData(2, 9)).^2+(RawData(1, 4) - RawData(2, 4)).^2) * rv;
    for i = 2:length(RawData) - 1
        RawData(i, 11) = 0.5 * sqrt((RawData(i+1, 9) - RawData(i-1, 9)).^2+(RawData(i+1, 4) - RawData(i-1, 4)).^2) * rv;
    end
    RawData(end, 11) = sqrt((RawData(end, 9) - RawData(end-1, 9)).^2+(RawData(end, 4) - RawData(end, 4)).^2) * rv;
    %第12列是平滑后速度阈值判断 (velocity < 5 cm/s for > 3 s) excluded.
    n = find(RawData(:, 11) < 5); %找速度小于5cm/s
    nn = -[n(1); n] + [n; n(end)];
    nn2 = find(nn > 1);
    nn3 = -[1; nn2] + [nn2; nn2(end)];
    nn4 = find(nn3 > 89); %找到持续时间长于3s的序列
    RawData(:, 12) = ones(length(RawData), 1);
    for i = 1:length(nn4)
        RawData(n(nn2(nn4(i))-nn3(nn4(i))):n(nn2(nn4(i))-1), 12) = 0;
    end
    %第12列以后是每一帧对应的每个neuron（对应列）的发放次数

    %% Neurons是一个元胞组，每一个元胞代表一个细胞（包含其放电的时间戳;对应视频帧数；对应trial数）
    Neurons = neuron;
    for k = 1:size(Neurons, 1)
        Neurons{k, 1}(:, 1) = Neurons{k, 1}(:, 1) * 30000;
        Neurons{k, 1}(:, 2) = b1 .* Neurons{k, 1}(:, 1) + b0;
        Neurons{k, 1}(:, 3) = round(Neurons{k, 1}(:, 2));
        for j = 1:length(y) - 1
            for i = 1:length(Neurons{k, 1})
                if Neurons{k, 1}(i, 3) >= 1 && Neurons{k, 1}(i, 3) < y(1)
                    Neurons{k, 1}(i, 4) = 1;
                end
                if Neurons{k, 1}(i, 3) >= y(j) && Neurons{k, 1}(i, 3) < y(j+1)
                    Neurons{k, 1}(i, 4) = j + 1;
                end
                if Neurons{k, 1}(i, 3) >= y(end)
                    Neurons{k, 1}(i, 4) = length(y) + 1;
                end
            end
        end
    end

    %% 统计单帧的spike发放总数，每个细胞单独计算
    %yi1：第i个细胞，有spike发放的对应帧数；
    %yi2：第i个细胞，有spike发放的对应帧数的发放次数；
    %     [y11,y12]=trailseq(spiketimestamp{1,1}(spiketimestamp{1,1}(:,3)>0,3));
    %     y11=y11';
    %     y12=y12';

    % 节取Trial:1:42之间的神经元发放
    NeuronsV = cell(size(neuron, 1), 1);
    for p = 1:size(Neurons, 1)
        %         NeuronsV{p,1}=zeros(size(Neurons{p,1},1),2);
        [NeuronsV{p, 1}(:, 1), NeuronsV{p, 1}(:, 2)] = trailseq(Neurons{p, 1}(Neurons{p, 1}(:, 3) > 0, 3));

        %[NeuronsV{p,1}(:,1),NeuronsV{p,1}(:,2)]=trailseq(Neurons{p,1}(Neurons{p,1}(:,3)>0&Neurons{p,1}(:,4)<43,3));
        %         eval(['[y',num2str(p),'1,y',num2str(p),'2]=trailseq(spiketimestamp{',num2str(p),',1}(spiketimestamp{',num2str(p),',1}(:,3)>0,3));'])
    end

    %% 第7列以后是每一帧对应的每个neuron（对应列）的发放次数

    for j = 1:size(Neurons, 1)
        for ip = 1:length(RawData)
            for il = 1:length(NeuronsV{j, 1})

                if ip == NeuronsV{j, 1}(il, 1)
                    RawData(ip, j+12) = NeuronsV{j, 1}(il, 2); %第12列以后
                end
            end
        end
    end

    save([path,'Organized\', VideoFiles(o).name(1:end-4), '_RawData.mat'], 'RawData');
    save([path,'Organized\', VideoFiles(o).name(1:end-4), '_RegParam.mat'], 'b0', 'b1'); %Regression parameters

end
toc

%%
