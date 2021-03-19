%***************************************
%Author: Chaoqun Wang
%Date: 2019-10-15
%***************************************
%% 流程初始化
clc
clear all; close all;
x_I=10; y_I=700;           % 设置初始点
x_G=700; y_G=700;       % 设置目标点
Thr=50;                 % 设置目标点阈值
Delta= 30;              % 设置扩展步长
%% 建树初始化
T.v(1).x = x_I;         % T是我们要做的树，v是节点，这里先把起始点加入到T里面来
T.v(1).y = y_I; 
T.v(1).xPrev = x_I;     % 起始节点的父节点仍然是其本身
T.v(1).yPrev = y_I;
T.v(1).dist=0;          % 从父节点到该节点的距离，这里可取欧氏距离
T.v(1).indPrev = 0;     %父节点的索引
%% 开始构建树——作业部分
figure(1);
ImpRgb=imread('newmap.png');
Imp=rgb2gray(ImpRgb);
imshow(Imp)
xL=size(Imp,1);%地图x轴长度
yL=size(Imp,2);%地图y轴长度
hold on
plot(x_I, y_I, 'ro', 'MarkerSize',10, 'MarkerFaceColor','r');
plot(x_G, y_G, 'go', 'MarkerSize',10, 'MarkerFaceColor','g');% 绘制起点和目标点
count=1;  %计算点的数量，一开始只有根节点
for iter = 1:3000
    %Step 1: 在地图中随机采样一个点x_rand
    %提示：用（x_rand(1),x_rand(2)）表示环境中采样点的坐标
    x_rand = rand(1,2) * 800;  %在0到1之间生成随机数，生成一个1行2列的矩阵
    x_rand(1) = round(x_rand(1));
    x_rand(2) = round(x_rand(2));
    plot(x_rand(1), x_rand(2),'ro');
    %pause(0.1); 
    %Step 2: 遍历树，从树中找到最近邻近点x_near 
    %提示：x_near已经在树T里
    di = [];
    for i = 1:count  %count是当前建立的所有节点数
        dx = abs(x_rand(1) - T.v(i).x);
        dy = abs(x_rand(2) - T.v(i).y);
        d = (dx^2 + dy^2)^0.5;    %计算所有节点到新节点的距离
        di = [di d];   %记录下所有节点到新节点的距离
    end
    minimum = min(di);
    [x,y]=find(di==minimum);%find返回符合条件的行数列数，di的列数就是near节点加入的顺序(index)
    
    x_near=[T.v(y(1)).x T.v(y(1)).y];  %用y(1)就是节点在数中的索引  
    plot(x_near(1), x_near(2),'bx')
    previous_index = y(1);

    
    %Step 3: 扩展得到x_new节点
    %提示：注意使用扩展步长Delta
    dx = abs(x_rand(1) - x_near(1));
    dy = abs(x_rand(2) - x_near(2));
    d = (dx^2 + dy^2)^0.5;
    dx = x_rand(1) - x_near(1);
    dy = x_rand(2) - x_near(2);
    %Delta是步长
    dx = round(dx * Delta/d);
    dy = round(dy * Delta/d);
    x_new=[(x_near(1) + dx) (x_near(2) + dy)];
    %pause(2); 
    %检查节点是否是collision-free
    if ~collisionChecking(x_near,x_new,Imp) 
       continue;
    end
    
    %Step 4: 将x_new插入树T 
    %提示：新节点x_new的父节点是x_near
    count=count+1;
    T.v(count).x = x_new(1);         
    T.v(count).y = x_new(2); 
    T.v(count).xPrev = x_near(1);    
    T.v(count).yPrev = x_near(2);
    T.v(count).dist=d;   
    T.v(count).indPrev = previous_index;
    
    %Step 5:检查是否到达目标点附近 
    %提示：注意使用目标点阈值Thr，若当前节点和终点的欧式距离小于Thr，则跳出当前for循环
    dx = abs(x_new(1) - x_G);
    dy = abs(x_new(2) - y_G);
    d = (dx^2 + dy^2)^0.5;
    plot([x_near(1), x_new(1)],[x_near(2), x_new(2)], '--k');
    hold on;
    if d <= Thr
   %Step 6:将x_near和x_new之间的路径画出来
   %提示 1：使用plot绘制，因为要多次在同一张图上绘制线段，所以每次使用plot后需要接上hold on命令
   %提示 2：在判断终点条件弹出for循环前，记得把x_near和x_new之间的路径画出来
%         plot([x_G, x_new(1)],[y_G, x_new(2)]);
%         hold on;
        break;
   end
   %pause(0.2); %暂停0.1s，使得RRT扩展过程容易观察
end
%% 路径已经找到，反向查询
if iter < 2000
    path.pos(1).x = x_G; path.pos(1).y = y_G;
    path.pos(2).x = T.v(end).x; path.pos(2).y = T.v(end).y;
    pathIndex = T.v(end).indPrev; % 终点加入路径
    j=0;
    while 1
        path.pos(j+3).x = T.v(pathIndex).x;
        path.pos(j+3).y = T.v(pathIndex).y;
        pathIndex = T.v(pathIndex).indPrev;
        if pathIndex == 1
            break
        end
        j=j+1;
    end  % 沿终点回溯到起点
    path.pos(end+1).x = x_I; path.pos(end).y = y_I; % 起点加入路径
    for j = 2:length(path.pos)
        plot([path.pos(j).x; path.pos(j-1).x;], [path.pos(j).y; path.pos(j-1).y], 'b', 'Linewidth', 3);
    end
else
    disp('Error, no path found!');
end


