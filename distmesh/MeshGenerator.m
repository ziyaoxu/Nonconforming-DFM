%%  MeshGenerator：一个简单的2D三角网格生成器
%   Coded by  : 徐梓尧
%   Date      : 创建于2015-8-17
%               最后一次修改于2015-8-17
% 输出简介：
% 主要生成四个txt文件：MG_nodes.txt，MG_elements.txt,MG_boundary_nodes.txt和MG_boundary_edges.txt；
% 分别用于存放节点的几何坐标，三角元的三顶点标号，边界点逻辑判断数组和边界边端点标号。
% MG_nodes.txt：行数为节点总数，列数为2；第i行两列分别对应第i个节点的x、y坐标；
% MG_elements.txt：行数为三角元总数，列数为3；第i行三列分别对应节点局部编号为1、2、3的
% 三个节点的总体编号，三节点按逆时针排列；
% MG_boundary_nodes.txt：行数为节点总数；第i行为0表示这是一个内点，为1表示这是一个边界点；
% MG_boundary_edges.txt：行数为边界线元总数，列数为2；第i行两列表示第i边的两端点标号。
%%  [1]生成节点数据p和三角元数据t：
clc,clear
tic
% 网格剖分示例：
%--------------------------------------------------------------------------
% Example: (Rectangle with circular hole, refined at circle boundary)
%    fd=@(p) drectangle(p,0,1,0,1);
%    [p,t]=distmesh2d(fd,@huniform,0.05,[0,0;1,1],[0,0;1,0;1,1;0,1]);
% Example: (Rectangle with circular hole, refined at circle boundary)
    fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
    fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
    [p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
%--------------------------------------------------------------------------
% 数据优化(虽然一般情况下没什么卵用)Remove duplicated/unused nodes and fix element orientation.
[p,t]=fixmesh(p,t);
% 带宽情况：
bend1=max(t,[],2)-min(t,[],2);
ave_bend1=mean(bend1);
% 节点编码顺序可视化：
hold on
pic.h=plot(p(1,1),p(1,2),'ro');
for ii=1:size(p,1)
    set(pic.h,'XData',p(1:ii,1), 'YData',p(1:ii,2));
    drawnow 
end
set(pic.h,'XData',[], 'YData',[]); % 擦除
drawnow 
% 单元编码顺序可视化：
hold on
[pc,r]=circumcenter(p,t);
pic.h=plot(pc(1,1),pc(1,2),'r.','markersize',10);
for ii=1:size(t,1)
    set(pic.h,'XData',pc(1:ii,1), 'YData',pc(1:ii,2),'markersize',10);
    drawnow 
end
set(pic.h,'XData',[], 'YData',[]); % 擦除
drawnow 

%%  [2]数据保存
% 空格作为间隔符
% 写入节点信息：
  output_unit = fopen ( 'MG_nodes.txt', 'wt' );
  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'Error!' );
  end
  for ii = 1 : size(p,1)
    for jj = 1 : size(p,2)
      fprintf ( output_unit, '  %14f', p(ii,jj) );
    end
    fprintf ( output_unit, '\n' );
  end
  fclose ( output_unit );
% 写入三角元信息：
  output_unit = fopen ( 'MG_elements.txt', 'wt' );
  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'Error!' );
  end
  for ii = 1 : size(t,1)
    for jj = 1 : size(t,2)
      fprintf ( output_unit, '  %12d', t(ii,jj) );
    end
    fprintf ( output_unit, '\n' );
  end
  fclose ( output_unit );
%dlmwrite('MG_nodes.txt', p, ' ');
%dlmwrite('MG_elements.txt', t, ' ');
%%  缩短带宽，逆时针序（该步可跳过）
% 这一步一般可以跳过，因为
% 一般情况下，原始生成的节点和面元数据质量已经足够好，经过处理带宽反而略增，除非中间有洞的区域
% 才会有好效果，应该有所选择的使用；此外逆时针序更是鸡肋前面已经调过一次序了，这里无非是起到一个
% 检验的作用，其实完全不必加这一块。
triangulation_rcm('MG')  % 调点序，缩带宽，并保存于相应文件（rcm）
triangulation_orient('MG') % 逆时针
p2=load('MG_rcm_nodes.txt');
t2=load('MG_rcm_elements.txt');
% 带宽情况：
bend2=max(t2,[],2)-min(t2,[],2);
ave_bend2=mean(bend2);
% 节点编码顺序可视化：
hold on
pic2.h=plot(p2(1,1),p2(1,2),'bo');
for ii=1:size(p2,1)
    set(pic2.h,'XData',p2(1:ii,1), 'YData',p2(1:ii,2));
    drawnow 
end
set(pic2.h,'XData',[], 'YData',[]);% 擦除
drawnow 
% 单元编码顺序可视化：
hold on
[pc2,r2]=circumcenter(p2,t2);
pic.h=plot(pc2(1,1),pc2(1,2),'b.','markersize',10);
for ii=1:size(t2,1)
    set(pic.h,'XData',pc2(1:ii,1), 'YData',pc2(1:ii,2),'markersize',10);
    drawnow 
end
set(pic.h,'XData',[], 'YData',[]); % 擦除
drawnow 

%%  寻边界点并保存(这里，以及此后的处理流程，都假定上一步缩带宽被跳过；所处理数据都是原始的)
triangulation_boundary_nodes ( 'MG' )% 寻找并保存边界点
boundary_nodes=load('MG_boundary_nodes.txt');
% 边界点可视化
hold on
for ii=1:length(boundary_nodes) % 这里的boundary_nodes实际上是个与总结点等长度的“逻辑”数组
   if boundary_nodes(ii)==1 % 如果某一位置为1，则说明该相应位置的节点为边界点，否则为内点。
       plot(p(ii,1),p(ii,2),'r.')
   end
end
%%  寻边界线元并保存（共三列，第一列是边界元类型：1-外边界，2-内边界，第二三列为端点标号）
triangulation_boundary_edges ( 'MG' )% 寻找并保存边界线元
boundary_edges=load('MG_boundary_edges.txt');
% 边界可视化
hold on
for ii=1:size(boundary_edges,1)
   if boundary_edges(ii,1)==1% B的第一列等于1，代表外边界类型
      plot([p(boundary_edges(ii,2),1);p(boundary_edges(ii,3),1)],...
      [p(boundary_edges(ii,2),2);p(boundary_edges(ii,3),2)],'ro-','LineWidth',2)
   end
   if boundary_edges(ii,1)==2% B的第一列等于2，表示内边界类型
      plot([p(boundary_edges(ii,2),1);p(boundary_edges(ii,3),1)],...
      [p(boundary_edges(ii,2),2);p(boundary_edges(ii,3),2)],'b*-','LineWidth',2)
   end
   drawnow
end
toc
