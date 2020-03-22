%%  MeshGenerator��һ���򵥵�2D��������������
%   Coded by  : ����Ң
%   Date      : ������2015-8-17
%               ���һ���޸���2015-8-17
% �����飺
% ��Ҫ�����ĸ�txt�ļ���MG_nodes.txt��MG_elements.txt,MG_boundary_nodes.txt��MG_boundary_edges.txt��
% �ֱ����ڴ�Žڵ�ļ������꣬����Ԫ���������ţ��߽���߼��ж�����ͱ߽�߶˵��š�
% MG_nodes.txt������Ϊ�ڵ�����������Ϊ2����i�����зֱ��Ӧ��i���ڵ��x��y���ꣻ
% MG_elements.txt������Ϊ����Ԫ����������Ϊ3����i�����зֱ��Ӧ�ڵ�ֲ����Ϊ1��2��3��
% �����ڵ�������ţ����ڵ㰴��ʱ�����У�
% MG_boundary_nodes.txt������Ϊ�ڵ���������i��Ϊ0��ʾ����һ���ڵ㣬Ϊ1��ʾ����һ���߽�㣻
% MG_boundary_edges.txt������Ϊ�߽���Ԫ����������Ϊ2����i�����б�ʾ��i�ߵ����˵��š�
%%  [1]���ɽڵ�����p������Ԫ����t��
clc,clear
tic
% �����ʷ�ʾ����
%--------------------------------------------------------------------------
% Example: (Rectangle with circular hole, refined at circle boundary)
%    fd=@(p) drectangle(p,0,1,0,1);
%    [p,t]=distmesh2d(fd,@huniform,0.05,[0,0;1,1],[0,0;1,0;1,1;0,1]);
% Example: (Rectangle with circular hole, refined at circle boundary)
    fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
    fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
    [p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
%--------------------------------------------------------------------------
% �����Ż�(��Ȼһ�������ûʲô����)Remove duplicated/unused nodes and fix element orientation.
[p,t]=fixmesh(p,t);
% ���������
bend1=max(t,[],2)-min(t,[],2);
ave_bend1=mean(bend1);
% �ڵ����˳����ӻ���
hold on
pic.h=plot(p(1,1),p(1,2),'ro');
for ii=1:size(p,1)
    set(pic.h,'XData',p(1:ii,1), 'YData',p(1:ii,2));
    drawnow 
end
set(pic.h,'XData',[], 'YData',[]); % ����
drawnow 
% ��Ԫ����˳����ӻ���
hold on
[pc,r]=circumcenter(p,t);
pic.h=plot(pc(1,1),pc(1,2),'r.','markersize',10);
for ii=1:size(t,1)
    set(pic.h,'XData',pc(1:ii,1), 'YData',pc(1:ii,2),'markersize',10);
    drawnow 
end
set(pic.h,'XData',[], 'YData',[]); % ����
drawnow 

%%  [2]���ݱ���
% �ո���Ϊ�����
% д��ڵ���Ϣ��
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
% д������Ԫ��Ϣ��
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
%%  ���̴�����ʱ���򣨸ò���������
% ��һ��һ�������������Ϊ
% һ������£�ԭʼ���ɵĽڵ����Ԫ���������Ѿ��㹻�ã�������������������������м��ж�������
% �Ż��к�Ч����Ӧ������ѡ���ʹ�ã�������ʱ������Ǽ���ǰ���Ѿ�����һ�����ˣ������޷�����һ��
% ��������ã���ʵ��ȫ���ؼ���һ�顣
triangulation_rcm('MG')  % ����������������������Ӧ�ļ���rcm��
triangulation_orient('MG') % ��ʱ��
p2=load('MG_rcm_nodes.txt');
t2=load('MG_rcm_elements.txt');
% ���������
bend2=max(t2,[],2)-min(t2,[],2);
ave_bend2=mean(bend2);
% �ڵ����˳����ӻ���
hold on
pic2.h=plot(p2(1,1),p2(1,2),'bo');
for ii=1:size(p2,1)
    set(pic2.h,'XData',p2(1:ii,1), 'YData',p2(1:ii,2));
    drawnow 
end
set(pic2.h,'XData',[], 'YData',[]);% ����
drawnow 
% ��Ԫ����˳����ӻ���
hold on
[pc2,r2]=circumcenter(p2,t2);
pic.h=plot(pc2(1,1),pc2(1,2),'b.','markersize',10);
for ii=1:size(t2,1)
    set(pic.h,'XData',pc2(1:ii,1), 'YData',pc2(1:ii,2),'markersize',10);
    drawnow 
end
set(pic.h,'XData',[], 'YData',[]); % ����
drawnow 

%%  Ѱ�߽�㲢����(����Լ��˺�Ĵ������̣����ٶ���һ�����������������������ݶ���ԭʼ��)
triangulation_boundary_nodes ( 'MG' )% Ѱ�Ҳ�����߽��
boundary_nodes=load('MG_boundary_nodes.txt');
% �߽����ӻ�
hold on
for ii=1:length(boundary_nodes) % �����boundary_nodesʵ�����Ǹ����ܽ��ȳ��ȵġ��߼�������
   if boundary_nodes(ii)==1 % ���ĳһλ��Ϊ1����˵������Ӧλ�õĽڵ�Ϊ�߽�㣬����Ϊ�ڵ㡣
       plot(p(ii,1),p(ii,2),'r.')
   end
end
%%  Ѱ�߽���Ԫ�����棨�����У���һ���Ǳ߽�Ԫ���ͣ�1-��߽磬2-�ڱ߽磬�ڶ�����Ϊ�˵��ţ�
triangulation_boundary_edges ( 'MG' )% Ѱ�Ҳ�����߽���Ԫ
boundary_edges=load('MG_boundary_edges.txt');
% �߽���ӻ�
hold on
for ii=1:size(boundary_edges,1)
   if boundary_edges(ii,1)==1% B�ĵ�һ�е���1��������߽�����
      plot([p(boundary_edges(ii,2),1);p(boundary_edges(ii,3),1)],...
      [p(boundary_edges(ii,2),2);p(boundary_edges(ii,3),2)],'ro-','LineWidth',2)
   end
   if boundary_edges(ii,1)==2% B�ĵ�һ�е���2����ʾ�ڱ߽�����
      plot([p(boundary_edges(ii,2),1);p(boundary_edges(ii,3),1)],...
      [p(boundary_edges(ii,2),2);p(boundary_edges(ii,3),2)],'b*-','LineWidth',2)
   end
   drawnow
end
toc
