clc,clear
% set the geometry and triangulation parameters:
xmin = 0 ; xmax = 1 ; % domain 
ymin = 0 ; ymax = 1 ; % domain
% x direction is divided into N parts ,y direction is divided into M parts:
Cell_N = 1 ; Cell_M = 1 ; 
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ;
[coordinates,elements4,EtoEmap] = RectangularMesh( xmin,xmax,ymin,ymax,Cell_N,Cell_M );% periodic boundary condition
% plot the triangulation and fractures
figure;surf(reshape(coordinates(:,1),Cell_N+1,Cell_M+1),reshape(coordinates(:,2),Cell_N+1,Cell_M+1),zeros(Cell_N+1,Cell_M+1));
colormap('white');axis image;axis off; ax=axis;axis(ax*1.001); view([0,90]); hold on;
% heart curve:
%N_partition = 10007; t_start = 0; t_end = 2*pi; 
%t_partition = linspace(t_start,t_end,N_partition)';
%x_partition = -(13*cos(t_partition)-5*cos(2*t_partition)-2*cos(3*t_partition)-cos(4*t_partition));
%y_partition = 16*sin(t_partition).^3; 
%x_partition = x_partition/50+1/2;
%y_partition = y_partition/50+1/2;
% circle
%N_partition = 10007; t_start = 0; t_end = 2*pi; 
%t_partition = linspace(t_start,t_end,N_partition)';
%Radius = 0.25; x_center = 0.5+(rand-0.5)*1e-4; y_center = 0.5+(rand-0.5)*1e-4;
%x_partition = Radius*cos(t_partition) + x_center; y_partition = Radius*sin(t_partition) + y_center;
% semicircles:
N_partition = 10007; t_start = 0; t_end = 2*pi; 
t_partition = linspace(t_start,t_end,N_partition)';
x_partition = 3/8-1/8*cos(t_partition); 
x_partition(t_partition>pi)=5/8+1/8*cos(t_partition(t_partition>pi));
y_partition = 1/2+1/8*sin(t_partition); 
% draw the parametric curve:
hold on; plot(x_partition,y_partition,'r-','LineWidth',2)

text(-0.08,0.5,'p_D=1','HorizontalAlignment','center','Color','k','FontSize',16)
text(1.08,0.5,'p_D=0','HorizontalAlignment','center','Color','k','FontSize',16)
text(0.5,-0.045,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',16)
text(0.5,1.045,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',16)