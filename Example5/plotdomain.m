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
colormap('white');axis image;axis off; ax=axis;axis(ax*1.001);view([0,90]);
%text(-0.05,0.5,'p=1','HorizontalAlignment','center','Color','k','FontSize',14)
%text(1.05,0.5,'p=0','HorizontalAlignment','center','Color','k','FontSize',14)
%text(0.5,-0.04,'q=0','HorizontalAlignment','center','Color','k','FontSize',14)
%text(0.5,1.04,'q=0','HorizontalAlignment','center','Color','k','FontSize',14)
% ParametersofFracture: x_c,y_c,length,theta,width,permeability_nu,permeability_sigma,xa,ya,xb,yb,# of elements passed by this fracture. 
NumberofFractures = 6; 
ParametersofFracture = zeros(NumberofFractures, 11); 
ParametersofFracture(:,1:4) = [ 0.5, 0.5, 1, 0;...%... % position of fractures
                                0.5, 0.5, 1, pi/2;
                                0.75, 0.75, 0.5, 0;
                                0.75, 0.75, 0.5, pi/2;
                                0.625, 0.625, 0.25, 0;
                                0.625, 0.625, 0.25, pi/2];
% plot the triangulation and fractures
ParametersofFracture(:,7) = ParametersofFracture(:,1) - 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,8) = ParametersofFracture(:,2) - 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
ParametersofFracture(:,9) = ParametersofFracture(:,1) + 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,10) = ParametersofFracture(:,2) + 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
figure;surf(reshape(coordinates(:,1),Cell_N+1,Cell_M+1),reshape(coordinates(:,2),Cell_N+1,Cell_M+1),zeros(Cell_N+1,Cell_M+1));
colormap('white');axis image;axis off; ax=axis;axis(ax*1.001); view([0,90]); hold on;
plot(ParametersofFracture(:,[7,9])',ParametersofFracture(:,[8,10])',...
   'r-','LineWidth',2); hold off; 
text(-0.1,0.5,'q_N=-1','HorizontalAlignment','center','Color','k','FontSize',16)
text(1.1,0.5,'p_D=1','HorizontalAlignment','center','Color','k','FontSize',16)
text(0.5,-0.045,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',16)
text(0.5,1.05,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',16)
