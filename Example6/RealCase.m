function RealCase
% time : 2019.6.26 - 2019.6.26
% author : xuziyao
% this program is code for Continuous Finite Element Method for
% Poisson equation with delta diffusion coefficient on rectangular meshes.
% PDE model: 
% -div(D*grad(p)) = f, where D contains the delta fracture coefficient
% Neumann & Dirichlet boundary condition 
% upper and lower boundary: Neumann
% left and right boundary: Dirichlet
% Basis function is the bilinear Lagrange basis on the [0,1]*[0,1] reference square.
% phi_1(r,s) = (1-r)(1-s)
% phi_2(r,s) = r(1-s)
% phi_3(r,s) = rs
% phi_4(r,s) = (1-r)s
% Gradrphi_1 = -(1-s); Gradsphi_1 = -(1-r).
% Gradrphi_2 =   1-s ; Gradsphi_2 =    -r.
% Gradrphi_3 =    s  ; Gradsphi_3 =     r.
% Gradrphi_4 =   -s  ; Gradsphi_4 =   1-r .
% FINISHED
clc,clear
format long
% set the geometry and triangulation parameters:
xmin = 0 ; xmax = 700 ; % domain 
ymin = 0 ; ymax = 600 ; % domain
% x direction is divided into N parts ,y direction is divided into M parts:
Cell_N = 70*1.5 ; Cell_M = 60*1.5 ; 
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ;
[coordinates,elements4,~] = RectangularMesh( xmin,xmax,ymin,ymax,Cell_N,Cell_M );% periodic boundary condition
EtoEmap = Rectangulation_neighbor_rectangles( elements4 ); 
Tarea =  hx*hy; Elenth = [hx;hy;hx;hy];
Jacobimat = [1/hx,0;0,1/hy];%Jacobimat=[Dr/Dx,Dr/Dy ; Ds/Dx,Ds/Dy] 
boundary_edge = boundedges_Rec( elements4 ) ;
neumann = zeros(size(boundary_edge)) ; neumann_count = 0;  
% ParametersofFracture: x_c,y_c,length,theta,width,permeability,xa,ya,xb,yb,# of elements passed by this fracture. 
%FID,START_X,START_Y,END_X,END_Y
FractureData = ...
[1,269.611206,152.05243,356.9240112,310.14123;
2,249.5117187,514.990780001,272.218872,470.97082;
3,258.3590698,515.574580001,271.9851684,490.9682;
4,270.6622924,524.702640001,269.1347046,147.78143;
5,355.8302002,348.479800001,337.5810733205,600;
6,366.9730835,338.132990001,426.9185141723,600;
7,198.237915,222.724420001,175.1561889,597.603030001;
8,151.2785034,261.724610001,154.4623059774,600;
9,29.5026855,300.724610001,96.3599853,514.82739;
10,386.0808105,33.3621800002,440.585083,275.191830001;
11,459.6350708,40.2413900001,461.751709,204.812620001;
12,297.180603,237.62103,468.1018066,40.2413900001;
13,312.5264892,272.01678,417.3016967,140.7832;
14,330.5181884,298.47522,439.5266723,156.6582;
15,340.5723877,320.70019,367.5598755,286.304380001;
16,492.9725952,312.762820001,576.5811157,419.6546;
17,505.6726684,309.05859,576.0520019,405.367190001;
18,537.4227905,297.94598,623.3187866,376.68463;
19,322.5338745,380.76941,521.8778076,593.552180001;
20,344.9320678,481.56122,409.8867798,503.959410001;
21,371.8098755,468.12219,510.6787109,383.009210001;
22,432.2849731,510.678830001,642.8280029,374.04999;
23,527.528634971,600,700,473.015615092;
24,0,333.73321,441.2443847,0;
25,13.4389038,342.692380001,347.171875,595.791990001;
26,22.3981933,450.203790001,311.3347778,291.176630001;
27,26.8778076,506.199220001,199.343811,400.92779;
28,44.7963867,528.597410001,365.0905151,342.692380001;
29,378.5294189,309.095210001,512.918518,116.470640001;
30,461.4027099,253.099610001,530.8370971,134.38922;
31,347.171875,374.04999,640.5881958,253.099610001;
32,490.5203857,268.77844,564.4343872,145.58844;
33,47.0361938,181.425410001,53.7556152,306.85541;
34,382.4152832,424.151000001,447.8997192,371.76343;
35,587.9967651,394.78222,549.1029663,362.635190001;
36,589.9812011,393.59161,527.6716919,313.8194;
37,597.125,378.90722,533.6248169,295.960200001;
38,533.6248169,448.75738,453.8527832,326.91638;
39,511.7966919,461.85419,489.5715942,395.17901;
40,565.3748779,425.34161,483.6184692,315.40698;
41,534.4185791,407.482240001,467.3466186,315.803830001;
42,627.2874756,527.3388,574.8999023,498.763610001;
43,644.3532104,519.00439,586.4093017,490.03241;
44,655.8626098,502.335630001,602.6812133,476.53863;
45,415.355896,585.679380001,391.9401855,561.47003;
46,417.3402099,578.535580001,397.8933105,554.326230001;
47,403.0526733,592.029420001,382.0183105,561.86682;
48,495.1278686,505.113580001,468.1403198,481.30121;
49,533.6248169,254.84381,420.9121093,159.196590001;
50,508.6217041,221.10943,441.152771,159.59363;
51,418.5308838,229.04681,312.961914,93.3154300004;
52,362.5714111,174.6748,322.883789,120.69983;
53,357.8088989,216.3468,295.102478,114.74658;
54,402.2589111,283.41882,366.1433105,226.66559;
55,337.5681762,253.256220001,374.4776001,211.18744;
56,386.7808838,264.765620001,509.8123169,101.25281;
57,473.2996826,278.65643,561.0092163,144.909240001;
58,471.7122192,253.653200001,554.6593017,129.034240001;
59,559.0249023,219.125,567.3593139,153.64044;
60,567.7561035,214.759400001,573.7092895,162.37182;
61,574.8999023,215.553040001,579.6624145,173.88104;
62,557.0404663,285.006410001,600.6968994,325.48761;
63,565.3748779,283.022030001,607.0468139,323.503230001];
FractureData(:,1)=[]; FractureData = (FractureData-300)*(1-1e-10)+300;
NumberofFractures = size(FractureData,1); 
ParametersofFracture = zeros(NumberofFractures, 11); 
ParametersofFracture(:,7:10) = FractureData;
ParametersofFracture(:,1) =  (FractureData(:,1)+FractureData(:,3))/2;% position of fractures
ParametersofFracture(:,2) =  (FractureData(:,2)+FractureData(:,4))/2;% position of fractures
ParametersofFracture(:,3) =  sqrt((FractureData(:,1)-FractureData(:,3)).^2+(FractureData(:,2)-FractureData(:,4)).^2);% position of fractures
ParametersofFracture(:,4) =  atan2(FractureData(:,4)-FractureData(:,2),FractureData(:,3)-FractureData(:,1));% position of fractures
ParametersofFracture(:,5) = 1e-2; % thickness of fractures
ParametersofFracture(:,6)= 1e-8 ; % permeability of fractures
% plot the triangulation and fractures
figure;surf(reshape(coordinates(:,1),Cell_N+1,Cell_M+1),reshape(coordinates(:,2),Cell_N+1,Cell_M+1),zeros(Cell_N+1,Cell_M+1));
colormap('white');axis image;view([0,90]);hold on
plot(ParametersofFracture(:,[7,9])',ParametersofFracture(:,[8,10])',...
   'r-','LineWidth',1.5); hold off; drawnow;
% (2) compute basis data
HighOrderDegree = 1 ; % bilinear Lagrange basis 
HighOrderNp = (HighOrderDegree+1)^2; % degree of freedom in each element for high order polynomial
Gauss_num = 6 ; % the number of Gauss quadrature point in [-1,1]
Lobatto_num = 6 ; % the number of Lobatto quadrature point in [-1,1]
[Gauss_x,Gauss_w] = JacobiGQ(0,0,Gauss_num-1); 
Gauss_x = 0 + ( Gauss_x - (-1) ) / 2 ; % Gauss_x now is in [0,1].
Gauss_w = Gauss_w / 2 ;% sum(Gauss_w) = 1
% Tensor product of 1D quadrature rule
[quad_lamda1,quad_lamda2]=meshgrid(Gauss_x,Gauss_x); 
quad_lamda1 = quad_lamda1(:);
quad_lamda2 = quad_lamda2(:);
quad_w = Gauss_w*Gauss_w';
quad_w = quad_w(:);
basisP = zeros( HighOrderNp , size(quad_w,1) ) ; % basis function's value in quadrature point
GradrbasisP = zeros( HighOrderNp , size(quad_w,1) ) ; % gradient's value in quadrature point
GradsbasisP = zeros( HighOrderNp , size(quad_w,1) ) ; % gradient's value in quadrature point
boundaryP = zeros(HighOrderNp,Gauss_num,4);
GradrboundaryP = zeros(HighOrderNp,Gauss_num,4);
GradsboundaryP = zeros(HighOrderNp,Gauss_num,4);
% compute the basis function and it's gradient's value on quadrature point:
basisP(1,:) = (1-quad_lamda1).*(1-quad_lamda2) ;
basisP(2,:) =    quad_lamda1 .*(1-quad_lamda2) ;
basisP(3,:) =    quad_lamda1 .*   quad_lamda2  ;
basisP(4,:) = (1-quad_lamda1).*   quad_lamda2 ;
GradrbasisP(1,:) = -(1-quad_lamda2); 
GradrbasisP(2,:) =  (1-quad_lamda2); 
GradrbasisP(3,:) =     quad_lamda2 ; 
GradrbasisP(4,:) =    -quad_lamda2 ; 
GradsbasisP(1,:) = -(1-quad_lamda1);
GradsbasisP(2,:) =    -quad_lamda1;
GradsbasisP(3,:) =     quad_lamda1;
GradsbasisP(4,:) =   1-quad_lamda1;
boundaryP(1,:,1) =  1-Gauss_x; GradrboundaryP(1,:,1) = -1; GradsboundaryP(1,:,1) = -(1-Gauss_x); % lower edge 
boundaryP(1,:,2) =  0; GradrboundaryP(1,:,2) = -(1-Gauss_x); GradsboundaryP(1,:,2) = 0; % right edge
boundaryP(1,:,3) =  0; GradrboundaryP(1,:,3) = 0; GradsboundaryP(1,:,3) = -Gauss_x; % upper edge 
boundaryP(1,:,4) =  Gauss_x; GradrboundaryP(1,:,4) = -Gauss_x; GradsboundaryP(1,:,4) = -1; % left edge 
boundaryP(2,:,1) =  Gauss_x; GradrboundaryP(1,:,1) = 1; GradsboundaryP(1,:,1) = -Gauss_x; % lower edge 
boundaryP(2,:,2) =  1-Gauss_x; GradrboundaryP(1,:,2) = 1-Gauss_x; GradsboundaryP(1,:,2) = -1; % right edge
boundaryP(2,:,3) =  0; GradrboundaryP(1,:,3) = 0; GradsboundaryP(1,:,3) = -(1-Gauss_x); % upper edge 
boundaryP(2,:,4) =  0; GradrboundaryP(1,:,4) = Gauss_x; GradsboundaryP(1,:,4) = 0; % left edge 
boundaryP(3,:,1) =  0; GradrboundaryP(1,:,1) = 0; GradsboundaryP(1,:,1) = Gauss_x; % lower edge 
boundaryP(3,:,2) =  Gauss_x; GradrboundaryP(1,:,2) = Gauss_x; GradsboundaryP(1,:,2) = 1; % right edge
boundaryP(3,:,3) =  1-Gauss_x; GradrboundaryP(1,:,3) = 1; GradsboundaryP(1,:,3) = 1-Gauss_x; % upper edge 
boundaryP(3,:,4) =  0; GradrboundaryP(1,:,4) = 1-Gauss_x; GradsboundaryP(1,:,4) = 0; % left edge 
boundaryP(4,:,1) =  0; GradrboundaryP(1,:,1) = 0; GradsboundaryP(1,:,1) = (1-Gauss_x); % lower edge 
boundaryP(4,:,2) =  0; GradrboundaryP(1,:,2) = -Gauss_x; GradsboundaryP(1,:,2) = 0; % right edge
boundaryP(4,:,3) =  Gauss_x; GradrboundaryP(1,:,3) = -1; GradsboundaryP(1,:,3) = Gauss_x; % upper edge 
boundaryP(4,:,4) =  1-Gauss_x; GradrboundaryP(1,:,4) = -(1-Gauss_x); GradsboundaryP(1,:,4) = 1; % left edge 
% compute fractures and Barriers basic data
[ParametersofFracture,FracturesPath,FracturesLength,basisPonFracture,GradnubasisPonFracture,GradsigmabasisPonFracture] = ...
    SetFractureBasisData1(coordinates,elements4,ParametersofFracture,Gauss_x);
% define stiffness matrix and right hand side vector
b = zeros(size(coordinates,1),1);% right hand side of the scheme
u = zeros(size(coordinates,1),1);% solution vector
IK = zeros(size( elements4,2 )*HighOrderNp^2+size(FracturesPath,1)*HighOrderNp^2,1);
JK = zeros(size( elements4,2 )*HighOrderNp^2+size(FracturesPath,1)*HighOrderNp^2,1);
VK = zeros(size( elements4,2 )*HighOrderNp^2+size(FracturesPath,1)*HighOrderNp^2,1);
for jj = 1 : size( elements4,2 ) % assembly K & b 
xy = [coordinates(elements4(1,jj),1)+hx*quad_lamda1,coordinates(elements4(1,jj),2)+hy*quad_lamda2];
% assembly stiffness matrix
stima = Tarea * ( ... 
    (GradrbasisP*Jacobimat(1,1)+GradsbasisP*Jacobimat(2,1))...
        *diag(km(xy).*quad_w)...
        *(GradrbasisP*Jacobimat(1,1)+GradsbasisP*Jacobimat(2,1))' ...
      + ...
    (GradrbasisP*Jacobimat(1,2)+GradsbasisP*Jacobimat(2,2))...
        *diag(km(xy).*quad_w)...
        *(GradrbasisP*Jacobimat(1,2)+GradsbasisP*Jacobimat(2,2))' ) ; 
[J_e,I_e] = meshgrid(elements4(:,jj)); I_index = I_e(:); J_index = J_e(:);
IK((jj-1)*HighOrderNp^2+1:jj*HighOrderNp^2)=I_index;
JK((jj-1)*HighOrderNp^2+1:jj*HighOrderNp^2)=J_index;
VK((jj-1)*HighOrderNp^2+1:jj*HighOrderNp^2)=stima(:);
% assembly right hand side matrix 
b(elements4(:,jj)) = b(elements4(:,jj)) + Tarea*basisP*(f(xy).*quad_w) ;
% Neumann conditions
for ii = 1 : 4  % the ii-th edge of jj-th element
    vertex1 = ii ; 
    vertex2 = mod(ii,4)+1 ;
    edge_vec =  coordinates(elements4(vertex2,jj),:)' - ...
        coordinates(elements4(vertex1,jj),:)' ;
    outer_vec = - [0,-1;1,0]*(edge_vec/norm(edge_vec));
    E = EtoEmap(ii,jj) ; % the index of neighborhood of element jj
    if E == -1
        if abs(outer_vec(2))>1/2 % up and down boundary, Neumann boundary condition
        neumann_count = neumann_count + 1; neumann(neumann_count,:) = elements4([vertex1;vertex2],jj) ;   
        Boundary_xy = (1-Gauss_x)*coordinates(elements4(vertex1,jj),:) + Gauss_x*coordinates(elements4(vertex2,jj),:);
        b(elements4([vertex1;vertex2],jj)) = b(elements4([vertex1;vertex2],jj)) + ...
            Elenth(ii)*boundaryP([vertex1;vertex2],:,ii)*(g(Boundary_xy).*Gauss_w);
        end
    end
end
end

for k = 1 : NumberofFractures
    for ell = 1 : ParametersofFracture(k,11)
    CurrentIndex = sum(ParametersofFracture(1:k-1,11)) + ell;
    jj = FracturesPath( CurrentIndex );
    [J_e,I_e] = meshgrid(elements4(:,jj));I_index = I_e(:); J_index = J_e(:); 
    stima = ParametersofFracture(k,5)*ParametersofFracture(k,6)*FracturesLength( CurrentIndex )*...
        GradnubasisPonFracture((CurrentIndex-1)*4+1:CurrentIndex*4, 1:end-2)*diag(Gauss_w)*...
            GradnubasisPonFracture((CurrentIndex-1)*4+1:CurrentIndex*4, 1:end-2)';
    IK(size( elements4,2 )*HighOrderNp^2 + ( (CurrentIndex-1)*HighOrderNp^2+1:CurrentIndex*HighOrderNp^2 )) = I_index;
    JK(size( elements4,2 )*HighOrderNp^2 + ( (CurrentIndex-1)*HighOrderNp^2+1:CurrentIndex*HighOrderNp^2 )) = J_index;
    VK(size( elements4,2 )*HighOrderNp^2 + ( (CurrentIndex-1)*HighOrderNp^2+1:CurrentIndex*HighOrderNp^2 )) = stima(:);
   end
end
K = sparse(IK,JK,VK,size(coordinates,1),size(coordinates,1));
dirichlet = setdiff(boundary_edge, neumann, 'rows') ; 
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
% Dirichlet conditions
u(unique(dirichlet)) = p(coordinates(unique(dirichlet),:));
b = b - K * u;
% Computation of the solution
u(FreeNodes) = K(FreeNodes,FreeNodes) \ b(FreeNodes);

% post-processing----------------------------------------------------------
X_meshgrid = reshape(coordinates(:,1),Cell_N+1,Cell_M+1);
Y_meshgrid = reshape(coordinates(:,2),Cell_N+1,Cell_M+1);
U_meshgrid = reshape(u,Cell_N+1,Cell_M+1);

figure; 
surf(reshape(coordinates(:,1),Cell_N+1,Cell_M+1),reshape(coordinates(:,2),Cell_N+1,Cell_M+1),reshape(u,Cell_N+1,Cell_M+1))
shading interp;colormap('jet');caxis([0,1013250]);cb=colorbar;cb.Position = cb.Position.*[1.07,1,1,1];
hold on; plot3(ParametersofFracture(:,[7,9])',ParametersofFracture(:,[8,10])',(max(u)+1)*ones(size(ParametersofFracture(:,[8,10])')),...
   'w-','LineWidth',0.05); hold off; axis equal; axis off;ax=axis;axis(ax*1.001); view([0,90]);
drawnow

dof = length(FreeNodes)
sparsity = nnz(K(FreeNodes,FreeNodes))/dof^2*1000
conditionalnumber = condest(K(FreeNodes,FreeNodes))

% -------------------------------------------------------------------------
function km = km(xy)
x = xy(:,1);
y = xy(:,2);
km = (1e-14)*ones(size(x)); 
end
function f = f(xy)
x = xy(:,1);
y = xy(:,2);
%f = 2*cos(x).*cos(y);
f = zeros(size(x));
end
function p = p(xy)
x = xy(:,1);
y = xy(:,2);
p = ((700-x)/700)*1013250; 
end
function g = g(xy)
x = xy(:,1);
y = xy(:,2);
g = zeros(size(x)); 
end

end

