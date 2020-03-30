function CircleFractures
% time : 2020.1.26 - 2020.1.26
% author : xuziyao
% this program is code for Continuous Finite Element Method for
% Poisson equation with delta diffusion coefficient on rectangular meshes with
% curved 2d fractures.
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
xmin = 0 ; xmax = 1 ; % domain 
ymin = 0 ; ymax = 1 ; % domain
% x direction is divided into N parts ,y direction is divided into M parts:
Cell_N = 60 ; Cell_M =  Cell_N; 
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ;
[coordinates,elements4,EtoEmap] = RectangularMesh( xmin,xmax,ymin,ymax,Cell_N,Cell_M );% periodic boundary condition
EtoEmap = Rectangulation_neighbor_rectangles( elements4 ); 
Tarea =  hx*hy; Elenth = [hx;hy;hx;hy];
Jacobimat = [1/hx,0;0,1/hy];%Jacobimat=[Dr/Dx,Dr/Dy ; Ds/Dx,Ds/Dy] 
boundary_edge = boundedges_Rec( elements4 ) ;
neumann = zeros(size(boundary_edge)) ; neumann_count = 0;  
figure;surf(reshape(coordinates(:,1),Cell_N+1,Cell_M+1),reshape(coordinates(:,2),Cell_N+1,Cell_M+1),zeros(Cell_N+1,Cell_M+1));
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
% Compute curved fractures and Barriers basic data
[NumberofFractures,NumberofElementsPassed,FracturesPath,FracturesLength,basisPonFracture,GradnubasisPonFracture,GradsigmabasisPonFracture,...
    PermeabilitynuonFracture,ThicknessonFracture,ds_dtonFracture] = ...
    SetCurvedFractureBasisData(coordinates,elements4,Gauss_x);
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
        if abs(outer_vec(2))>0.5 %  up and down boundary, Neumann boundary condition
        neumann_count = neumann_count + 1; neumann(neumann_count,:) = elements4([vertex1;vertex2],jj) ;   
        Boundary_xy = (1-Gauss_x)*coordinates(elements4(vertex1,jj),:) + Gauss_x*coordinates(elements4(vertex2,jj),:);
        b(elements4([vertex1;vertex2],jj)) = b(elements4([vertex1;vertex2],jj)) + ...
            Elenth(ii)*boundaryP([vertex1;vertex2],:,ii)*(g(Boundary_xy).*Gauss_w);
        end
    end
end

end

for k = 1 : NumberofFractures
    for ell = 1 : NumberofElementsPassed(k)
    CurrentIndex = sum(NumberofElementsPassed(1:k-1)) + ell;
    jj = FracturesPath( CurrentIndex );
    [J_e,I_e] = meshgrid(elements4(:,jj));I_index = I_e(:); J_index = J_e(:); 
    stima = FracturesLength( CurrentIndex )*...
        GradnubasisPonFracture((CurrentIndex-1)*4+1:CurrentIndex*4, 1:end-2)*...
        diag(Gauss_w.*ThicknessonFracture(CurrentIndex,1:end-2)'.*PermeabilitynuonFracture(CurrentIndex,1:end-2)'.*ds_dtonFracture(CurrentIndex,1:end-2)')*...
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
%dof = length(FreeNodes)
%conditionalnumber=condest(K(FreeNodes,FreeNodes))
X_meshgrid = reshape(coordinates(:,1),Cell_M+1,Cell_N+1);
Y_meshgrid = reshape(coordinates(:,2),Cell_M+1,Cell_N+1);
U_meshgrid = reshape(u,Cell_M+1,Cell_N+1);
close all;
hold on;
surf(X_meshgrid,Y_meshgrid,U_meshgrid)
colormap(parula);xlim([0,1]); ylim([0 1]);view(50,30);%shading flat;%caxis([-100000, 1200000]); 
contour(X_meshgrid,Y_meshgrid,U_meshgrid)
hold off;

save('data')

[X_slice0_90, Y_slice0_90] = meshgrid(linspace(xmin*(1+1e-7),xmax*(1-1e-7),1000),0.45);
R_slice0_90 = sqrt((X_slice0_90 - X_slice0_90(1)).^2+(Y_slice0_90-Y_slice0_90(1)).^2);
U_slice0_90 = interp2(X_meshgrid',Y_meshgrid',U_meshgrid',X_slice0_90,Y_slice0_90,'linear');
save('Data','X_slice0_90','Y_slice0_90','R_slice0_90','U_slice0_90','-append') 
figure;plot(R_slice0_90,U_slice0_90,'LineWidth',1.5);

[X_slice1_90, Y_slice1_90] = meshgrid(linspace(xmin*(1+1e-7),xmax*(1-1e-7),1000),0.5);
R_slice1_90 = sqrt((X_slice1_90 - X_slice1_90(1)).^2+(Y_slice1_90-Y_slice1_90(1)).^2);
U_slice1_90 = interp2(X_meshgrid',Y_meshgrid',U_meshgrid',X_slice1_90,Y_slice1_90,'linear');
save('Data','X_slice1_90','Y_slice1_90','R_slice1_90','U_slice1_90','-append') 
figure;plot(R_slice1_90,U_slice1_90,'LineWidth',1.5);

[X_slice2_90, Y_slice2_90] = meshgrid(0.3,linspace(ymin*(1+1e-7),ymax*(1-1e-7),1000));
R_slice2_90 = sqrt((X_slice2_90 - X_slice2_90(1)).^2+(Y_slice2_90-Y_slice2_90(1)).^2);
U_slice2_90 = interp2(X_meshgrid',Y_meshgrid',U_meshgrid',X_slice2_90,Y_slice2_90,'linear');
save('Data','X_slice2_90','Y_slice2_90','R_slice2_90','U_slice2_90','-append') 
figure;plot(R_slice2_90,U_slice2_90,'LineWidth',1.5);

[X_slice3_90, Y_slice3_90] = meshgrid(0.4,linspace(ymin*(1+1e-7),ymax*(1-1e-7),1000));
R_slice3_90 = sqrt((X_slice3_90 - X_slice3_90(1)).^2+(Y_slice3_90-Y_slice3_90(1)).^2);
U_slice3_90 = interp2(X_meshgrid',Y_meshgrid',U_meshgrid',X_slice3_90,Y_slice3_90,'linear');
save('Data','X_slice3_90','Y_slice3_90','R_slice3_90','U_slice3_90','-append') 
figure;plot(R_slice3_90,U_slice3_90,'LineWidth',1.5);

[X_slice4_90, Y_slice4_90] = meshgrid(0.2,linspace(ymin*(1+1e-7),ymax*(1-1e-7),1000));
R_slice4_90 = sqrt((X_slice4_90 - X_slice4_90(1)).^2+(Y_slice4_90-Y_slice4_90(1)).^2);
U_slice4_90 = interp2(X_meshgrid',Y_meshgrid',U_meshgrid',X_slice4_90,Y_slice4_90,'linear');
save('Data','X_slice4_90','Y_slice4_90','R_slice4_90','U_slice4_90','-append') 
figure;plot(R_slice4_90,U_slice4_90,'LineWidth',1.5);


% interpolate the data (then draw it):
%[X_slice, Y_slice] = meshgrid(linspace(xmin*(1+1e-7),xmax*(1-1e-7),1000),500);
%U_slice = interp2(reshape(coordinates(:,1),Cell_N+1,Cell_M+1)',reshape(coordinates(:,2),Cell_N+1,Cell_M+1)'...
%    ,reshape(u,Cell_N+1,Cell_M+1)',X_slice,Y_slice,'linear');
%figure;plot(X_slice,U_slice,'LineWidth',1.5); xlim([0,700]); ylim([0, 1200000]);
%hold on;
%hold off;




%[X_slice, Y_slice] = meshgrid(625,linspace(ymin*(1+1e-7),ymax*(1-1e-7),1000));
%U_slice = interp2(reshape(coordinates(:,1),Cell_N+1,Cell_M+1)',reshape(coordinates(:,2),Cell_N+1,Cell_M+1)'...
%    ,reshape(u,Cell_N+1,Cell_M+1)',X_slice,Y_slice,'linear');
%figure;plot(Y_slice,U_slice,'LineWidth',1.5); xlim([0,600]); ylim([0, 900000]);
%hold on;
%hold off;


% interpolate the data (then draw it):
%Interpx = linspace( xmin+1e-4,xmax-1e-4,100);Interpy = linspace( ymin+1e-4,ymax-1e-4,100);[XX,YY] = meshgrid(Interpx,Interpy);
%Interpu_FEM = tri2grid( coordinates' , elements4 , u , Interpx , Interpy ); 
%figure; h=surf(XX,YY,Interpu_FEM);view(10,40);xlim([xmin,xmax]); ylim([ymin ymax]) %set(h, 'LineStyle','none');
%save('FEM_DeltaFractureSteadyState','Interpu_FEM','hx')


% -------------------------------------------------------------------------
function km = km(xy)
x = xy(:,1);
y = xy(:,2);
km = ones(size(x)); 
end
function f = f(xy)
x = xy(:,1);
y = xy(:,2);
f = zeros(size(x));
end
function p = p(xy)
x = xy(:,1);
y = xy(:,2);
p = 1-x; 
end
function g = g(xy)
x = xy(:,1);
y = xy(:,2);
g = zeros(size(x));
end

end