function SingleFracture3D
% time : 2020.2.15 - 2020.3.8
% author : xuziyao
% this program is code for Continuous Finite Element Method for
% Poisson equation on cubic meshes with fracture faces.
% PDE model: 
% -div(D*grad(p)) = f, 
% Neumann & Dirichlet boundary condition 
% left and right boundary: Dirichlet
% others: Neumann 
% Basis function is the trilinear Lagrange basis on the (r,s,t)\in [0,1]*[0,1]*[0,1] reference cube.
% phi_1(r,s,t) = (1-r)(1-s)(1-t)
% phi_2(r,s,t) = r(1-s)(1-t)
% phi_3(r,s,t) = rs(1-t)
% phi_4(r,s,t) = (1-r)s(1-t)
% phi_5(r,s,t) = (1-r)(1-s)t
% phi_6(r,s,t) = r(1-s)t
% phi_7(r,s,t) = rst
% phi_8(r,s,t) = (1-r)st
% Gradrphi_1 = -(1-s)(1-t); Gradsphi_1 = -(1-r)(1-t); Gradtphi_1 = -(1-r)(1-s)
% Gradrphi_2 =  (1-s)(1-t); Gradsphi_2 =     -r(1-t); Gradtphi_2 = -r(1-s)
% Gradrphi_3 =    s  (1-t); Gradsphi_3 =      r(1-t); Gradtphi_3 = -rs
% Gradrphi_4 =   -s  (1-t); Gradsphi_4 =  (1-r)(1-t); Gradtphi_4 = -(1-r)s
% Gradrphi_5 = -(1-s)t; Gradsphi_5 = -(1-r)t; Gradtphi_5 = (1-r)(1-s)
% Gradrphi_6 =  (1-s)t; Gradsphi_6 =     -rt; Gradtphi_6 = r(1-s)
% Gradrphi_7 =    s  t; Gradsphi_7 =      rt; Gradtphi_7 = rs
% Gradrphi_8 =   -s  t; Gradsphi_8 =  (1-r)t; Gradtphi_8 = (1-r)s
% FINISHED
clc,clear
format long
% set the geometry and triangulation parameters:
xmin = 0 ; xmax = 1 ; % domain 
ymin = 0 ; ymax = 1 ; % domain
zmin = 0 ; zmax = 1 ; % domain
% x direction is divided into N parts ,y direction is divided into M parts, z direction is divided into L parts
Cell_N = 50 ; Cell_M = Cell_N ; Cell_L = Cell_N; 
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ; hz = ( zmax - zmin ) / Cell_L ;
[coordinates,elements,EtoEmap] = CubicMesh( xmin,xmax,ymin,ymax,zmin,zmax,Cell_N,Cell_M,Cell_L );% nonperiodic/periodic boundary condition
Tvolumn =  hx*hy*hz; Farea = [hx*hy;hx*hy;hy*hz;hy*hz;hz*hx;hz*hx]; 
Jacobimat = [1/hx,0,0 ; 0,1/hy,0 ; 0,0,1/hz];%Jacobimat=[Dr/Dx,Dr/Dy,Dr/Dz ; Ds/Dx,Ds/Dy,Ds/Dz ; Dt/Dx,Dt/Dy,Dt/Dz] 
FacetoLocalNodes = [1,2,3,4; 5,6,7,8; 1,4,5,8; 2,3,6,7; 1,2,5,6; 3,4,7,8];
neumann = zeros( 2*( Cell_N*Cell_M + Cell_M*Cell_L + Cell_L*Cell_N) , 4) ; neumann_count = 0;
boundary_faces = zeros( 2*( Cell_N*Cell_M + Cell_M*Cell_L + Cell_L*Cell_N) , 4) ; boundary_count = 0;
% ParametersofFracture: width,permeability,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,sigma1,sigma2,sigma3,...
%                       # of elements passed by this fracture. 
NumberofFractures = 1; 
ParametersofFracture = zeros(NumberofFractures, 18); 
ParametersofFracture(:,1) = 1e-3; % width of fractures
ParametersofFracture(:,2) = 1e8; % permeability of fractures
ParametersofFracture(:,3:14) = [0.25,0.2,0.1,  0.25,0.7,0.7,  0.75,0.7,0.7,  0.75,0.2,0.1]...
                              + (rand(NumberofFractures,12)-0)*1e-5; % position of fractures
% (2) compute basis data
HighOrderDegree = 1 ; % bilinear Lagrange basis 
HighOrderNp = (HighOrderDegree+1)^3; % degree of freedom in each element for high order polynomial
Gauss_num = 3 ; % the number of Gauss quadrature point in [-1,1]
Lobatto_num = 6 ; % the number of Lobatto quadrature point in [-1,1]
[Gauss_x,Gauss_w] = JacobiGQ(0,0,Gauss_num-1); 
Gauss_x = 0 + ( Gauss_x - (-1) ) / 2 ; % Gauss_x now is in [0,1].
Gauss_w = Gauss_w / 2 ;% sum(Gauss_w) = 1
% Tensor product of 1D quadrature rule
[quad_lamda1,quad_lamda2,quad_lamda3]=meshgrid(Gauss_x,Gauss_x,Gauss_x); 
quad_lamda1 = quad_lamda1(:);
quad_lamda2 = quad_lamda2(:);
quad_lamda3 = quad_lamda3(:);
quad_w = kron(Gauss_w',Gauss_w*Gauss_w');
quad_w = quad_w(:);
basisP = zeros( HighOrderNp , size(quad_w,1) ) ; % basis function's value in quadrature point
GradrbasisP = zeros( HighOrderNp , size(quad_w,1) ) ; % gradient's value in quadrature point
GradsbasisP = zeros( HighOrderNp , size(quad_w,1) ) ; % gradient's value in quadrature point
GradtbasisP = zeros( HighOrderNp , size(quad_w,1) ) ; % gradient's value in quadrature point
boundaryP = zeros(HighOrderNp,Gauss_num^2,6);
GradrboundaryP = zeros(HighOrderNp,Gauss_num^2,6);
GradsboundaryP = zeros(HighOrderNp,Gauss_num^2,6);
GradtboundaryP = zeros(HighOrderNp,Gauss_num^2,6);
% compute the basis function and it's gradient's value on quadrature point:
basisP(1,:) = (1-quad_lamda1).*(1-quad_lamda2).*(1-quad_lamda3) ;
basisP(2,:) =    quad_lamda1 .*(1-quad_lamda2).*(1-quad_lamda3) ;
basisP(3,:) =    quad_lamda1 .*   quad_lamda2 .*(1-quad_lamda3) ;
basisP(4,:) = (1-quad_lamda1).*   quad_lamda2 .*(1-quad_lamda3) ;
basisP(5,:) = (1-quad_lamda1).*(1-quad_lamda2).*quad_lamda3 ;
basisP(6,:) =    quad_lamda1 .*(1-quad_lamda2).*quad_lamda3 ;
basisP(7,:) =    quad_lamda1 .*   quad_lamda2 .*quad_lamda3 ;
basisP(8,:) = (1-quad_lamda1).*   quad_lamda2 .*quad_lamda3 ;
GradrbasisP(1,:) = -(1-quad_lamda2).*(1-quad_lamda3); 
GradrbasisP(2,:) =  (1-quad_lamda2).*(1-quad_lamda3); 
GradrbasisP(3,:) =     quad_lamda2 .*(1-quad_lamda3) ; 
GradrbasisP(4,:) =    -quad_lamda2 .*(1-quad_lamda3) ; 
GradrbasisP(5,:) = -(1-quad_lamda2).*quad_lamda3; 
GradrbasisP(6,:) =  (1-quad_lamda2).*quad_lamda3; 
GradrbasisP(7,:) =     quad_lamda2 .*quad_lamda3 ; 
GradrbasisP(8,:) =    -quad_lamda2 .*quad_lamda3; 

GradsbasisP(1,:) = -(1-quad_lamda1).*(1-quad_lamda3);
GradsbasisP(2,:) =    -quad_lamda1.*(1-quad_lamda3);
GradsbasisP(3,:) =     quad_lamda1.*(1-quad_lamda3);
GradsbasisP(4,:) =   (1-quad_lamda1).*(1-quad_lamda3);
GradsbasisP(5,:) = -(1-quad_lamda1).*quad_lamda3;
GradsbasisP(6,:) =    -quad_lamda1.*quad_lamda3;
GradsbasisP(7,:) =     quad_lamda1.*quad_lamda3;
GradsbasisP(8,:) =  (1-quad_lamda1).*quad_lamda3;

GradtbasisP(1,:) = -(1-quad_lamda1).*(1-quad_lamda2);
GradtbasisP(2,:) = -quad_lamda1.*(1-quad_lamda2);
GradtbasisP(3,:) = -quad_lamda1.*quad_lamda2;
GradtbasisP(4,:) = -(1-quad_lamda1).*quad_lamda2;
GradtbasisP(5,:) = (1-quad_lamda1).*(1-quad_lamda2);
GradtbasisP(6,:) = quad_lamda1.*(1-quad_lamda2);
GradtbasisP(7,:) = quad_lamda1.*quad_lamda2;
GradtbasisP(8,:) = (1-quad_lamda1).*quad_lamda2;

face_w = Gauss_w*Gauss_w'; face_w=face_w(:);
quad_face_r = zeros(6,Gauss_num^2);
quad_face_s = zeros(6,Gauss_num^2);
quad_face_t = zeros(6,Gauss_num^2);
for ii = 1 : 6
if ii==1
[face_r,face_s,face_t] = meshgrid(Gauss_x,Gauss_x,0);    
elseif ii==2
[face_r,face_s,face_t] = meshgrid(Gauss_x,Gauss_x,1);
elseif ii==3
[face_r,face_s,face_t] = meshgrid(0,Gauss_x,Gauss_x);    
elseif ii==4
[face_r,face_s,face_t] = meshgrid(1,Gauss_x,Gauss_x);    
elseif ii==5
[face_r,face_s,face_t] = meshgrid(Gauss_x,0,Gauss_x);    
elseif ii==6
[face_r,face_s,face_t] = meshgrid(Gauss_x,1,Gauss_x);    
end
face_r=face_r(:);face_s=face_s(:);face_t=face_t(:);
quad_face_r(ii,:) = face_r;
quad_face_s(ii,:) = face_s;
quad_face_t(ii,:) = face_t;
boundaryP(1,:,ii) = (1-face_r).*(1-face_s).*(1-face_t);
boundaryP(2,:,ii) = face_r.*(1-face_s).*(1-face_t);
boundaryP(3,:,ii) = face_r.*face_s.*(1-face_t);
boundaryP(4,:,ii) = (1-face_r).*face_s.*(1-face_t);
boundaryP(5,:,ii) = (1-face_r).*(1-face_s).*face_t;
boundaryP(6,:,ii) = face_r.*(1-face_s).*face_t;
boundaryP(7,:,ii) = face_r.*face_s.*face_t;
boundaryP(8,:,ii) = (1-face_r).*face_s.*face_t;

GradrboundaryP(1,:,ii) = -(1-face_s).*(1-face_t);
GradrboundaryP(2,:,ii) = (1-face_s).*(1-face_t);
GradrboundaryP(3,:,ii) = face_s.*(1-face_t);
GradrboundaryP(4,:,ii) = -face_s.*(1-face_t);
GradrboundaryP(5,:,ii) = -(1-face_s).*face_t;
GradrboundaryP(6,:,ii) = (1-face_s).*face_t;
GradrboundaryP(7,:,ii) = face_s.*face_t;
GradrboundaryP(8,:,ii) = -face_s.*face_t;

GradsboundaryP(1,:,ii) = -(1-face_r).*(1-face_t);
GradsboundaryP(2,:,ii) = -face_r.*(1-face_t);
GradsboundaryP(3,:,ii) = face_r.*(1-face_t);
GradsboundaryP(4,:,ii) = (1-face_r).*(1-face_t);
GradsboundaryP(5,:,ii) = -(1-face_r).*face_t;
GradsboundaryP(6,:,ii) = -face_r.*face_t;
GradsboundaryP(7,:,ii) = face_r.*face_t;
GradsboundaryP(8,:,ii) = (1-face_r).*face_t;

GradtboundaryP(1,:,ii) = -(1-face_r).*(1-face_s);
GradtboundaryP(2,:,ii) = -face_r.*(1-face_s);
GradtboundaryP(3,:,ii) = -face_r.*face_s;
GradtboundaryP(4,:,ii) = -(1-face_r).*face_s;
GradtboundaryP(5,:,ii) = (1-face_r).*(1-face_s);
GradtboundaryP(6,:,ii) = face_r.*(1-face_s);
GradtboundaryP(7,:,ii) = face_r.*face_s;
GradtboundaryP(8,:,ii) = (1-face_r).*face_s;
end
% compute fractures and Barriers basic data
[ParametersofFracture,FracturesPath,FracturesArea,basisPonFracture, ...
    GradxbasisPonFracture,GradybasisPonFracture,GradzbasisPonFracture,GradsigmabasisPonFracture] = ...
    SetFractureBasisData3d(coordinates,elements,ParametersofFracture,Gauss_x);
% define stiffness matrix and right hand side vector
b = zeros(size(coordinates,1),1);% right hand side of the scheme
u = zeros(size(coordinates,1),1);% solution vector
IK = zeros(size( elements,2 )*HighOrderNp^2 + size(FracturesPath,1)*HighOrderNp^2,1);
JK = zeros(size( elements,2 )*HighOrderNp^2 + size(FracturesPath,1)*HighOrderNp^2,1);
VK = zeros(size( elements,2 )*HighOrderNp^2 + size(FracturesPath,1)*HighOrderNp^2,1);
for jj = 1 : size( elements,2 ) % assembly K & b 
    jj
xyz = [coordinates(elements(1,jj),1)+hx*quad_lamda1,coordinates(elements(1,jj),2)+hy*quad_lamda2,coordinates(elements(1,jj),3)+hz*quad_lamda3];
% assembly stiffness matrix
stima = Tvolumn * ( ... 
    (GradrbasisP*Jacobimat(1,1)+GradsbasisP*Jacobimat(2,1)+GradtbasisP*Jacobimat(3,1))...
        *diag(km(xyz).*quad_w)...
        *(GradrbasisP*Jacobimat(1,1)+GradsbasisP*Jacobimat(2,1)+GradtbasisP*Jacobimat(3,1))' ...
      + ...
    (GradrbasisP*Jacobimat(1,2)+GradsbasisP*Jacobimat(2,2)+GradtbasisP*Jacobimat(3,2))...
        *diag(km(xyz).*quad_w)...
        *(GradrbasisP*Jacobimat(1,2)+GradsbasisP*Jacobimat(2,2)+GradtbasisP*Jacobimat(3,2))'  ...
      + ...
    (GradrbasisP*Jacobimat(1,3)+GradsbasisP*Jacobimat(2,3)+GradtbasisP*Jacobimat(3,3))...
        *diag(km(xyz).*quad_w)...
        *(GradrbasisP*Jacobimat(1,3)+GradsbasisP*Jacobimat(2,3)+GradtbasisP*Jacobimat(3,3))' ) ;
[J_e,I_e] = meshgrid(elements(:,jj)); I_index = I_e(:); J_index = J_e(:);
IK((jj-1)*HighOrderNp^2+1:jj*HighOrderNp^2)=I_index;
JK((jj-1)*HighOrderNp^2+1:jj*HighOrderNp^2)=J_index;
VK((jj-1)*HighOrderNp^2+1:jj*HighOrderNp^2)=stima(:);
% assembly right hand side matrix 
b(elements(:,jj)) = b(elements(:,jj)) + Tvolumn*basisP*(f(xyz).*quad_w) ;
% Neumann conditions
for ii = 1 : 6  % the ii-th face of jj-th element
    if ii==1
        outer_vec = [0;0;-1];
    elseif ii==2
        outer_vec = [0;0;1];
    elseif ii==3
        outer_vec = [-1;0;0];
    elseif ii==4
        outer_vec = [1;0;0];
    elseif ii==5
        outer_vec = [0;-1;0];
    else
        outer_vec = [0;1;0];
    end
    E = EtoEmap(ii,jj) ; % the index of neighborhood of element jj
    if E == -1
        boundary_count = boundary_count + 1; boundary_faces(boundary_count,:) = elements(FacetoLocalNodes(ii,:),jj) ;
        if abs(outer_vec(1))<0.5 %  Neumann boundary condition
        neumann_count = neumann_count + 1; neumann(neumann_count,:) = elements(FacetoLocalNodes(ii,:),jj) ;   
        Boundary_xyz = ...
            [coordinates(elements(1,jj),1)+hx*quad_face_r(ii,:)',coordinates(elements(1,jj),2)+hy*quad_face_s(ii,:)',coordinates(elements(1,jj),3)+hz*quad_face_t(ii,:)'];
        b(elements(FacetoLocalNodes(ii,:),jj)) = b(elements(FacetoLocalNodes(ii,:),jj)) + ...
            Farea(ii)*boundaryP(FacetoLocalNodes(ii,:),:,ii)*(g(Boundary_xyz,outer_vec).*face_w);
        end
    end
end
end
% fractures 
for k = 1 : NumberofFractures
    for ell = 1 : ParametersofFracture(k,18)
    CurrentIndex = sum(ParametersofFracture(1:k-1,18)) + ell;
    jj = FracturesPath( CurrentIndex );
    [J_e,I_e] = meshgrid(elements(:,jj)); I_index = I_e(:); J_index = J_e(:); 
    stima = ParametersofFracture(k,1)*ParametersofFracture(k,2)*FracturesArea( CurrentIndex ) * ( ...
         GradxbasisPonFracture((CurrentIndex-1)*8+[1:8],:)*diag(face_w)*GradxbasisPonFracture((CurrentIndex-1)*8+[1:8],:)' ...
       + GradybasisPonFracture((CurrentIndex-1)*8+[1:8],:)*diag(face_w)*GradybasisPonFracture((CurrentIndex-1)*8+[1:8],:)' ...
       + GradzbasisPonFracture((CurrentIndex-1)*8+[1:8],:)*diag(face_w)*GradzbasisPonFracture((CurrentIndex-1)*8+[1:8],:)' ...
       - GradsigmabasisPonFracture((CurrentIndex-1)*8+[1:8],:)*diag(face_w)*GradsigmabasisPonFracture((CurrentIndex-1)*8+[1:8],:)');
    IK(size( elements,2 )*HighOrderNp^2 + ( (CurrentIndex-1)*HighOrderNp^2+1:CurrentIndex*HighOrderNp^2 )) = I_index;
    JK(size( elements,2 )*HighOrderNp^2 + ( (CurrentIndex-1)*HighOrderNp^2+1:CurrentIndex*HighOrderNp^2 )) = J_index;
    VK(size( elements,2 )*HighOrderNp^2 + ( (CurrentIndex-1)*HighOrderNp^2+1:CurrentIndex*HighOrderNp^2 )) = stima(:);
   end
end

K = sparse(IK,JK,VK,size(coordinates,1),size(coordinates,1));
dirichlet = setdiff(boundary_faces, neumann, 'rows') ; 
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
% Dirichlet conditions
u(unique(dirichlet)) = p(coordinates(unique(dirichlet),:));
b = b - K * u;
% Computation of the solution
u(FreeNodes) = K(FreeNodes,FreeNodes) \ b(FreeNodes);

% 3d slices
[XX,YY,ZZ] = meshgrid( linspace(xmin,xmax,Cell_N+1), linspace(ymin,ymax,Cell_M+1), linspace(zmin,zmax,Cell_L+1) );
VV = reshape(u,Cell_N+1,Cell_M+1,Cell_L+1);
VV = permute(VV,[2,1,3]);
xslice = [];   
yslice = [];
zslice = [0.1,0.3,0.5,0.7];
figure;slice(XX,YY,ZZ,VV,xslice,yslice,zslice);
colormap('jet');
shading interp
xlim([xmin,xmax]); ylim([ymin,ymax]); zlim([zmin,zmax]);

save('SingleFractureData')
function km = km(xyz)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
km = ones(size(x)); 
end
function f = f(xyz)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
f = zeros(size(x));
end
function p = p(xyz)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
p = 1-x;
end
function g = g(xyz,n)
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
g = zeros(size(x));
end









end