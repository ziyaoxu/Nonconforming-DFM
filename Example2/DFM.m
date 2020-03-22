function DFM
% time : 2019.6.27 - 2019.6.27
% author : xuziyao
% this program is code for Continuous Finite Element Method for
% Discrete Fracture Model on triangular meshes.
% PDE model: 
% -div(D*grad(p)) = f, where D contains the fracture
% Neumann & Dirichlet boundary condition 
% upper and lower boundary: Neumann
% left and right boundary: Dirichlet
% Basis function is linear Lagrange basis on the triangular elements.
% FINISHED
clc,clear
format long
% set the geometry and triangulation parameters:
Radius = 1;
h0 = 0.1; Geps = h0 * 1e-5 ;
% ParametersofFracture: x_c,y_c,length,theta,width,permeability,xa,ya,xb,yb,# of elements passed by this fracture. 
NumberofFractures = 3; 
ParametersofFracture = zeros(NumberofFractures, 11); 
ParametersofFracture(:,1:4) = [ 0, 0, Radius, 0;... % position of fractures
                                0, 0, Radius, pi/3;...
                                0, 0, Radius, -pi/3];
ParametersofFracture(:,5) = 1e-4; % width of fractures
ParametersofFracture(:,6)= [1e4;2e4;3e4] ; % permeability of fractures
% plot the triangulation and fractures
ParametersofFracture(:,7) = ParametersofFracture(:,1) - 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,8) = ParametersofFracture(:,2) - 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
ParametersofFracture(:,9) = ParametersofFracture(:,1) + 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,10) = ParametersofFracture(:,2) + 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
NumberofFractureNodes = sum( int64(ParametersofFracture(:,3)/h0+1) );
CoordofFractureNodes = zeros(NumberofFractureNodes,2);
for k = 1 : NumberofFractures
    CoordofFractureNodes(sum(int64(ParametersofFracture(1:(k-1),3)/h0+1))+1:sum(int64(ParametersofFracture(1:k,3)/h0+1)),:) = ...
        [ linspace(ParametersofFracture(k,7),ParametersofFracture(k,9),int64(ParametersofFracture(k,3)/h0+1));...
          linspace(ParametersofFracture(k,8),ParametersofFracture(k,10),int64(ParametersofFracture(k,3)/h0+1)) ]';
end
pfix = CoordofFractureNodes; % Fixed node positions (NFIXx2)
% Uniform Mesh on Unit Circle
fd=@(p) sqrt(sum(p.^2,2))-Radius;
[coordinates,elements3]=distmesh2d(fd,@huniform,h0,[-Radius,-Radius;Radius,Radius],pfix); elements3 = elements3';
EtoEmap = triangulation_neighbor_triangles ( elements3 ); 
[ Tarea , Elenth , Jacobimat ] = ComputeGeometryData( elements3 ,coordinates );%Jacobimat=[Dr/Dx,Dr/Dy ; Ds/Dx,Ds/Dy] 
boundary_edge = boundedges(coordinates,elements3') ;
%neumann = boundary_edge ; dirichlet = setdiff(boundary_edge, neumann, 'rows') ; 
neumann = zeros(size(boundary_edge)) ; neumann_count = 0;  
EtoFmap = zeros(size(elements3,1)*NumberofFractures,size(elements3,2)); 
for k = 1 : NumberofFractures
for jj = 1 : size(elements3,2)
    for ii = 1 : 3
        vertex1 = mod(ii,3)+1 ; 
        vertex2 = mod(ii+1,3)+1 ;
        if all(ismember([elements3(vertex1,jj);elements3(vertex2,jj)],...
                find( donesegment(coordinates,[ParametersofFracture(k,7:8);ParametersofFracture(k,9:10)])<Geps )))
            EtoFmap((k-1)*3+ii,jj) = 1 ;
        end
    end
end
end
figure;pdeplot(coordinates',elements3,'EdgeColor','black');
hold on; plot(ParametersofFracture(:,[7,9])',ParametersofFracture(:,[8,10])',...
   'r-','LineWidth',1); hold off; 
axis equal; axis off;ax=axis;axis(ax*1.001); view([0,90]);drawnow;
% (2) compute basis data
HighOrderDegree = 1 ; % polynomial's order in each element 
Gauss_num = 6 ; % the number of Gauss quadrature point in [-1,1]
Lobatto_num = 6 ; % the number of Lobatto quadrature point in [-1,1]
[Gauss_x,Gauss_w] = JacobiGQ(0,0,Gauss_num-1); 
Gauss_x = 0 + ( Gauss_x - (-1) ) / 2 ; % Gauss_x now is in [0,1].
Gauss_w = Gauss_w / 2 ;% sum(Gauss_w) = 1
[ quad_r , quad_s , quad_w ] = quad_rule1 ( 10 ) ;% Hammer quadrature rule
quad_num = size(quad_w,1) ; % the number of quadrature point in reference triangular 
basisP = zeros( 3 , quad_num ) ; % basis function's value in quadrature point
GradrbasisP = zeros( 3 , quad_num ) ; % gradient's value on quadrature point
GradsbasisP = zeros( 3 , quad_num ) ; % gradient's value on quadrature point
% compute the basis function and it's gradient's value on quadrature point:
basisP(1,:) = quad_r' ;
basisP(2,:) = quad_s' ;
basisP(3,:) = (1-quad_r-quad_s)' ;
GradrbasisP(1,:) = 1 ; 
GradrbasisP(2,:) = 0 ; 
GradrbasisP(3,:) = -1 ; 
GradsbasisP(1,:) = 0;
GradsbasisP(2,:) = 1;
GradsbasisP(3,:) = -1;
% ±ß½ç»ý·Ö¾ØÕó
quad_rs = zeros( Gauss_num,2,3 ) ;
% coordinate(rs) of Gauss points of edge A1A2 in reference triangular.
quad_rs(:,:,1) = ones(size(Gauss_x))*[0,1] + Gauss_x*([0,0]-[0,1]) ; 
% coordinate(rs) of Gauss points of edge A2A3 in reference triangular.
quad_rs(:,:,2) = ones(size(Gauss_x))*[0,0]+ Gauss_x*([1,0]-[0,0]) ; 
% coordinate(rs) of Gauss points of edge A3A1 in reference triangular.
quad_rs(:,:,3) = ones(size(Gauss_x))*[1,0] + Gauss_x*([0,1]-[1,0])  ; 
boundaryP = zeros(3,Gauss_num,3);
GradrboundaryP= zeros(3,Gauss_num,3);
GradsboundaryP= zeros(3,Gauss_num,3);
for ii = 1 : 3
    boundaryP(1,:,ii) = quad_rs(:,1,ii)';  
    boundaryP(2,:,ii) = quad_rs(:,2,ii)';  
    boundaryP(3,:,ii) = (1-quad_rs(:,1,ii)-quad_rs(:,2,ii))';   
    GradrboundaryP(1,:,ii) = 1 * ones(1,Gauss_num); 
    GradrboundaryP(2,:,ii) = 0 * ones(1,Gauss_num); 
    GradrboundaryP(3,:,ii) = -1 * ones(1,Gauss_num); 
    GradsboundaryP(1,:,ii) = 0 * ones(1,Gauss_num); 
    GradsboundaryP(2,:,ii) = 1 * ones(1,Gauss_num); 
    GradsboundaryP(3,:,ii) = -1 * ones(1,Gauss_num);  
end
% define stiffness metrix and right hand side vector
K = sparse(size(coordinates,1),size(coordinates,1)); % stiffness matrix
b = sparse(size(coordinates,1),1);% right hand side of the scheme
u_DFM = sparse(size(coordinates,1),1);% solution vector
for jj = 1 : size( elements3,2 ) % assembly K & b 
xy = quad_r * coordinates(elements3(1,jj),:)... 
    + quad_s * coordinates(elements3(2,jj),:)... 
    +(1-quad_r-quad_s)*coordinates(elements3(3,jj),:) ; 
% assembly stiffness matrix
K(elements3(:,jj),elements3(:,jj)) = K(elements3(:,jj),elements3(:,jj)) + ...
         Tarea(jj)* ( ... 
    (GradrbasisP*Jacobimat(1,1,jj)+GradsbasisP*Jacobimat(2,1,jj))...
        *diag(km(xy).*quad_w)...
        *(GradrbasisP*Jacobimat(1,1,jj)+GradsbasisP*Jacobimat(2,1,jj))' ...
     + ...
    (GradrbasisP*Jacobimat(1,2,jj)+GradsbasisP*Jacobimat(2,2,jj))...
        *diag(km(xy).*quad_w)...
        *(GradrbasisP*Jacobimat(1,2,jj)+GradsbasisP*Jacobimat(2,2,jj))' ) ; 
% assembly right hand side matrix 
b(elements3(:,jj)) = b(elements3(:,jj)) + Tarea(jj)*basisP*(f(xy).*quad_w) ;
% Neumann conditions
for ii = 1 : 3  % the ii-th edge of jj-th element
    vertex1 = mod(ii,3)+1 ; 
    vertex2 = mod(ii+1,3)+1 ;
    edge_vec =  coordinates(elements3(vertex2,jj),:)' - ...
        coordinates(elements3(vertex1,jj),:)' ;
    outer_vec = - [0,-1;1,0]*(edge_vec/norm(edge_vec));
    E = EtoEmap(ii,jj) ; % the index of neighborhood of element jj
    if E == -1
        if abs(outer_vec(2))<-1/2 % up and down boundary, Neumann boundary condition, B1,B2 (p) take inner value
        neumann_count = neumann_count + 1; neumann(neumann_count,:) = elements3([vertex1;vertex2],jj) ;   
        Boundary_xy = (1-Gauss_x)*coordinates(elements3(vertex1,jj),:) + Gauss_x*coordinates(elements3(vertex2,jj),:);
        b(elements3([vertex1;vertex2],jj)) = b(elements3([vertex1;vertex2],jj)) + ...
            Elenth(ii,jj)*boundaryP([vertex1;vertex2],:,ii)*(g(Boundary_xy).*Gauss_w);
        end
    end
end
end
% Fracture Term Assembly
GmulG = [1,-1;-1,1] ; % Stiffness matrix on referents line elements
for k = 1 : NumberofFractures
for jj = 1 : size(elements3,2) 
    if any(EtoFmap((k-1)*3+1:k*3,jj))
    for ii = 1 : 3
        vertex1 = mod(ii,3)+1 ; 
        vertex2 = mod(ii+1,3)+1 ;
        if EtoFmap((k-1)*3+ii,jj)~=0
            K(elements3([vertex1;vertex2],jj),elements3([vertex1;vertex2],jj)) = ...
                K(elements3([vertex1;vertex2],jj),elements3([vertex1;vertex2],jj)) + ...
                ParametersofFracture(k,5)/2*ParametersofFracture(k,6)* GmulG  ...
                / Elenth(ii,jj);
                
        end
    end
    end
end
end
dirichlet = setdiff(boundary_edge, neumann, 'rows') ; 
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
% Dirichlet conditions
u_DFM(unique(dirichlet)) = p(coordinates(unique(dirichlet),:));
b = b - K * u_DFM;
% Computation of the solution
u_DFM(FreeNodes) = K(FreeNodes,FreeNodes) \ b(FreeNodes);
figure; showsurf( coordinates , elements3 , u_DFM );xlim([-Radius,Radius]); ylim([-Radius Radius])
colormap('jet'); cb=colorbar;cb.Position = cb.Position.*[1.06,1,1,1];
hold on; plot3(ParametersofFracture(:,[7,9])',ParametersofFracture(:,[8,10])',(max(u_DFM)+1)*ones(size(ParametersofFracture(:,[8,10])')),...
   'w-','LineWidth',0.05); hold off; axis equal; axis off;ax=axis;axis(ax*1.001); view([0,90]);

save('DFMData','Radius','NumberofFractures','ParametersofFracture','coordinates','elements3','u_DFM')
% -------------------------------------------------------------------------

end

