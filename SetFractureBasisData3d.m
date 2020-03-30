function [ParametersofFracture,FracturesPath,FracturesArea,basisPonFracture, ...
    GradxbasisPonFracture,GradybasisPonFracture,GradzbasisPonFracture,GradsigmabasisPonFracture] =...
    SetFractureBasisData3d(coordinates,elements,ParametersofFracture,quad_x)
quad_num = size(quad_x,1);
% Global geometric information:
xmax = max(coordinates(:,1)); xmin = min(coordinates(:,1));
ymax = max(coordinates(:,2)); ymin = min(coordinates(:,2));
zmax = max(coordinates(:,3)); zmin = min(coordinates(:,3));
% x direction is divided into N parts ,y direction is divided into M parts, z direction is divided into L parts
Cell_N = (elements(4,1)-elements(1,1)-1) ; 
Cell_M = (elements(5,1)-elements(1,1)) / (Cell_N+1) - 1 ;
Cell_L = size(elements,2) / (Cell_N*Cell_M) ; 
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ; hz = ( zmax - zmin ) / Cell_L ;
% set fractures data:
NumberofFractures = size(ParametersofFracture,1);
% ParametersofFracture: width,permeability,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,sigma1,sigma2,sigma3,...
%                       # of elements passed by this fracture. 
for k = 1 : NumberofFractures
    coor_vertex1 = ParametersofFracture(k,3:5)';
    coor_vertex2 = ParametersofFracture(k,6:8)';
    coor_vertex3 = ParametersofFracture(k,9:11)';
    coor_vertex4 = ParametersofFracture(k,12:14)';
    sigma = cross((coor_vertex1-coor_vertex2) , (coor_vertex3-coor_vertex2));
    sigma = sigma/norm(sigma,2);
    ParametersofFracture(k,15:17) = sigma;
    if abs(sigma(3))<0.1 % perpendicular to XOY
        [Ix,Jy,~,~,~] = SetFractureBasisData2d(xmin,xmax,ymin,ymax,Cell_N,Cell_M,hx,hy,...
            coor_vertex1(1),coor_vertex1(2),coor_vertex2(1),coor_vertex2(2),quad_x);
        Kz = [ceil((coor_vertex2(3)-zmin)/hz):ceil((coor_vertex3(3)-zmin)/hz)]';
        ParametersofFracture(k,18) = size(Kz,1)*size(Ix,1);
    elseif abs(sigma(1))<0.1 % perpendicular to YOZ
        [Jy,Kz,sy,tz,FracturesLength] = SetFractureBasisData2d(ymin,ymax,zmin,zmax,Cell_M,Cell_L,hy,hz,...
            coor_vertex1(2),coor_vertex1(3),coor_vertex2(2),coor_vertex2(3),quad_x);
        Ix = [ceil((coor_vertex2(1)-xmin)/hx):ceil((coor_vertex3(1)-xmin)/hx)]';
        ParametersofFracture(k,18) = size(Ix,1)*size(Jy,1);
    elseif abs(sigma(2))<0.1 % perpendicular to ZOX
        [Kz,Ix,tz,rx,FracturesLength] = SetFractureBasisData2d(zmin,zmax,xmin,xmax,Cell_L,Cell_N,hz,hx,...
            coor_vertex1(3),coor_vertex1(1),coor_vertex2(3),coor_vertex2(1),quad_x);
        Jy = [ceil((coor_vertex2(2)-ymin)/hy):ceil((coor_vertex3(2)-ymin)/hy)]';
        ParametersofFracture(k,18) = size(Jy,1)*size(Ix,1);
    else
        error('wrong')
    end
end

FracturesPath = zeros(sum(ParametersofFracture(:,18)),1);
FracturesArea = zeros(sum(ParametersofFracture(:,18)),1);
basisPonFracture = zeros(8*sum(ParametersofFracture(:,18)),quad_num^2);
GradxbasisPonFracture = zeros(8*sum(ParametersofFracture(:,18)),quad_num^2);
GradybasisPonFracture = zeros(8*sum(ParametersofFracture(:,18)),quad_num^2);
GradzbasisPonFracture = zeros(8*sum(ParametersofFracture(:,18)),quad_num^2);
GradsigmabasisPonFracture = zeros(8*sum(ParametersofFracture(:,18)),quad_num^2);

for k = 1 : NumberofFractures
    coor_vertex1 = ParametersofFracture(k,3:5)';
    coor_vertex2 = ParametersofFracture(k,6:8)';
    coor_vertex3 = ParametersofFracture(k,9:11)';
    coor_vertex4 = ParametersofFracture(k,12:14)';
    sigma = cross((coor_vertex1-coor_vertex2) , (coor_vertex3-coor_vertex2));
    sigma = sigma/norm(sigma,2);
    ParametersofFracture(k,15:17) = sigma;
    if abs(sigma(3))<0.1 % perpendicular to XOY
        [Ix,Jy,rx,sy,FracturesLength] = SetFractureBasisData2d(xmin,xmax,ymin,ymax,Cell_N,Cell_M,hx,hy,...
            coor_vertex1(1),coor_vertex1(2),coor_vertex2(1),coor_vertex2(2),quad_x);
        Kz = [ceil((coor_vertex2(3)-zmin)/hz):ceil((coor_vertex3(3)-zmin)/hz)]';
        I = repmat(Ix,size(Kz,1),1);
        J = repmat(Jy,size(Kz,1),1);
        K = kron(Kz,ones(size(Ix,1),1));
        r = repmat(rx,size(Kz,1),quad_num);
        s = repmat(sy,size(Kz,1),quad_num);
        if size(Kz,1)==1
            tz = ((coor_vertex2(3)-zmin)/hz-Kz(1)+1) + (coor_vertex3(3)-coor_vertex2(3))/hz * kron(quad_x',ones(1,quad_num));
            t = repmat(tz,size(Ix,1),1);
            FracturesA = repmat(FracturesLength,size(Kz,1),1);
            FracturesA = FracturesA *(coor_vertex3(3)-coor_vertex2(3));
        elseif size(Kz,1)==2
            tz1 = ((coor_vertex2(3)-zmin)/hz-Kz(1)+1) + (Kz(1)*hz-coor_vertex2(3))/hz * kron(quad_x',ones(1,quad_num));
            tz2 = 0 + (coor_vertex3(3)-(Kz(end)-1)*hz)/hz * kron(quad_x',ones(1,quad_num));
            t = [repmat(tz1,size(Ix,1),1);repmat(tz2,size(Ix,1),1)];
            FracturesA = repmat(FracturesLength,size(Kz,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Kz(1)*hz-coor_vertex2(3));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(3)-(Kz(end)-1)*hz);
        else
            tz1 = ((coor_vertex2(3)-zmin)/hz-Kz(1)+1) + (Kz(1)*hz-coor_vertex2(3))/hz * kron(quad_x',ones(1,quad_num));
            tz2 = 0 + (coor_vertex3(3)-(Kz(end)-1)*hz)/hz * kron(quad_x',ones(1,quad_num));
            t = [repmat(tz1,size(Ix,1),1);repmat(kron(quad_x',ones(1,quad_num)),size(Ix,1)*(size(Kz,1)-2),1);repmat(tz2,size(Ix,1),1)];
            FracturesA = repmat(FracturesLength,size(Kz,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Kz(1)*hz-coor_vertex2(3));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(3)-(Kz(end)-1)*hz);
            FracturesA(size(Ix,1)+1:end-size(Ix,1),: ) = FracturesA(size(Ix,1)+1:end-size(Ix,1),: ) *hz;
        end
    elseif abs(sigma(1))<0.1 % perpendicular to YOZ
        [Jy,Kz,sy,tz,FracturesLength] = SetFractureBasisData2d(ymin,ymax,zmin,zmax,Cell_M,Cell_L,hy,hz,...
            coor_vertex1(2),coor_vertex1(3),coor_vertex2(2),coor_vertex2(3),quad_x);
        Ix = [ceil((coor_vertex2(1)-xmin)/hx):ceil((coor_vertex3(1)-xmin)/hx)]';
        J = repmat(Jy,size(Ix,1),1);
        K = repmat(Kz,size(Ix,1),1);
        I = kron(Ix,ones(size(Jy,1),1));
        s = repmat(sy,size(Ix,1),quad_num);
        t = repmat(tz,size(Ix,1),quad_num);
        if size(Ix,1)==1
            rx = ((coor_vertex2(1)-xmin)/hx-Ix(1)+1) + (coor_vertex3(1)-coor_vertex2(1))/hx * kron(quad_x',ones(1,quad_num));
            r = repmat(rx,size(Jy,1),1);
            FracturesA = repmat(FracturesLength,size(Ix,1),1);
            FracturesA = FracturesA *(coor_vertex3(1)-coor_vertex2(1));
        elseif size(Ix,1)==2
            rx1 = ((coor_vertex2(1)-xmin)/hx-Ix(1)+1) + (Ix(1)*hx-coor_vertex2(1))/hx * kron(quad_x',ones(1,quad_num));
            rx2 = 0 + (coor_vertex3(1)-(Ix(end)-1)*hx)/hx * kron(quad_x',ones(1,quad_num));
            r = [repmat(rx1,size(Jy,1),1);repmat(rx2,size(Jy,1),1)];
            FracturesA = repmat(FracturesLength,size(Ix,1),1);
            FracturesA(1:size(Jy,1),: ) = FracturesA(1:size(Jy,1),: ) *(Ix(1)*hx-coor_vertex2(1));
            FracturesA(end-size(Jy,1)+1:end,: ) = FracturesA(end-size(Jy,1)+1:end,: ) *(coor_vertex3(1)-(Ix(end)-1)*hx);
        else
            rx1 = ((coor_vertex2(1)-xmin)/hx-Ix(1)+1) + (Ix(1)*hx-coor_vertex2(1))/hx * kron(quad_x',ones(1,quad_num));
            rx2 = 0 + (coor_vertex3(1)-(Ix(end)-1)*hx)/hx * kron(quad_x',ones(1,quad_num));
            r = [repmat(rx1,size(Jy,1),1);repmat(kron(quad_x',ones(1,quad_num)),size(Jy,1)*(size(Ix,1)-2),1);repmat(rx2,size(Jy,1),1)];
            FracturesA = repmat(FracturesLength,size(Ix,1),1);
            FracturesA(1:size(Jy,1),: ) = FracturesA(1:size(Jy,1),: ) *(Ix(1)*hx-coor_vertex2(1));
            FracturesA(end-size(Jy,1)+1:end,: ) = FracturesA(end-size(Jy,1)+1:end,: ) *(coor_vertex3(1)-(Ix(end)-1)*hx);
            FracturesA(size(Jy,1)+1:end-size(Jy,1),: ) = FracturesA(size(Jy,1)+1:end-size(Jy,1),: ) *hx;
        end
    elseif abs(sigma(2))<0.1 % perpendicular to ZOX
        [Kz,Ix,tz,rx,FracturesLength] = SetFractureBasisData2d(zmin,zmax,xmin,xmax,Cell_L,Cell_N,hz,hx,...
            coor_vertex1(3),coor_vertex1(1),coor_vertex2(3),coor_vertex2(1),quad_x);
        Jy = [ceil((coor_vertex2(2)-ymin)/hy):ceil((coor_vertex3(2)-ymin)/hy)]';
        K = repmat(Kz,size(Jy,1),1);
        I = repmat(Ix,size(Jy,1),1);
        J = kron(Jy,ones(size(Kz,1),1));
        t = repmat(tz,size(Jy,1),quad_num);
        r = repmat(rx,size(Jy,1),quad_num);
        if size(Jy,1)==1
            sy = ((coor_vertex2(2)-ymin)/hy-Jy(1)+1) + (coor_vertex3(2)-coor_vertex2(2))/hy * kron(quad_x',ones(1,quad_num));
            s = repmat(sy,size(Ix,1),1);
            FracturesA = repmat(FracturesLength,size(Jy,1),1);
            FracturesA = FracturesA *(coor_vertex3(2)-coor_vertex2(2));
        elseif size(Jy,1)==2
            sy1 = ((coor_vertex2(2)-ymin)/hy-Jy(1)+1) + (Jy(1)*hy-coor_vertex2(2))/hy * kron(quad_x',ones(1,quad_num));
            sy2 = 0 + (coor_vertex3(2)-(Jy(end)-1)*hy)/hy * kron(quad_x',ones(1,quad_num));
            s = [repmat(sy1,size(Ix,1),1);repmat(sy2,size(Ix,1),1)];
            FracturesA = repmat(FracturesLength,size(Jy,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Jy(1)*hy-coor_vertex2(2));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(2)-(Jy(end)-1)*hy);
        else
            sy1 = ((coor_vertex2(2)-ymin)/hy-Jy(1)+1) + (Jy(1)*hy-coor_vertex2(2))/hy * kron(quad_x',ones(1,quad_num));
            sy2 = 0 + (coor_vertex3(2)-(Jy(end)-1)*hy)/hy * kron(quad_x',ones(1,quad_num));
            s = [repmat(sy1,size(Ix,1),1);repmat(kron(quad_x',ones(1,quad_num)),size(Ix,1)*(size(Jy,1)-2),1);repmat(sy2,size(Ix,1),1)];
            FracturesA = repmat(FracturesLength,size(Jy,1),1);
            FracturesA(1:size(Ix,1),: ) = FracturesA(1:size(Ix,1),: ) *(Jy(1)*hy-coor_vertex2(2));
            FracturesA(end-size(Ix,1)+1:end,: ) = FracturesA(end-size(Ix,1)+1:end,: ) *(coor_vertex3(2)-(Jy(end)-1)*hy);
            FracturesA(size(Ix,1)+1:end-size(Ix,1),: ) = FracturesA(size(Ix,1)+1:end-size(Ix,1),: ) *hy;
        end
    else
        error('wrong')
    end
FracturesPath(sum(ParametersofFracture(1:k-1,18))+1:sum(ParametersofFracture(1:k,18)),1) = ...
    (K-1)*Cell_N*Cell_M + (J-1)*Cell_N + I;
FracturesArea(sum(ParametersofFracture(1:k-1,18))+1:sum(ParametersofFracture(1:k,18)),1) = ...
    FracturesA;
% 0 ---
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*(1-s).*(1-t);
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[2:8:8*ParametersofFracture(k,18)],:) = ...
    r.*(1-s).*(1-t);
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[3:8:8*ParametersofFracture(k,18)],:) = ...
    r.*s.*(1-t);
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[4:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*s.*(1-t);
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[5:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*(1-s).*t;
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[6:8:8*ParametersofFracture(k,18)],:) = ...
    r.*(1-s).*t;
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[7:8:8*ParametersofFracture(k,18)],:) = ...
    r.*s.*t;
basisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[8:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*s.*t;
% 1 ---
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-s).*(1-t)/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[2:8:8*ParametersofFracture(k,18)],:) = ...
     (1-s).*(1-t)/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[3:8:8*ParametersofFracture(k,18)],:) = ...
    s.*(1-t)/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[4:8:8*ParametersofFracture(k,18)],:) = ...
    -s.*(1-t)/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[5:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-s).*t/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[6:8:8*ParametersofFracture(k,18)],:) = ...
    (1-s).*t/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[7:8:8*ParametersofFracture(k,18)],:) = ...
    s.*t/hx;
GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[8:8:8*ParametersofFracture(k,18)],:) = ...
    -s.*t/hx;
% 2 ---
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-r).*(1-t)/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[2:8:8*ParametersofFracture(k,18)],:) = ...
    -r.*(1-t)/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[3:8:8*ParametersofFracture(k,18)],:) = ...
    r.*(1-t)/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[4:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*(1-t)/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[5:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-r).*t/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[6:8:8*ParametersofFracture(k,18)],:) = ...
    -r.*t/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[7:8:8*ParametersofFracture(k,18)],:) = ...
    r.*t/hy;
GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[8:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*t/hy;
% 3 ---
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-r).*(1-s)/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[2:8:8*ParametersofFracture(k,18)],:) = ...
    -r.*(1-s)/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[3:8:8*ParametersofFracture(k,18)],:) = ...
    -r.*s/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[4:8:8*ParametersofFracture(k,18)],:) = ...
    -(1-r).*s/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[5:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*(1-s)/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[6:8:8*ParametersofFracture(k,18)],:) = ...
    r.*(1-s)/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[7:8:8*ParametersofFracture(k,18)],:) = ...
    r.*s/hz;
GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[8:8:8*ParametersofFracture(k,18)],:) = ...
    (1-r).*s/hz;
% 4 
GradsigmabasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8*ParametersofFracture(k,18)],:) = ...
    sigma(1)*GradxbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8*ParametersofFracture(k,18)],:)...
   +sigma(2)*GradybasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8*ParametersofFracture(k,18)],:)...
   +sigma(3)*GradzbasisPonFracture(8*sum(ParametersofFracture(1:k-1,18))+[1:8*ParametersofFracture(k,18)],:) ;
end


function [Ix,Jy,r,s,FracturesLength] = SetFractureBasisData2d(xmin,xmax,ymin,ymax,Cell_N,Cell_M,hx,hy,xa,ya,xb,yb,quad_x)
[coordinates2,elements4,EtoEmap] = RectangularMesh( xmin,xmax,ymin,ymax,Cell_N,Cell_M );        
StartEndElements = zeros(1,2);
StartEndElements(1) = ( ceil((ya-ymin)/hy)-1 )*Cell_N + ceil((xa-xmin)/hx);
StartEndElements(2) = ( ceil((yb-ymin)/hy)-1 )*Cell_N + ceil((xb-xmin)/hx);
CurrentElement = StartEndElements(1); 
NumElePassed = 1; 
pp = 0;
while ( CurrentElement ~= StartEndElements(2) )
[~, ~, ii_int] = polyxpoly(coordinates2(elements4([1;2;3;4;1],CurrentElement),1),coordinates2(elements4([1;2;3;4;1],CurrentElement),2)...
   ,[xa,xb],[ya,yb]); edge_int = ii_int(:,1);
qq = setdiff(edge_int,pp);
NextElement = EtoEmap(qq,CurrentElement);
pp = find(EtoEmap(:,NextElement) == CurrentElement) ; 
CurrentElement = NextElement;
NumElePassed = NumElePassed + 1;
end
Ix = zeros(NumElePassed,1);
Jy = zeros(NumElePassed,1);
r = zeros(NumElePassed,size(quad_x,1));
s = zeros(NumElePassed,size(quad_x,1));
FracturesLength = zeros(NumElePassed,1);
CurrentElement = StartEndElements(1); CurrentIndex = 1;
Ix(CurrentIndex) = mod(CurrentElement-1,Cell_N)+1;
Jy(CurrentIndex) = ceil(CurrentElement/Cell_N);
pp = 0; pplamda = ([xa;ya]-coordinates2(elements4(1,CurrentElement),:)')./[hx;hy];
while ( CurrentElement ~= StartEndElements(2) )
    [x_int, y_int, ii_int] = polyxpoly(coordinates2(elements4([1;2;3;4;1],CurrentElement),1),coordinates2(elements4([1;2;3;4;1],CurrentElement),2)...
       ,[xa,xb],[ya,yb]); edge_int = ii_int(:,1);
    qq = setdiff(edge_int,pp); qqlamda = ([x_int(edge_int==qq);y_int(edge_int==qq)]-coordinates2(elements4(1,CurrentElement),:)')./[hx;hy];
      if length(edge_int) == 1 
            FracturesLength(CurrentIndex) = ...
                sqrt( (x_int-xa)^2 + (y_int-ya)^2 );
      elseif length(edge_int) == 2 
            FracturesLength(CurrentIndex) = sqrt( diff(x_int)^2 + diff(y_int)^2 );
      else
          error('intersection points is not one or two');
      end
    % compute basis data (basisPonFracture,GradnubasisPonFracture) in the following:
    r(CurrentIndex,:) = pplamda(1) + (qqlamda(1)-pplamda(1))*[quad_x] ; 
    s(CurrentIndex,:) = pplamda(2) + (qqlamda(2)-pplamda(2))*[quad_x] ; 
    NextElement = EtoEmap(qq,CurrentElement);
    pp = find(EtoEmap(:,NextElement) == CurrentElement); pplamda = ([x_int(edge_int==qq);y_int(edge_int==qq)]-coordinates2(elements4(1,NextElement),:)')./[hx;hy];
    CurrentElement = NextElement; CurrentIndex = CurrentIndex + 1;
Ix(CurrentIndex) = mod(CurrentElement-1,Cell_N)+1;
Jy(CurrentIndex) = ceil(CurrentElement/Cell_N);
end
[x_int, y_int, ~] = polyxpoly(coordinates2(elements4([1;2;3;4;1],CurrentElement),1),coordinates2(elements4([1;2;3;4;1],CurrentElement),2)...
   ,[xa,xb],[ya,yb]);    
qqlamda = ([xb;yb]-coordinates2(elements4(1,CurrentElement),:)')./[hx;hy];
FracturesLength(CurrentIndex) = ...
       sqrt( (x_int-xb)^2 + (y_int-yb)^2 );
% compute basis data (basisPonFracture,GradnubasisPonFracture) in the following:
r(CurrentIndex,:) = pplamda(1) + (qqlamda(1)-pplamda(1))*[quad_x] ; 
s(CurrentIndex,:) = pplamda(2) + (qqlamda(2)-pplamda(2))*[quad_x] ; 
end

end

