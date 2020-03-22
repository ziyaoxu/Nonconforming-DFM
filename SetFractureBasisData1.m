function [ParametersofFracture,FracturesPath,FracturesLength,basisPonFracture,GradnubasisPonFracture,GradsigmabasisPonFracture] = ...
    SetFractureBasisData0(coordinates,elements4,ParametersofFracture,quad_x)
% time : 2019.3.20 - 2019.3.20
% modified on 2019.6.26, the basis data on tip of fractures is no longer zero. 
% author : xuziyao
% on rectangular meshes
% quad_x is distributed in [0,1] as the quadrature points along the fractures.
% Basis function is the bilinear Lagrange basis on the [0,1]*[0,1] reference square.
% phi_1(r,s) = (1-r)(1-s)
% phi_2(r,s) = r(1-s)
% phi_3(r,s) = rs
% phi_4(r,s) = (1-r)s
% Gradrphi_1 = -(1-s); Gradsphi_1 = -(1-r).
% Gradrphi_2 =   1-s ; Gradsphi_2 =    -r.
% Gradrphi_3 =    s  ; Gradsphi_3 =     r.
% Gradrphi_4 =   -s  ; Gradsphi_4 =   1-r .
ParametersofFracture(:,11) = 0; % the number of elements which is passed by the k-th fracture
EtoEmap = Rectangulation_neighbor_rectangles( elements4 ); 
hx = ( max(coordinates(:,1)) - min(coordinates(:,1)) )/ (elements4(4,1)-elements4(1,1)-1) ; 
hy = ( max(coordinates(:,2)) - min(coordinates(:,2)) )/ (size(elements4,2)/(elements4(4,1)-elements4(1,1)-1)) ; 
NumberofFractures = size(ParametersofFracture,1);
% ParametersofFracture: x_c,y_c,length,theta,width,permeability,xa,ya,xb,yb,
% # of elements which are passed through by this fracture. 
ParametersofFracture(:,7) = ParametersofFracture(:,1) - 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,8) = ParametersofFracture(:,2) - 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
ParametersofFracture(:,9) = ParametersofFracture(:,1) + 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,10) = ParametersofFracture(:,2) + 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
StartEndElements = zeros(NumberofFractures,2);
for jj = 1 : size(elements4,2) % compute the start and end elements of each fracture's path
    isStartElement = inpolygon( ParametersofFracture(:,7),ParametersofFracture(:,8),...
        coordinates(elements4([1;2;3;4;1],jj),1),coordinates(elements4([1;2;3;4;1],jj),2) );
    isEndElement   = inpolygon( ParametersofFracture(:,9),ParametersofFracture(:,10),...
        coordinates(elements4([1;2;3;4;1],jj),1),coordinates(elements4([1;2;3;4;1],jj),2) );
    StartEndElements(isStartElement,1) = jj;
    StartEndElements(isEndElement,2) = jj;
end
for k = 1 : NumberofFractures % compute ParametersofFracture(:,11), the number of elements along each fracture's path
    CurrentElement = StartEndElements(k,1); 
    ParametersofFracture(k,11) = ParametersofFracture(k,11) + 1;
    pp = 0;
    while ( CurrentElement ~= StartEndElements(k,2) )
        [~, ~, ii_int] = polyxpoly(coordinates(elements4([1;2;3;4;1],CurrentElement),1),coordinates(elements4([1;2;3;4;1],CurrentElement),2)...
           ,ParametersofFracture(k,[7,9]),ParametersofFracture(k,[8,10])); edge_int = ii_int(:,1);
        qq = setdiff(edge_int,pp);
        NextElement = EtoEmap(qq,CurrentElement);
        pp = find(EtoEmap(:,NextElement) == CurrentElement) ; 
        CurrentElement = NextElement;
        ParametersofFracture(k,11) = ParametersofFracture(k,11) + 1;
    end
end
FracturesPath = zeros(sum(ParametersofFracture(:,11)),1);
FracturesLength = zeros(sum(ParametersofFracture(:,11)),1);
basisPonFracture = zeros( 4*length(FracturesPath) , length(quad_x) + 2 );
GradnubasisPonFracture = zeros( 4*length(FracturesPath) , length(quad_x) + 2 );
GradsigmabasisPonFracture = zeros( 4*length(FracturesPath) , length(quad_x) + 2 );
CurrentIndex = 0;
for k = 1 : NumberofFractures
% compute the basis data for each element which is passed through by the
% k-th fracture:
    CurrentElement = StartEndElements(k,1); CurrentIndex = CurrentIndex + 1;
    FracturesPath(CurrentIndex) = CurrentElement;
    pp = 0; pplamda = ([ParametersofFracture(k,7);ParametersofFracture(k,8)]-coordinates(elements4(1,CurrentElement),:)')./[hx;hy];
    while ( CurrentElement ~= StartEndElements(k,2) )
        [x_int, y_int, ii_int] = polyxpoly(coordinates(elements4([1;2;3;4;1],CurrentElement),1),coordinates(elements4([1;2;3;4;1],CurrentElement),2)...
           ,ParametersofFracture(k,[7,9]),ParametersofFracture(k,[8,10])); edge_int = ii_int(:,1);
        qq = setdiff(edge_int,pp); qqlamda = ([x_int(edge_int==qq);y_int(edge_int==qq)]-coordinates(elements4(1,CurrentElement),:)')./[hx;hy];
          if length(edge_int) == 1 
                FracturesLength(CurrentIndex) = ...
                    sqrt( (x_int-ParametersofFracture(k,7))^2 + (y_int-ParametersofFracture(k,8))^2 );
          elseif length(edge_int) == 2 
                FracturesLength(CurrentIndex) = sqrt( diff(x_int)^2 + diff(y_int)^2 );
          else
              error('intersection points is not one or two');
          end
        % compute basis data (basisPonFracture,GradnubasisPonFracture) in the following:
        quad_lamda1 = pplamda(1) + (qqlamda(1)-pplamda(1))*[quad_x;0;1] ; 
        quad_lamda2 = pplamda(2) + (qqlamda(2)-pplamda(2))*[quad_x;0;1] ; 
        basisPonFracture( (CurrentIndex-1)*4 + 1 , :) = (1-quad_lamda1).*(1-quad_lamda2);
        basisPonFracture( (CurrentIndex-1)*4 + 2 , :) = quad_lamda1.*(1-quad_lamda2);
        basisPonFracture( (CurrentIndex-1)*4 + 3 , :) = quad_lamda1.*quad_lamda2;
        basisPonFracture( (CurrentIndex-1)*4 + 4 , :) = (1-quad_lamda1).*quad_lamda2;
        GradnubasisPonFracture( (CurrentIndex-1)*4 + 1 , :) = ...
            [(-(1-quad_lamda2))/hx,(-(1-quad_lamda1))/hy]*[cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))];
        GradnubasisPonFracture( (CurrentIndex-1)*4 + 2 , :) = ...
            [( (1-quad_lamda2))/hx,(   -quad_lamda1)/hy]*[cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))];
        GradnubasisPonFracture( (CurrentIndex-1)*4 + 3 , :) = ...
            [(    quad_lamda2 ) /hx,(    quad_lamda1)/hy]*[cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))];
        GradnubasisPonFracture( (CurrentIndex-1)*4 + 4 , :) = ...
            [(   -quad_lamda2 )/hx,( (1-quad_lamda1))/hy]*[cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))];
        GradsigmabasisPonFracture( (CurrentIndex-1)*4 + 1 , :) = ...
            [(-(1-quad_lamda2))/hx,(-(1-quad_lamda1))/hy]*[cos(ParametersofFracture(k,4)+pi/2);sin(ParametersofFracture(k,4)+pi/2)];
        GradsigmabasisPonFracture( (CurrentIndex-1)*4 + 2 , :) = ...
            [( (1-quad_lamda2))/hx,(   -quad_lamda1)/hy]*[cos(ParametersofFracture(k,4)+pi/2);sin(ParametersofFracture(k,4)+pi/2)];
        GradsigmabasisPonFracture( (CurrentIndex-1)*4 + 3 , :) = ...
            [(    quad_lamda2 ) /hx,(    quad_lamda1)/hy]*[cos(ParametersofFracture(k,4)+pi/2);sin(ParametersofFracture(k,4)+pi/2)];
        GradsigmabasisPonFracture( (CurrentIndex-1)*4 + 4 , :) = ...
            [(   -quad_lamda2 )/hx,( (1-quad_lamda1))/hy]*[cos(ParametersofFracture(k,4)+pi/2);sin(ParametersofFracture(k,4)+pi/2)];
        NextElement = EtoEmap(qq,CurrentElement);
        pp = find(EtoEmap(:,NextElement) == CurrentElement); pplamda = ([x_int(edge_int==qq);y_int(edge_int==qq)]-coordinates(elements4(1,NextElement),:)')./[hx;hy];
        CurrentElement = NextElement; CurrentIndex = CurrentIndex + 1;
        FracturesPath(CurrentIndex) = CurrentElement;
    end
    [x_int, y_int, ~] = polyxpoly(coordinates(elements4([1;2;3;4;1],CurrentElement),1),coordinates(elements4([1;2;3;4;1],CurrentElement),2)...
       ,ParametersofFracture(k,[7,9]),ParametersofFracture(k,[8,10]));    
    qqlamda = ([ParametersofFracture(k,9);ParametersofFracture(k,10)]-coordinates(elements4(1,CurrentElement),:)')./[hx;hy];
    FracturesLength(CurrentIndex) = ...
           sqrt( (x_int-ParametersofFracture(k,9))^2 + (y_int-ParametersofFracture(k,10))^2 );
        % compute basis data (basisPonFracture,GradnubasisPonFracture) in the following:
        quad_lamda1 = pplamda(1) + (qqlamda(1)-pplamda(1))*[quad_x;0;1] ; 
        quad_lamda2 = pplamda(2) + (qqlamda(2)-pplamda(2))*[quad_x;0;1] ; 
        basisPonFracture( (CurrentIndex-1)*4 + 1 , :) = (1-quad_lamda1).*(1-quad_lamda2);
        basisPonFracture( (CurrentIndex-1)*4 + 2 , :) = quad_lamda1.*(1-quad_lamda2);
        basisPonFracture( (CurrentIndex-1)*4 + 3 , :) = quad_lamda1.*quad_lamda2;
        basisPonFracture( (CurrentIndex-1)*4 + 4 , :) = (1-quad_lamda1).*quad_lamda2;
        GradnubasisPonFracture( (CurrentIndex-1)*4 + 1 , :) = ...
            [(-(1-quad_lamda2))/hx,(-(1-quad_lamda1))/hy]*[cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))];
        GradnubasisPonFracture( (CurrentIndex-1)*4 + 2 , :) = ...
            [( (1-quad_lamda2))/hx,(   -quad_lamda1)/hy]*[cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))];
        GradnubasisPonFracture( (CurrentIndex-1)*4 + 3 , :) = ...
            [(    quad_lamda2 ) /hx,(    quad_lamda1)/hy]*[cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))];
        GradnubasisPonFracture( (CurrentIndex-1)*4 + 4 , :) = ...
            [(   -quad_lamda2 )/hx,( (1-quad_lamda1))/hy]*[cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))];
        GradsigmabasisPonFracture( (CurrentIndex-1)*4 + 1 , :) = ...
            [(-(1-quad_lamda2))/hx,(-(1-quad_lamda1))/hy]*[cos(ParametersofFracture(k,4)+pi/2);sin(ParametersofFracture(k,4)+pi/2)];
        GradsigmabasisPonFracture( (CurrentIndex-1)*4 + 2 , :) = ...
            [( (1-quad_lamda2))/hx,(   -quad_lamda1)/hy]*[cos(ParametersofFracture(k,4)+pi/2);sin(ParametersofFracture(k,4)+pi/2)];
        GradsigmabasisPonFracture( (CurrentIndex-1)*4 + 3 , :) = ...
            [(    quad_lamda2 ) /hx,(    quad_lamda1)/hy]*[cos(ParametersofFracture(k,4)+pi/2);sin(ParametersofFracture(k,4)+pi/2)];
        GradsigmabasisPonFracture( (CurrentIndex-1)*4 + 4 , :) = ...
            [(   -quad_lamda2 )/hx,( (1-quad_lamda1))/hy]*[cos(ParametersofFracture(k,4)+pi/2);sin(ParametersofFracture(k,4)+pi/2)];
%        basisPonFracture( (CurrentIndex-ParametersofFracture(k,11))*4+1 :...
%            (CurrentIndex-ParametersofFracture(k,11)+1)*4, end-1 ) = 0;
%        GradnubasisPonFracture( (CurrentIndex-ParametersofFracture(k,11))*4+1 :...
%            (CurrentIndex-ParametersofFracture(k,11)+1)*4, end-1 ) = 0;
%        GradsigmabasisPonFracture( (CurrentIndex-ParametersofFracture(k,11))*4+1 :...
%            (CurrentIndex-ParametersofFracture(k,11)+1)*4, end-1 ) = 0;
%        basisPonFracture((CurrentIndex-1)*4+1 : CurrentIndex*4, end) = 0;
%        GradnubasisPonFracture((CurrentIndex-1)*4+1 : CurrentIndex*4, end) = 0;
%        GradsigmabasisPonFracture((CurrentIndex-1)*4+1 : CurrentIndex*4, end) = 0;
end
end