function [NumberofFractures,NumberofElementsPassed,FracturesPath,FracturesLength,...
    basisPonFracture,GradnubasisPonFracture,GradsigmabasisPonFracture,...
    PermeabilitynuonFracture,ThicknessonFracture,ds_dtonFracture] = ...
    SetCurvedFractureBasisData(coordinates,elements4,quad_x)
% time : 2020.1.26 - 2020.1.28
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
% Global geometric information:
EtoEmap = Rectangulation_neighbor_rectangles( elements4 );
xmax = max(coordinates(:,1)); xmin = min(coordinates(:,1));
ymax = max(coordinates(:,2)); ymin = min(coordinates(:,2));
% x direction is divided into N parts ,y direction is divided into M parts:
Cell_N = (elements4(4,1)-elements4(1,1)-1) ; Cell_M = size(elements4,2)/Cell_N ;
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ;
% set fractures data:
NumberofFractures = 0;
NumberofElementsPassed = zeros(NumberofFractures,1);
% list of large prime numbers used to choose N :
% 1009,2003,3001,5003,10007,20011,30011,40009,50021,60013,70001,80021,90001,100003
% ////////////////////////parametric curve 1://////////////////////////////
NumberofFractures = NumberofFractures + 1;
N_partition = 100003; t_start = 0; t_end = 2*pi; 
t_partition = linspace(t_start,t_end,N_partition)';
Radius = 0.25; x_center = 0.5+(rand-0.5)*1e-4; y_center = 0.5+(rand-0.5)*1e-4;
x_partition = Radius*cos(t_partition) + x_center; y_partition = Radius*sin(t_partition) + y_center;
%//////////////////////////////////////////////////////////////////////////
% draw the parametric curve:
hold on; plot(x_partition,y_partition,'r')
Ele_located = ( ceil((y_partition-ymin)/hy)-1 )*Cell_N + ceil((x_partition-xmin)/hx);
where_inter = find( diff([-1; Ele_located])~=0 );
Ele_passed = Ele_located(where_inter);
NumberofElementsPassed(NumberofFractures,1) = length(Ele_passed);
FracturesPath( sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + [1:length(Ele_passed)] , 1 ) = ...
    Ele_passed;
% calculate the following:
FracturesLength(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + [1:length(Ele_passed)],1) = 0;
basisPonFracture(4*sum(NumberofElementsPassed(1:NumberofFractures-1,1))+[1:4*length(Ele_passed)] , length(quad_x) + 2) = 0;
GradnubasisPonFracture(4*sum(NumberofElementsPassed(1:NumberofFractures-1,1))+[1:4*length(Ele_passed)] , length(quad_x) + 2) = 0;
GradsigmabasisPonFracture(4*sum(NumberofElementsPassed(1:NumberofFractures-1,1))+[1:4*length(Ele_passed)] , length(quad_x) + 2) = 0;
PermeabilitynuonFracture(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + [1:length(Ele_passed)] , length(quad_x) + 2) = 0;
ThicknessonFracture(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + [1:length(Ele_passed)] , length(quad_x) + 2) = 0;
ds_dtonFracture(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + [1:length(Ele_passed)] , length(quad_x) + 2) = 0;
for ell = 1:length(Ele_passed) 
    if ell == 1
        edge_in = -1;
        edge_out = find( EtoEmap(:,Ele_passed(ell)) == Ele_passed(ell+1) );
        time_in = t_start;
    elseif ell == length(Ele_passed)
        edge_in  = find( EtoEmap(:,Ele_passed(ell)) == Ele_passed(ell-1) );
        edge_out = -1;
        time_out = t_end;
    else
        edge_in  = find( EtoEmap(:,Ele_passed(ell)) == Ele_passed(ell-1) );
        edge_out = find( EtoEmap(:,Ele_passed(ell)) == Ele_passed(ell+1) );
    end
% ////////////////////////parametric curve 1://////////////////////////////
    if (edge_in == 1) || (edge_in == 3) 
        % t=y_inverse(y0)
        theta = asin( ((coordinates(elements4(1,Ele_passed(ell)),2)+(edge_in-1)/2*hy) - y_center) / Radius );
        t_candidates = [theta;pi-theta;2*pi+theta];
        time_in = t_candidates(t_candidates>=t_partition(where_inter(ell)-1) & t_candidates<=t_partition(where_inter(ell)));
    end
    if (edge_in == 2) || (edge_in == 4)
        % t=x_inverse(x0)
        theta = acos( ((coordinates(elements4(1,Ele_passed(ell)),1)+(4-edge_in)/2*hx) - x_center) / Radius );
        t_candidates = [theta;2*pi-theta];
        time_in = t_candidates(t_candidates>=t_partition(where_inter(ell)-1) & t_candidates<=t_partition(where_inter(ell)));
    end
    if (edge_out == 1) || (edge_out == 3) 
        % t=y_inverse(y0)
        theta = asin( ((coordinates(elements4(1,Ele_passed(ell)),2)+(edge_out-1)/2*hy) - y_center) / Radius );
        t_candidates = [theta;pi-theta;2*pi+theta];
        time_out = t_candidates(t_candidates>=t_partition(where_inter(ell+1)-1) & t_candidates<=t_partition(where_inter(ell+1)));
    end
    if (edge_out == 2) || (edge_out == 4)
        % t=x_inverse(x0)
        theta = acos( ((coordinates(elements4(1,Ele_passed(ell)),1)+(4-edge_out)/2*hx) - x_center) / Radius );
        t_candidates = [theta;2*pi-theta];
        time_out = t_candidates(t_candidates>=t_partition(where_inter(ell+1)-1) & t_candidates<=t_partition(where_inter(ell+1)));
    end
t_sampling = time_in + (time_out - time_in)*[0;quad_x;1];
x_sampling = Radius*cos(t_sampling) + x_center; 
y_sampling = Radius*sin(t_sampling) + y_center;
quad_lamda1 = (x_sampling - coordinates(elements4(1,Ele_passed(ell)),1))/hx;
quad_lamda2 = (y_sampling - coordinates(elements4(1,Ele_passed(ell)),2))/hy;
PermeabilitynuonFracture(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + ell , 1 : length(quad_x) + 2) = ...
    1e7 * ones(size(  t_sampling  ));
ThicknessonFracture(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + ell , 1 : length(quad_x) + 2) = ...
    0.005 * ones(size(  t_sampling  ));
ds_dtonFracture(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + ell , 1 : length(quad_x) + 2) = ...
    Radius * ones(size(  t_sampling  ));
FractureNu = [-sin(t_sampling) , cos(t_sampling)];
FractureSigma = - ([0,-1;1,0] * FractureNu')';
%//////////////////////////////////////////////////////////////////////////
FracturesLength(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + ell,1) = time_out - time_in;   
basisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 1 , 1 : length(quad_x) + 2) = ...
    (1-quad_lamda1).*(1-quad_lamda2);
basisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 2 , 1 : length(quad_x) + 2) = ...
    quad_lamda1.*(1-quad_lamda2);
basisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 3 , 1 : length(quad_x) + 2) = ...
    quad_lamda1.*quad_lamda2;
basisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 4 , 1 : length(quad_x) + 2) = ...
    (1-quad_lamda1).*quad_lamda2;

GradnubasisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 1 , 1 : length(quad_x) + 2) = ...
    sum([(-(1-quad_lamda2))/hx,(-(1-quad_lamda1))/hy].*FractureNu,2);
GradnubasisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 2 , 1 : length(quad_x) + 2) = ...
    sum([( (1-quad_lamda2))/hx,(   -quad_lamda1)/hy].*FractureNu,2);
GradnubasisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 3 , 1 : length(quad_x) + 2) = ...
    sum([(    quad_lamda2 ) /hx,(    quad_lamda1)/hy].*FractureNu,2);
GradnubasisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 4 , 1 : length(quad_x) + 2) = ...
    sum([(   -quad_lamda2 )/hx,( (1-quad_lamda1))/hy].*FractureNu,2);

GradsigmabasisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 1 , 1 : length(quad_x) + 2) = ...
    sum([(-(1-quad_lamda2))/hx,(-(1-quad_lamda1))/hy].*FractureSigma,2);
GradsigmabasisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 2 , 1 : length(quad_x) + 2) = ...
    sum([( (1-quad_lamda2))/hx,(   -quad_lamda1)/hy].*FractureSigma,2);
GradsigmabasisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 3 , 1 : length(quad_x) + 2) = ...
    sum([(    quad_lamda2 ) /hx,(    quad_lamda1)/hy].*FractureSigma,2);
GradsigmabasisPonFracture(4*(sum(NumberofElementsPassed(1:NumberofFractures-1,1))+ell-1) + 4 , 1 : length(quad_x) + 2) = ...
    sum([(   -quad_lamda2 )/hx,( (1-quad_lamda1))/hy].*FractureSigma,2);
end


end