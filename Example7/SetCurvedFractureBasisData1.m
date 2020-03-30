function [NumberofFractures,NumberofElementsPassed,FracturesPath,FracturesLength,...
    basisPonFracture,GradnubasisPonFracture,GradsigmabasisPonFracture,...
    PermeabilitynuonFracture,ThicknessonFracture,ds_dtonFracture] = ...
    SetCurvedFractureBasisData1(coordinates,elements4,quad_x)
% heart curve
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
% parametric equation 1
N_partition = 100003; t_start = 0; t_end = 2*pi; 
t_partition = linspace(t_start,t_end,N_partition)';
x_partition = -(13*cos(t_partition)-5*cos(2*t_partition)-2*cos(3*t_partition)-cos(4*t_partition))/50+1/2;
y_partition = (16*sin(t_partition).^3)/50+1/2; 
% parametric equation 2

% draw the parametric curve:
% hold on; plot(x_partition,y_partition,'r')
Ele_located = ( ceil((y_partition-ymin)/hy)-1 )*Cell_N + ceil((x_partition-xmin)/hx);
when_inter = find( diff([-1; Ele_located])~=0 );
Ele_passed = Ele_located(when_inter);
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
for ell = 1 : length(Ele_passed)
    if ell ==1
    time_in = t_start;
    time_out = t_partition(when_inter(ell+1)-1);
    elseif ell == length(Ele_passed)
    time_in = t_partition(when_inter(ell));
    time_out = t_end;
    else
    time_in = t_partition(when_inter(ell));
    time_out = t_partition(when_inter(ell+1)-1);
    end
t_sampling = time_in + (time_out - time_in)*[0;quad_x;1];
x_sampling = -(13*cos(t_sampling)-5*cos(2*t_sampling)-2*cos(3*t_sampling)-cos(4*t_sampling))/50+1/2; 
y_sampling = (16*sin(t_sampling).^3)/50+1/2;
quad_lamda1 = (x_sampling - coordinates(elements4(1,Ele_passed(ell)),1))/hx;
quad_lamda2 = (y_sampling - coordinates(elements4(1,Ele_passed(ell)),2))/hy;
PermeabilitynuonFracture(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + ell , 1 : length(quad_x) + 2) = ...
    1e7 * ones(size(  t_sampling  ));
ThicknessonFracture(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + ell , 1 : length(quad_x) + 2) = ...
    1e-2 * ones(size(  t_sampling  ));
dx_dt = -1/50*(-13*sin(t_sampling)+10*sin(2*t_sampling)+6*sin(3*t_sampling)+4*sin(4*t_sampling));
dy_dt = 16/50*(3*sin(t_sampling).^2 .* cos(t_sampling));
ds_dtonFracture(sum(NumberofElementsPassed(1:NumberofFractures-1,1)) + ell , 1 : length(quad_x) + 2) = ...
    sqrt(dx_dt.^2 + dy_dt.^2);
if sqrt(dx_dt.^2 + dy_dt.^2)<1e-10
FractureNu = [dx_dt*0 , dy_dt*0];
FractureSigma = - ([0,-1;1,0] * FractureNu')';    
else
FractureNu = [dx_dt , dy_dt]./sqrt(dx_dt.^2 + dy_dt.^2);
FractureSigma = - ([0,-1;1,0] * FractureNu')';    
end
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