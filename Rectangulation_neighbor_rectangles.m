function EtoEmap = Rectangulation_neighbor_rectangles( elements4 )
%RECTANGULATION_NEIGHBOR_RECTANGLES Summary of this function goes here
%   Detailed explanation goes here
% x direction is divided into N parts ,y direction is divided into M parts:
NM = size(elements4,2);
N = elements4(4,1)-elements4(1,1)-1;
M = NM/N;
EtoEmap=zeros(4,N*M);
for I=1:M 
    for J=1:N
        num = (I-1)*N+J; 
        EtoEmap(:,num)=[num-N;num+1;num+N;num-1];
    end 
end
EtoEmap(1,1:N)=-1;
EtoEmap(3,(M-1)*N+1:M*N)=-1;
EtoEmap(4,1:N:end)=-1;
EtoEmap(2,N:N:end)=-1;
end

