function [coordinates,elements4,EtoEmap] = RectangularMesh( xmin,xmax,ymin,ymax,N,M )
%RECTANGULARMESH Summary of this function goes here
%   Detailed explanation goes here
XX=linspace(xmin,xmax,N+1); YY=linspace(ymin,ymax,M+1); 
[X,Y]=meshgrid(XX,YY);      % 纵线横线交错生成方格网节点的x、y坐标，分别存入X、Y矩阵
coordinates=zeros((N+1)*(M+1),2); elements4=zeros(4,N*M); EtoEmap=zeros(4,N*M);
X=X'; Y=Y';    
coordinates(:,1)=X(:);coordinates(:,2)=Y(:); 
for I=1:M 
    for J=1:N
        num = (I-1)*N+J; 
        start=(I-1)*(N+1)+J;    
        elements4(:,num)=[start;start+1;start+N+2;start+N+1];
        EtoEmap(:,num)=[num-N;num+1;num+N;num-1];
    end 
end
EtoEmap(1,1:N)=-1;
EtoEmap(3,(M-1)*N+1:M*N)=-1;
EtoEmap(4,1:N:end)=-1;
EtoEmap(2,N:N:end)=-1;
% periodic boundary condition :
EtoEmap(1,1:N)=(M-1)*N+1:M*N;
EtoEmap(3,(M-1)*N+1:M*N)=1:N;
EtoEmap(4,1:N:end)=N:N:M*N;
EtoEmap(2,N:N:end)=1:N:M*N;
end

