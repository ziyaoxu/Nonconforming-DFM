function [coordinates,elements,EtoEmap] = CubicMesh( xmin,xmax,ymin,ymax,zmin,zmax,N,M,L )
% time : 2020.2.11 - 2020.2.11
% author : xuziyao
%   本函数目前仅适用于立方体区域的均匀立方体剖分.
XX=linspace(xmin,xmax,N+1); YY=linspace(ymin,ymax,M+1); ZZ=linspace(zmin,zmax,L+1);
[X,Y,Z]=meshgrid(XX,YY,ZZ);      % 纵线横线交错生成方格网节点的x,y,z坐标，分别存入X,Y,Z矩阵
X = permute(X,[2,1,3]); Y = permute(Y,[2,1,3]); Z = permute(Z,[2,1,3]); 
coordinates=zeros((N+1)*(M+1)*(L+1),3); elements=zeros(8,N*M*L); EtoEmap=zeros(6,N*M*L);
coordinates(:,1)=X(:);coordinates(:,2)=Y(:);coordinates(:,3)=Z(:); 
for K=1:L
    for I=1:M 
        for J=1:N
            num = (K-1)*N*M + (I-1)*N + J; 
            start=(K-1)*(N+1)*(M+1) + (I-1)*(N+1) + J;    
            elements(:,num)=[start;start+1;start+N+2;start+N+1;...
                start+(N+1)*(M+1);start+1+(N+1)*(M+1);start+N+2+(N+1)*(M+1);start+N+1+(N+1)*(M+1)];
            EtoEmap(:,num)=[num-N*M;num+N*M;num-1;num+1;num-N;num+N];
        end 
    end
end
EtoEmap(1,1:N*M)=-1;
EtoEmap(2,end-N*M+1:end)=-1;
EtoEmap(3,1:N:end)=-1;
EtoEmap(4,N:N:end)=-1;
EtoEmap(5,reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N],L,1)',[],1))=-1 ;
EtoEmap(6,reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N]+(M-1)*N,L,1)',[],1))=-1;
% periodic boundary condition :
isperiodic = 0;
if isperiodic == 1
EtoEmap(1,1:N*M)=(1:N*M)+(L-1)*N*M;
EtoEmap(2,end-N*M+1:end)=1:N*M;
EtoEmap(3,1:N:end)=(1:N:N*M*L)+N-1;
EtoEmap(4,N:N:end)=(1:N:N*M*L);
EtoEmap(5,reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N],L,1)',[],1))=...
    reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N]+(M-1)*N,L,1)',[],1);
EtoEmap(6,reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N]+(M-1)*N,L,1)',[],1))=...
    reshape(repmat([0:L-1]'*N*M,1,N)' + repmat([1:N],L,1)',[],1);
end

end

