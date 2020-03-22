function boundary_edge = boundedges_Rec(elements4)
% time : 2019.3.21 - 2019.3.21
% author : xuziyao
% Find boundary edges of rectangular meshes
NM = size(elements4,2);
N = elements4(4,1)-elements4(1,1)-1;
M = NM/N;
boundary_edge = zeros(2*N+2*M,2);
boundary_edge(1:N,:)=[1:N;2:N+1]';
boundary_edge(N+1:N+M,:)=(N+1)*[1:M;2:M+1]';
boundary_edge(N+M+1:2*N+M,:)=(N+1)*(M+1)-[0:N-1;1:N]';
boundary_edge(2*N+M+1:2*N+2*M,:)=1+(N+1)*[M:-1:1;M-1:-1:0]';
end
