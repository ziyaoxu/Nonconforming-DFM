function [ Tarea , Elenth , Jacobimat ] = ComputeGeometryData( elements3 ,coordinates )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Tarea=zeros(1,size(elements3,2)); 
Elenth=zeros(size(elements3));
Jacobimat = zeros(2,2,size(elements3,2)); %  [Dr/Dx,Dr/Dy ; Ds/Dx,Ds/Dy] 
for jj = 1 : size(elements3,2)
    for ii = 1 : 3
        vertex1 = mod(ii,3)+1 ; 
        vertex2 = mod(ii+1,3)+1 ;
        edge_vec =  coordinates(elements3(vertex2,jj),:)' - ...
                coordinates(elements3(vertex1,jj),:)' ;
        Elenth(ii,jj)=norm(edge_vec);
    end
    Tarea(jj) = det([1,1,1;coordinates(elements3(:,jj),:)'])/2 ;
    Jacobimat(:,:,jj) = [0,1,0;0,0,1]*...
        ([1,1,1;coordinates(elements3(:,jj),:)']\[0,0;1,0;0,1]);
end
end

