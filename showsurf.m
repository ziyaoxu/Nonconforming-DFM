function h = showsurf( coordinates , elements , u )
%SHOWSURF show the surf of S  
% S is the function value in the triangular mesh 
u = full(u);
if size(u,2)~=1
    error('data must be a vector ,not a matrix ! ')
end
if size(elements,2)>size(elements,1)
   elements = elements' ;  
end
% interpolate the data (then draw it) :
h = trisurf(elements,coordinates(:,1),coordinates(:,2),u,'facecolor','interp');
%set(h, 'LineStyle','none');
view(10,40);
xlabel('x'),ylabel('y')  % axis name 
h = gcf ; 
% title('Ñ¹Á¦·Ö²¼')
end