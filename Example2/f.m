function f = f(xy)
x = xy(:,1);
y = xy(:,2);
r = sqrt(x.^2+y.^2);
%f = 2*cos(x).*cos(y);
%f = zeros(size(x));
f = 0*2*exp(-r);
end
