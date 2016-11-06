function [a,phi,varargout] = xy2cos(x,y,f)
% [A,PHI] = XY2COS(X,Y,f) - returns the amplitude and phase of the least square fit to a cosine oscillation
%                       at frequencies f, in units of cycle per unit of x  
%                       y = a*cos(2*pi*f*x+phi) + e
%                       if f is a row vector of frequencies of size 1xn then A and PHI will be of size 1xn 
% [A,PHI,E] = XY2COS(X,Y,f) returns the residual E
% Shane Elipot 2016
  
  y2 = y(:);
  x2 = x(:);
  qg = find(~isnan(y));
  q = find(isnan(x2)|isnan(y2));
  x2(q) = [];
  y2(q) = [];
 
  G = [cos(2*pi*x2*f) -sin(2*pi*x2*f)];
  m = G\y2;
  m = reshape(m,2,size(f,2));
  
  phi = atan2(m(2,:),m(1,:));
  a = m(1,:)./cos(phi);
  
  if nargin == 3
    foo = NaN*ones(size(y));
    foo(qg) = y2 - G*m(:);
    varargout{1} = foo;
  end
  