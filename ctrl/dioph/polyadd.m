function R = polyadd(A,varargin)
% Adds two or more row vectors of possibly differing lengths.
% R(x) = A(x) + B(x) + C(x) + · · · .
np = length(varargin); nmax = length(A);
for i=1:np % find maximum order in the input argument list
nmax = max(nmax,length(varargin{i}));
end % if
R = [zeros(1,nmax-length(A)),A];
for i=1:np
varargin{i} = [zeros(1,nmax - length(varargin{i}),1), varargin{i}];
R = R+varargin{i};
end % if
return