function w = wij(m,n)
% WIJ    Produces swap matrix
%
%   A = WIJ(N) produces an N^2-by-N^2 swap matrix.
%   A = WIJ(M,N) produces an MN-by-MN swap matrix.

% $Id: wij.m 7 2010-04-09 12:11:31Z hsqi $

if nargin == 1
    n=m;
end

if m <= 1 || n <= 1
    w = 1;
    return;
end

d = m*n;
w = zeros(d);
for k = 1:d
    j = mod(k,n);
    if j == 0
        j = n;
    end
    i = (k-j)/n+1;
    w((j-1)*m+i,k) = 1;
end;
