function A = Bilinear2(m,n,x,y,h1,h2)
% This functions uses Bilinear interpolation 
% to produce a sparse matrix of weights
% ***********************************************************************
%   The difference between Bilinear2 and Bilinear is that Bilinear2 takes
%   into account the boundary and does NOT assume that the image is
%   embedded in a zero boundary.
% ***********************************************************************
%Input:  m,n   - size of the HR image
%        x,y   - displaced pixels
%        h1,h2 - pixel size%
%Output: A     - sparse matrix of coefficients

% Convert x and y to the coordinate system 1:m, 1:n
x = x/h1 + 1/2; y = y/h2 + 1/2;
% Vectorized version
%     Th = zeros(m*n,1);
    j=floor(x(:)); xi  = x(:)-j;  
    k=floor(y(:)); eta = y(:)-k;    
    ind1 = find(1<=j & j<m & 1<=k & k<n);
    jk = j(ind1) + (k(ind1)-1)*m;
    j1k = j(ind1) + 1+ (k(ind1)-1)*m;
    jk1 = j(ind1) + k(ind1)*m;
    j1k1 = j(ind1) + 1 + k(ind1)*m;    
    ii = [ind1;ind1;ind1;ind1];
    jj = [jk;j1k;jk1;j1k1];
    ss = [(1-xi(ind1)).*(1-eta(ind1)); (xi(ind1)).*(1-eta(ind1)); ...
         (1-xi(ind1)).*(eta(ind1));   (xi(ind1)).*(eta(ind1))];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind2 = find(j==m & 1<k & k<n);
    jk = j(ind2) + (k(ind2)-1)*m;
    jk1 = j(ind2) + k(ind2)*m;       
    ii = [ii; ind2;ind2];
    jj = [jj; jk;jk1];
    ss = [ss; (1-eta(ind2)); (eta(ind2))];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind3 = find(1<j & j<m & k==n);
    jk = j(ind3) + (k(ind3)-1)*m;
    j1k = j(ind3) + 1+ (k(ind3)-1)*m;   
    ii = [ii; ind3;ind3];
    jj = [jj; jk;j1k];
    ss = [ss; (1-xi(ind3)); (xi(ind3))];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind4 = find(j==m & k==n);
    jk = j(ind4) + (k(ind4)-1)*m;
    ii = [ii; ind4];
    jj = [jj; jk];
    ss = [ss; ind4./ind4];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ind5 = find(j==m & k==1);
    jk = j(ind5) + (k(ind5)-1)*m;
    ii = [ii; ind5];
    jj = [jj; jk];
    ss = [ss; ind5./ind5];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind6 = find(j==1 & k==n);
    jk = j(ind6) + (k(ind6)-1)*m;
    ii = [ii; ind6];
    jj = [jj; jk];
    ss = [ss; ind6./ind6];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    A = sparse(ii,jj,ss,n*m,n*m);
    return;
    
