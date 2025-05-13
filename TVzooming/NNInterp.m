function x = NNInterp(y,n)
[ny,nx] = size(y);
x = zeros(ny*n,nx*n);
for i=1:ny
    for j=1:nx
        %Average of adjacent pixels
        for k=1:n
            for l=1:n
                ii = k+(i-1)*n;
                jj = l+(j-1)*n;
                x(ii,jj) = y(i,j);
            end
        end
    end
end

