function Ax = LSPZoom(H,HT,Dss,x,n,m,tau,beta2)
xx = reshape(x,n,m);
Ax = 2*tau*HT(Dss.*H(xx))+beta2*xx;
Ax = Ax(:);
