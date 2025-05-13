function [u,txtr] = DnoisetauGB(x0,p,opts)
%%% image decomposition on clean image
%%% min || |u| ||_1+ tau ||(u+div g)-f||^2_2 + mu|| |g| ||_p
%%% Algorithm: ADMM-GB(B.S.He,M.Tao,X.M.Yuan,SIAM J Optim,22,313-340(2012)   
%%% Input: 
%%% x0  -- target image 
%%% p   -- p-norm of term || |g| ||_p, p can be {1,2,infty}
%%% opts - trade off parameters (\tau,\mu) and penalty parameters (\beta_i)
%%% Output: u   -- cartoon part
%%%         txtr - texture part 

tau   = opts.tau;    mu    = opts.mu; 
beta1 = opts.beta1;  beta2 = opts.beta2; beta3 = opts.beta3;
MaxIt = opts.MaxIt;
[n1,n2,n3] = size(x0);
%%%%%%%%%%%%%%%%% Periodic  boundary condtion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1h= zeros(n1,n2,n3); d1h(1,1,:) = -1; d1h(n1,1,:) = 1; d1h = fft2(d1h);
d2h= zeros(n1,n2,n3); d2h(1,1,:) = -1; d2h(1,n2,:) = 1; d2h = fft2(d2h);
Px = @(x) [x(2:n1,:,:)-x(1:n1-1,:,:); x(1,:,:)-x(n1,:,:)]; %%\nabla_1 x 
Py = @(x) [x(:,2:n2,:)-x(:,1:n2-1,:), x(:,1,:)-x(:,n2,:)]; %%\nabla_2 y
PTx= @(x) [x(n1,:,:)-x(1,:,:); x(1:n1-1,:,:)-x(2:n1,:,:)]; %%\nabla_1^T x 
PTy= @(x) [x(:,n2,:)-x(:,1,:), x(:,1:n2-1,:)-x(:,2:n2,:)]; %%\nabla_2^T y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% initinalization
u   = zeros(n1,n2,n3);   
x1  = zeros(n1,n2,n3);    x2 = zeros(n1,n2,n3); 
y   = zeros(n1,n2,n3); 
z1  = zeros(n1,n2,n3);    z2 = zeros(n1,n2,n3);  
txtr= zeros(n1,n2,n3); 
lbd11 = zeros(n1,n2,n3);  lbd12 = zeros(n1,n2,n3); 
lbd2  = zeros(n1,n2,n3); 
lbd31 = zeros(n1,n2,n3);  lbd32 = zeros(n1,n2,n3); 
MDu  = beta1*(abs(d1h).^2+abs(d2h).^2) + beta2;
MDy  = beta2 + 1;
Mg11 = beta2*abs(d1h).^2 + beta3;
Mg22 = beta2*abs(d2h).^2 + beta3;
Mg12 = beta2*d1h.*conj(d2h);
Mg21 = beta2*d2h.*conj(d1h);
Mdet = beta2*beta3*(abs(d1h).^2 + abs(d2h).^2) + beta3*beta3; 
HTx0 = x0;

for k = 1:MaxIt
  
    %%%% \tilde g
    Temp = beta2*(u-y) - lbd2;
    b1   = Px(Temp) + beta3*z1 + lbd31;
    b2   = Py(Temp) + beta3*z2 + lbd32;
    nume1= Mg22.*fft2(b1) - Mg12.*fft2(b2);
    nume2= Mg11.*fft2(b2) - Mg21.*fft2(b1);    
    gn1  = real( ifft2(nume1./Mdet) );
    gn2  = real( ifft2(nume2./Mdet) );
    
    %%%% \tilde u
    tep1 = beta1*x1 + lbd11 + beta2*gn1;
    tep2 = beta1*x2 + lbd12 + beta2*gn2;
    Temp = PTx(tep1) + PTy(tep2) + lbd2 + beta2*y;
    un   = real(ifft2(fft2(Temp)./MDu));
    
    %%%% \tilde x
    sk1  = Px(un) - lbd11/beta1;
    sk2  = Py(un) - lbd12/beta1;
    nsk  = sqrt(sk1.^2 + sk2.^2); nsk(nsk==0)=1;
    nsk  = max(1 - tau./(beta1*nsk),0);
    xn1  = nsk.*sk1;
    xn2  = nsk.*sk2;
    
    %%%% \tilde y
    txtr= -PTx(gn1)-PTy(gn2);
    Temp= HTx0 + beta2*(un+txtr) - lbd2;
    yn  = Temp/MDy;
    
    %%%% \tilde z
    tep1 = gn1 - lbd31/beta3;
    tep2 = gn2 - lbd32/beta3;
    [poz1,poz2] = Project(tep1,tep2,mu/beta3,p);
    zn1  = tep1 - poz1;
    zn2  = tep2 - poz2;
 
    %%%% update lagrange multipliers
    lbdn11 = lbd11 - beta1*(Px(un)-xn1);
    lbdn12 = lbd12 - beta1*(Px(un)-xn2);
    lbdn2  = lbd2  - beta2*(un+txtr-yn);
    lbdn31 = lbd31 - beta3*(gn1-zn1);
    lbdn32 = lbd32 - beta3*(gn2-zn2);
        
    %%%% correction step
    Temp= beta1*(PTx(xn1-x1) + PTy(xn2-x2)) + beta2*(yn-y);
    u   = un + real( ifft2(fft2(Temp)./MDu));
    x1  = xn1;      x2  = xn2;
    y   = yn;
    z1  = zn1;      z2  = zn2;
    lbd11 = lbdn11; lbd12 = lbdn12;
    lbd2  = lbdn2; 
    lbd31 = lbdn31; lbd32 = lbdn32;
    
    if mod(k,5)==0
        figure(10);   %%%% ÏÔÊ¾½á¹û
        subplot(221); imshow(x0,[]);    title('target image')
        subplot(222); imshow(u,[]);  title('cartoon')
        subplot(223); imshow(txtr-min(txtr(:)));  title('texture'); pause(0.1)
    end
end




