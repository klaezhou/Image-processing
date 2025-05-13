function [u,txtr,SNR]=Inpainttau(x0,S,p,opts,I)
% % % This function computing the image decomposition problem for 
% % % image with missing pixels, i.e., the model 
% % % min || |u| ||_1+ tau ||S(u+div g)-f||^2_2 + mu|| |g| ||_p
% % % where S represents the mask
% % % The algorithm is the ADMM 3-blocks (convergence is ambiguous)   
% % % Input: 
% % % x0  -- target image 
% % % p   -- p-norm of term || |g| ||_p, p can be {1,2,infty}
% % % opts - trade off parameters (\tau,\mu) and penalty parameters (\beta_i)
% % % Output: u   -- cartoon part
% % %         txtr - texture part 
% % %         SNR  - SNR of inpainted image
% % % Reference: M.K. Ng, X.M. Yuan, and W.X. Zhang, "Coupled Variational
% % % Image Decomposition and Restoration Model for Blurred Cartoon-Plus-Texture
% % % Images With Missing Pixels", IEEE TIP, 22(6),2013, pp. 2233-2246 
% % % Written by W.X. Zhang, wenxing84@gmail.com
tau  = opts.tau; 
mu   = opts.mu; 
beta1= opts.beta1;
beta2= opts.beta2; 
beta3= opts.beta3;
MaxIt= opts.MaxIt;


[n1,n2,n3] = size(x0);
%%%%%%%%%%%%%%%%% Periodic  boundary condtion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1h = zeros(n1,n2,n3); d1h(1,1,:) = -1; d1h(n1,1,:) = 1; d1h = fft2(d1h);
d2h = zeros(n1,n2,n3); d2h(1,1,:) = -1; d2h(1,n2,:) = 1; d2h = fft2(d2h);
Px  = @(x) [x(2:n1,:,:)-x(1:n1-1,:,:); x(1,:,:)-x(n1,:,:)]; %%\nabla_1 x 
Py  = @(x) [x(:,2:n2,:)-x(:,1:n2-1,:), x(:,1,:)-x(:,n2,:)]; %%\nabla_2 y
PTx = @(x) [x(n1,:,:)-x(1,:,:); x(1:n1-1,:,:)-x(2:n1,:,:)]; %%\nabla_1^T x 
PTy = @(x) [x(:,n2,:)-x(:,1,:), x(:,1:n2-1,:)-x(:,2:n2,:)]; %%\nabla_2^T y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% initinalization
ZERO = zeros(n1,n2,n3); 
u    = ZERO;
x1   = ZERO;   x2    = ZERO; 
y    = ZERO; 
g1   = ZERO;   g2    = ZERO; 
lbd11= ZERO;   lbd12 = ZERO;
lbd2 = ZERO; 
lbd31= ZERO;   lbd32 = ZERO;  

MDu  = beta1*(abs(d1h).^2+abs(d2h).^2) + beta2;
MDy  = beta2 + S;
Mg11 = beta2*abs(d1h).^2 + beta3;
Mg22 = beta2*abs(d2h).^2 + beta3;
Mg12 = beta2*d1h.*conj(d2h);
Mg21 = beta2*d2h.*conj(d1h);
Mdet = beta2*beta3*(abs(d1h).^2 + abs(d2h).^2) + beta3*beta3; 
STx0 = S.*x0;

Im   = norm(I(:));
SNR  = zeros(1,MaxIt);
txtr = zeros(size(x0));

for k = 1:MaxIt
        
    SNR(k) = 20*log10(Im/norm(u(:)+txtr(:)-I(:)));
    
    %%%% \tilde u
    tep1 = beta1*x1 + lbd11 + beta2*g1;
    tep2 = beta1*x2 + lbd12 + beta2*g2;
    Temp = PTx(tep1) + PTy(tep2) + lbd2 + beta2*y;
    un   = real(ifft2(fft2(Temp)./MDu));
  
    %%%% \tilde z
    tep1 = g1 + lbd31/beta3;
    tep2 = g2 + lbd32/beta3;
    [poz1,poz2] = Project(tep1,tep2,mu/beta3,p);
    zn1  = tep1 - poz1;
    zn2  = tep2 - poz2;
 
    %%%% \tilde x
    dxun = Px(un);
    dyun = Py(un);
    sk1  = dxun - lbd11/beta1;
    sk2  = dyun - lbd12/beta1;
    nsk  = sqrt(sk1.^2 + sk2.^2); nsk(nsk==0)=1;
    nsk  = max(1 - (tau/beta1)./nsk,0);
    xn1  = nsk.*sk1;
    xn2  = nsk.*sk2;
    
    %%%% \tilde g
    Temp = beta2*(un-y) - lbd2;
    b1   = Px(Temp) + beta3*zn1 - lbd31;
    b2   = Py(Temp) + beta3*zn2 - lbd32;
    nume1= Mg22.*fft2(b1) - Mg12.*fft2(b2);
    nume2= Mg11.*fft2(b2) - Mg21.*fft2(b1);    
    gn1  = real( ifft2(nume1./Mdet) );
    gn2  = real( ifft2(nume2./Mdet) );
    
    %%%% \tilde y
    txtr = -PTx(gn1)-PTy(gn2);
    Temp = STx0 + beta2*(un+txtr) - lbd2;
    yn   = Temp./MDy;
    
    %%%% update lagrange multipliers
    lbdn11 = lbd11 - beta1*(dxun-xn1);
    lbdn12 = lbd12 - beta1*(dyun-xn2);
    lbdn2  = lbd2  - beta2*(un+txtr-yn);
    lbdn31 = lbd31 - beta3*(zn1-gn1);
    lbdn32 = lbd32 - beta3*(zn2-gn2);
    
    
    %%%% update variable
    u    = un; 
    g1   = gn1;    g2 = gn2; 
    x1   = xn1;    x2 = xn2; 
    y    = yn; 
    z1   = zn1;    z2 = zn2;
    lbd11= lbdn11; lbd12 = lbdn12;
    lbd2 = lbdn2; 
    lbd31= lbdn31; lbd32 = lbdn32;
     
end



% imshow(u,'border','tigh'),rectangle('Position', [180 10 70 70],'EdgeColor','w','lineWidth',3)
% imshow(u(10:10+70,180:180+70,:),'border','tigh')
% imshow(x0,'border','tigh'),
% rectangle('Position', [390 20 60 80],'EdgeColor','r','lineWidth',3)
% rectangle('Position', [340 420 60 80],'EdgeColor','r','lineWidth',3)
% a=[u(20:20+80,390:390+60,:),ones(81,5,3),0.5+txtr(20:20+80,390:390+60,:)];
% b=[x0(420:420+80,340:340+60,:),u(420:420+80,340:340+60,:),0.5+txtr(420:420+80,340:340+60,:)];
% imshow(a,'border','tigh')

% imshow(x0,'border','tigh'),
% rectangle('Position', [20 60 60 80],'EdgeColor','r','lineWidth',2)
