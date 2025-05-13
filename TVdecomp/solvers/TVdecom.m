function [u,txtr] = TVdecom(x0,opts,Tol)
% % % This function computing the image decomposition problem for 
% % % clean image, i.e., the model 
% % % min || |u| ||_1+ tau ||(u+div g)-f||^2_2 + mu|| |g| ||_p
% % % The algorithm is the ADMM 3-blocks (convergence is ambiguous)   
% % % Input: 
% % % x0  -- target image 
% % % p   -- p-norm of term || |g| ||_p, p can be {1,2,infty}
% % % opts - trade off parameters (\tau,\mu) and penalty parameters (\beta_i)
% % % Output: u   -- cartoon part
% % %         txtr - texture part 

p     = 1;
tau   = opts.tau;  
mu    = opts.mu; 
beta  = opts.beta;  
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
ZERO = zeros(n1,n2,n3); 
u    = ZERO;
x1   = ZERO;  x2  = ZERO;
y    = ZERO;
g1   = ZERO;  g2  = ZERO;
txtr = ZERO; 
lbd11= ZERO; lbd12 = ZERO;
lbd2 = ZERO; 
lbd31= ZERO; lbd32 = ZERO;  
MDu  = beta*(abs(d1h).^2+abs(d2h).^2) + beta;
MDy  = beta + 1;
Mg11 = beta*abs(d1h).^2 + beta;
Mg22 = beta*abs(d2h).^2 + beta;
Mg12 = beta*d1h.*conj(d2h);
Mg21 = beta*d2h.*conj(d1h);
Mdet = beta*beta*(abs(d1h).^2 + abs(d2h).^2) + beta*beta; 
HTx0 = x0;

for k = 1:MaxIt    
    %%%% \tilde u
    tep1 = beta*x1 + lbd11 + beta*g1;
    tep2 = beta*x2 + lbd12 + beta*g2;
    Temp = PTx(tep1) + PTy(tep2) + lbd2 + beta*y;
    u    = real(ifft2(fft2(Temp)./MDu));
  
    %%%% \tilde z
    tep1 = g1 + lbd31/beta;
    tep2 = g2 + lbd32/beta;
    [poz1,poz2] = Project(tep1,tep2,mu/beta,p);
    z1  = tep1 - poz1;
    z2  = tep2 - poz2;
 
    %%%% \tilde x
    sk1 = Px(u)- lbd11/beta;
    sk2 = Py(u)- lbd12/beta;
    nsk = sqrt(sk1.^2 + sk2.^2); nsk(nsk==0)=1;
    nsk = max(1 - tau./(beta*nsk),0);
    x1  = nsk.*sk1;
    x2  = nsk.*sk2;
    
    %%%% \tilde g
    Temp = beta*(u-y) - lbd2;
    b1   = Px(Temp) + beta*z1 - lbd31;
    b2   = Py(Temp) + beta*z2 - lbd32;
    nume1= Mg22.*fft2(b1) - Mg12.*fft2(b2);
    nume2= Mg11.*fft2(b2) - Mg21.*fft2(b1);    
    g1   = real( ifft2(nume1./Mdet) );
    g2   = real( ifft2(nume2./Mdet) );
    
    %%%% \tilde y
    txtr = -PTx(g1)-PTy(g2);
    Temp = HTx0 + beta*(u+txtr) - lbd2;
    y    = Temp/MDy;
    
    %%%% update lagrange multipliers
    lbd11 = lbd11 - beta*(Px(u)-x1);
    lbd12 = lbd12 - beta*(Py(u)-x2);
    lbd2  = lbd2  - beta*(u+txtr-y);
    lbd31 = lbd31 - beta*(z1-g1);
    lbd32 = lbd32 - beta*(z2-g2);
    
    er=norm(beta*(Px(u)-x1))+norm(beta*(Py(u)-x2))
    
    if k>2 && er< Tol
        disp('successful');
        disp(k);
        break;    
    end
       
    subplot(221); imshow(x0,[]);   title('target image')
    subplot(223); imshow(u,[]);  title(['cartoon, Iter:',num2str(k)])
    subplot(224); imshow(txtr-min(txtr(:))); 
        
end




