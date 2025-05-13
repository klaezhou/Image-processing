function [x,PSNR] = TVinpaint(x0,h,S,opts,I)

beta = opts.beta;
tau  = opts.tau;
MaxIt= opts.MaxIt;  

[n1,n2,n3] = size(x0);

%%%%Periodic  boundary condtion
siz = size(h);
center = [fix(siz(1)/2+1),fix(siz(2)/2+1)];
P  = zeros(n1,n2,n3); 
for i =1:n3; P(1:siz(1),1:siz(2),i) = h; end
D  = fft2(circshift(P,1-center));
H  = @(x) real(ifft2(D.*fft2(x)));       %%%% Blur operator.  B x 
HT = @(x) real(ifft2(conj(D).*fft2(x))); %%%% Transpose of blur operator.

%%%%%%%%%%%%%%%%% Í¼ÏñµÄÌÝ¶È %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1h = zeros(n1,n2,n3); d1h(1,1,:) = -1; d1h(n1,1,:) = 1; d1h = fft2(d1h);
d2h = zeros(n1,n2,n3); d2h(1,1,:) = -1; d2h(1,n2,:) = 1; d2h = fft2(d2h);
Px  = @(x) [x(2:n1,:,:)-x(1:n1-1,:,:); x(1,:,:)-x(n1,:,:)]; %%\nabla_1 x 
Py  = @(x) [x(:,2:n2,:)-x(:,1:n2-1,:), x(:,1,:)-x(:,n2,:)]; %%\nabla_2 y
PTx = @(x) [x(n1,:,:)-x(1,:,:); x(1:n1-1,:,:)-x(2:n1,:,:)]; %%\nabla_1^T x 
PTy = @(x) [x(:,n2,:)-x(:,1,:), x(:,1:n2-1,:)-x(:,2:n2,:)]; %%\nabla_2^T y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Initialization
x  = x0;
y1 = zeros(size(x0));    y2 = zeros(size(x0));
v  = zeros(size(x0)); 
z  = zeros(size(x0));
lbd1 = zeros(size(x0));    
lbd21= zeros(size(x0));  lbd22 = zeros(size(x0));
lbd3 = zeros(size(x0));

Dt = beta*(abs(D).^2+abs(d1h).^2+abs(d2h).^2);
Ds = beta*(1+S); 

PSNR = zeros(1,MaxIt); 

for k=1:MaxIt 
    
%%%%%%%%%%%% x
Temp = H(lbd1+beta*v)+PTx(lbd21+beta*y1)+PTy(lbd22+beta*y2);
x   = real(ifft2(fft2(Temp)./Dt));

%%%%%%%%%%%% v
Temp = S.*(lbd3+beta*z) - lbd1 + beta*H(x);
v    = Temp./Ds;

%%%%%%%%%%%  u =(y,z)
%%% y
sk1 = Px(x)-lbd21/beta;
sk2 = Py(x)-lbd22/beta;
nsk = sqrt(sk1.^2+sk2.^2); nsk(nsk==0) = 1;
nsk = max(nsk-tau/beta,0)./nsk;
y1  = sk1.*nsk;
y2  = sk2.*nsk;

%%% \tilde z
Sv = S.*v;
z  = (x0-lbd3 + beta*Sv)/(1+beta);

%%% update lbd
lbd1  = lbd1 - beta*(H(x)-v);
lbd21 = lbd21- beta*(Px(x)-y1);
lbd22 = lbd22- beta*(Py(x)-y2);
lbd3  = lbd3 - beta*(Sv-z);

PSNR(k)  = psnr(x,I);

if mod(k,5)==0;
    figure(10); imshow(x); 
    title(['Iter:',num2str(k),'  PSNR:',num2str(PSNR(k))]);  pause(0.1);
end


end


