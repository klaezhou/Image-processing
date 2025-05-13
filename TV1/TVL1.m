function [x,PSNR]=TVL1(x0,h,opts,I)

[n1,n2,n3] = size(x0);

alpha = opts.alpha;   %%% 模型参数
beta1 = opts.beta1;  
beta2 = opts.beta2;  %%ADMM方法罚参数
MaxIt = opts.MaxIt;  %%%迭代次数

%%%%%%%%%%%周期边界条件：fourier域下的filter 
siz = size(h);
center = [fix(siz(1)/2+1),fix(siz(2)/2+1)];
P  = zeros(n1,n2,n3); 
for i =1:n3; P(1:siz(1),1:siz(2),i) = h; end
D  = fft2(circshift(P,1-center));
H  = @(x) real(ifft2(D.*fft2(x)));       %%%% Bx 
HT = @(x) real(ifft2(conj(D).*fft2(x))); %%%% B^Tx

%%%%%%%%%%%%%%%%% 图像的梯度 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1h = zeros(n1,n2,n3); d1h(1,1,:) = -1; d1h(n1,1,:) = 1; d1h = fft2(d1h);
d2h = zeros(n1,n2,n3); d2h(1,1,:) = -1; d2h(1,n2,:) = 1; d2h = fft2(d2h);
Px  = @(x) [x(2:n1,:,:)-x(1:n1-1,:,:); x(1,:,:)-x(n1,:,:)]; %%\nabla_1 x 
Py  = @(x) [x(:,2:n2,:)-x(:,1:n2-1,:), x(:,1,:)-x(:,n2,:)]; %%\nabla_2 y
PTx = @(x) [x(n1,:,:)-x(1,:,:); x(1:n1-1,:,:)-x(2:n1,:,:)]; %%\nabla_1^T x 
PTy = @(x) [x(:,n2,:)-x(:,1,:), x(:,1:n2-1,:)-x(:,2:n2,:)]; %%\nabla_2^T y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Initializations y=(w,z), x and lbd.
x    = x0;
yw1  = Px(x0); yw2=Py(x0); yz = H(x); 
lbdw1= zeros(n1,n2,n3); 
lbdw2= zeros(n1,n2,n3); 
lbdz = zeros(n1,n2,n3); 

MatD = beta1*(abs(d1h).^2+abs(d2h).^2) + beta2*abs(D).^2;  %%%系数矩阵

PSNR = zeros(MaxIt,1); 
for k= 1:MaxIt
  %%%  x^{k+1} 
  Temp = PTx(beta1*yw1-lbdw1)+PTy(beta1*yw2-lbdw2)+HT(beta2*yz-lbdz);
  x    = ifft2(fft2(Temp)./MatD); 
  
  %%% y^{k+1}  
  sk1 = Px(x) + lbdw1/beta2;
  sk2 = Py(x) + lbdw2/beta2;
  nsk = sqrt(sk1.^2+sk2.^2); nsk(nsk==0)=1;
  nsk = max(1-1./(beta2*nsk),0);
  yw1 = nsk.*sk1;
  yw2 = nsk.*sk2;  
  
  %%%%% z onto the L_1-ball
  yz = H(x) + lbdz/beta2 - x0; 
  yz = x0 + ProjectionWL1(yz,alpha);  
 
  %%% lbd^{k+1}
  lbdw1 = lbdw1 + beta1*(Px(x)-yw1);
  lbdw2 = lbdw2 + beta1*(Py(x)-yw2);  
  lbdz  = lbdz  + beta2*(H(x)-yz);  
  
  %%%%%%%%%%%%%%%%%%%
  PSNR(k) = psnr(x,I);   %%%计算信噪比
    
  figure(10);  imshow(x); 
  title(['Iter:',num2str(k),'  PSNR:',num2str(PSNR(k))]);  pause(0.1);
  
  
end


function [z,sigma]=ProjectionWL1(xx,alpha)
x = xx(:);
ax= abs(x(:));
n = numel(x);
M = sum(ax);
if M<=alpha;  z=xx; sigma=0;  return; end
if alpha<=0
    z=zeros(size(xx));  sigma=infty;   return;
end
y=sort(ax); E1=M; E=E1-n*y(1); i=1;
while E>alpha && i<n
    i = i+1;
    E1= E1-y(i-1);
    Ep= E;
    E = E1-(n-i+1)*y(i);
end
if (i>1)
  a=(y(i)-y(i-1))*alpha;
  b=E*y(i-1)-abs(Ep)*y(i);
  sigma=(a+b)/(E-abs(Ep));
else
  sigma=y(1)*(M-alpha)/(M-E);    
end
z=x(:); K = find (ax<sigma, 1);
if (~isempty(K));  z(ax<sigma)=0; end
K=find(ax>=sigma);
if (~isempty(K));  z(K)=x(K)-sigma*sign(x(K)); end
z=reshape(z,size(xx));




