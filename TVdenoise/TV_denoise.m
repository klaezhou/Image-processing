

function [u,PSNR,Time,Itr] = TV_denoise(x0,opts,I)
%%%% This function solve the unconstrained model for image deblurring
%%%% min_u |\nabla u|_1 + 0.5 |u-u0|_2^2
%%%%  u -- original image
%%%%  u0 -- observed image
%%%% Alternating direction method of multiplier is applied to sove the above model


[n1,n2,n3] = size(x0);  %% 图片大小
beta = opts.beta;       %% ADMM罚参数
mu   = opts.mu;         %% 模型参数
MaxIt= opts.MaxIt;      %% 最大迭代次数

%%%%%%%%%%%%%%%%% 周期边界条件 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1h = zeros(n1,n2,n3); d1h(1,1,:) = -1; d1h(n1,1,:) = 1; d1h = fft2(d1h);
d2h = zeros(n1,n2,n3); d2h(1,1,:) = -1; d2h(1,n2,:) = 1; d2h = fft2(d2h);
Px  = @(x) [x(2:n1,:,:)-x(1:n1-1,:,:); x(1,:,:)-x(n1,:,:)]; %%\nabla_1 x 
Py  = @(x) [x(:,2:n2,:)-x(:,1:n2-1,:), x(:,1,:)-x(:,n2,:)]; %%\nabla_2 y
PTx = @(x) [x(n1,:,:)-x(1,:,:); x(1:n1-1,:,:)-x(2:n1,:,:)]; %%\nabla_1^T x 
PTy = @(x) [x(:,n2,:)-x(:,1,:), x(:,1:n2-1,:)-x(:,2:n2,:)]; %%\nabla_2^T y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Intionalization
u    = x0; 
v1   = Px(u);
v2   = Py(u);
lbd11= zeros(n1,n2,n3);
lbd12= zeros(n1,n2,n3);
MDu  = mu + beta*(abs(d1h).^2+abs(d2h).^2);  %%%线性方程组
HTx0 = mu*x0;
PSNR = zeros(1,MaxIt); 

Time = cputime;  %%%记录时间
k    = 1; 
while k <= MaxIt 
    
    %%% step 1:  u -- 子问题
    Temp= PTx(beta*v1+lbd11) + PTy(beta*v2+lbd12) + HTx0;
    un  = real(ifft2(fft2(Temp)./MDu));    
    %%% step 2:  v --子问题
    dxun= Px(un);
    dyun= Py(un);
    sk1 = dxun - lbd11/beta;
    sk2 = dyun - lbd12/beta;
    nsk = sqrt(sk1.^2 + sk2.^2); nsk(nsk==0)=1;
    nsk = max(1-1./(beta*nsk),0);
    vn1 = sk1.*nsk;
    vn2 = sk2.*nsk;         
    %%%% 更新\lambda
    lbdn11 = lbd11 - beta*(dxun - vn1);
    lbdn12 = lbd12 - beta*(dyun - vn2);    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    PSNR(k) = psnr(u,I);  %%% 计算信噪比
    stopic = norm(u(:)-un(:))/norm(un(:));
    Time = cputime - Time;
    
    figure(10);
    imshow(u); 
    title(['Iter:',num2str(k),   'PSNR:',num2str(PSNR(k))])
    pause(0.1);
    u = un; 
    v1 =vn1; v2 = vn2; 
    lbd11 = lbdn11; lbd12 =lbdn12;
    k = k+1;
end
Time = cputime-Time; Itr =k;

    

