


% This function computes:
% argmin_x |||Dx|||_1, ||SHx-x^0||_2<=alpha
% L(x,y,lbd) = ||w||_1 + chi_K(z) + <lbd, Ax-y> + beta/2 ||Ax-y||^2_2
% Using an ADM approach. Where s is a downsampling operator and H is a
% convolution operator (circular boundary conditions).
% 1.in the case of l_2 norm, the constraint is ||SHx-x^0||_2<=alpha
% 2. It is clear that TV is not a serious choice for image SR. 
% IN:
% - x0    : initial data.
% - MaxIt : iteration number.
% - H     : convolution operator (note: give it in spatial domain).
% - alpha : norm of the noise.
% - s     : upsampling factor (integer).
% OUT: 
% - x     : approximation of the solution.

function x = TVzoom(x0,h,opts)


%% Parameters 
alpha= opts.alpha;
beta = opts.beta;
s    = opts.s;
MaxIt= opts.MaxIt;

[n1,n2,n3] = size(x0);
n1  = n1*s; 
n2  = n2*s;

% %%%%%%%%%%%周期边界条件：fourier域下的filter 
siz = size(h);
center = [fix(siz(1)/2+1),fix(siz(2)/2+1)];
P  = zeros(n1,n2,n3); 
for i =1:n3; P(1:siz(1),1:siz(2),i) = h; end
D  = fft2(circshift(P,1-center));
H  = @(x) real(ifft2(D.*fft2(x)));       %%%% B x 
HT = @(x) real(ifft2(conj(D).*fft2(x))); %%%% B^Tx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 图像的梯度 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1h = zeros(n1,n2,n3); d1h(1,1,:) = -1; d1h(n1,1,:) = 1; d1h = fft2(d1h);
d2h = zeros(n1,n2,n3); d2h(1,1,:) = -1; d2h(1,n2,:) = 1; d2h = fft2(d2h);
Px  = @(x) [x(2:n1,:,:)-x(1:n1-1,:,:); x(1,:,:)-x(n1,:,:)]; %%\nabla_1 x 
Py  = @(x) [x(:,2:n2,:)-x(:,1:n2-1,:), x(:,1,:)-x(:,n2,:)]; %%\nabla_2 y
PTx = @(x) [x(n1,:,:)-x(1,:,:); x(1:n1-1,:,:)-x(2:n1,:,:)]; %%\nabla_1^T x 
PTy = @(x) [x(:,n2,:)-x(:,1,:), x(:,1:n2-1,:)-x(:,2:n2,:)]; %%\nabla_2^T y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

%% Initializations y=(w,z), x and lbd.
x  = x0;
z  = zeros(n1,n2,n3);              
w1 = zeros(n1,n2,n3);
w2 = zeros(n1,n2,n3);
lbd1  = zeros(n1,n2,n3);
lbd21 = zeros(n1,n2,n3);
lbd22 = zeros(n1,n2,n3);
MatD  = abs(d1h).^2 + abs(d2h).^2 + abs(D).^2;

for i = 1:MaxIt
    
    %%% x^{k+1}
    Temp = HT(z+lbd1/beta)+PTx(w1+lbd21/beta)+PTy(w2+lbd22/beta);
    x = real(ifft2(fft2(Temp)./MatD));  
    
    
    %%% y^{k+1}
    %Projection of z onto the l^2-ball
    Hx = H(x);
    z  = Hx - lbd1/beta;
    z  = ProjZ(z,x0,s,alpha);    
    
    % shrinkage on w
    sk1 = Px(x) - lbd21/beta;
    sk2 = Py(x) - lbd22/beta;
    nsk = sqrt(sk1.^2 + sk2.^2); nsk(nsk==0) = 1;
    nsk = max(1 - 1./(beta*nsk),0);
    w1  = sk1.*nsk;
    w2  = sk2.*nsk;
    
    %%% update lbd^{k+1}
    lbd1  = lbd1 + beta*(z-Hx);
    lbd21 = lbd21+ beta*(w1-Px(x));
    lbd22 = lbd22+ beta*(w2-Py(x));
    
    %%%%%%% 显示结果
    if mod(i,10)==0
        figure(10);imshow(x,[]); title(['Iter: ',num2str(i)]); pause(0.1)
    end    

end



function x = ProjZ(x,x0,s,alpha)   %%%% 在L2-ball上的投影
xp = x(1:s:end,1:s:end,:);
dst= norm(x0(:)-xp(:));
    if dst>alpha
        xp = x0+(alpha/dst)*(xp-x0);
    end
x(1:s:end,1:s:end,:) = xp ;
















% S= [ 1 0 0 0 0 0 0 0 0 0 0 0 0;
%      0 0 0 1 0 0 0 0 0 0 0 0 0;
%      0 0 0 0 0 0 1 0 0 0 0 0 0;
%      0 0 0 0 0 0 0 0 0 1 0 0 0;
%      0 0 0 0 0 0 0 0 0 0 0 0 1;
%       ];
% a =1:13; a=a';
% S*a
