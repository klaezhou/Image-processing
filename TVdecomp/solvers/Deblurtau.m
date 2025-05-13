%%% min || |u| ||_1+ tau ||(u+div g)-f||^2_2 + mu|| |g| ||_p
%%%    This model is for cartoon image denoising and image decomposition 
%%% motived by
%%% [1] I. Daubechies and G.Teschke, Variational image resotration by
%%%     means of wavelets: simultaneous decomposition deblurring 
%%%     and denosing. formula (1.1)


clc;   clear;   close all;
addpath ./../Images
I = double(imread('barbarargb.png'))/255; I = I(10:10+256,420:420+256,:);
I = double(imread('barbarargb.png'))/255; I = I(1:512,250:250+511,:);
% I = double(imread('weave.jpg'))/255; 


[n1,n2,n3] = size(I);
h = fspecial('disk',11);
% h = fspecial('gaussian',21,3);
x0= imfilter(I,h,'circular') ;

p = 1; 
tau=0.00005; mu=0.00001; beta1=0.01; beta2 =0.01; beta3 = 0.01;  %p=Inf for barbara
tau=0.00005; mu=0.00001; beta1=0.001; beta2 =beta1; beta3 =beta1;  %p=Inf for weave
tau=0.00005; mu=0.00001; beta1=0.005; beta2 =beta1; beta3 =beta1;
% tau=0.0002; mu=0.00001; beta1=0.1; beta2 =0.1; beta3 = 0.1;  %p=Inf
%%%%%%%%%%%%%%%%% Periodic  boundary condtion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siz = size(h);  center = [fix(siz(1)/2+1),fix(siz(2)/2+1)];
P   = zeros(n1,n2,n3); 
for i=1:n3;  P(1:siz(1),1:siz(2),i) = h;  end
D  = fft2(circshift(P,1-center)); 
H  = @(x) (ifft2(D.*fft2(x)));         %%%% Blur operator.  B x 
HT = @(x) (ifft2(conj(D).*fft2(x)));   %%%% Transpose of blur operator.
d1h  = zeros(n1,n2,n3); d1h(1,1,:) = -1; d1h(n1,1,:) = 1; d1h = fft2(d1h);
d2h  = zeros(n1,n2,n3); d2h(1,1,:) = -1; d2h(1,n2,:) = 1; d2h = fft2(d2h);
Px = @(x) [x(2:n1,:,:)-x(1:n1-1,:,:); x(1,:,:)-x(n1,:,:)]; %%\nabla_1 x 
Py = @(x) [x(:,2:n2,:)-x(:,1:n2-1,:), x(:,1,:)-x(:,n2,:)]; %%\nabla_2 y
PTx= @(x) [x(n1,:,:)-x(1,:,:); x(1:n1-1,:,:)-x(2:n1,:,:)]; %%\nabla_1^T x 
PTy= @(x) [x(:,n2,:)-x(:,1,:), x(:,1:n2-1,:)-x(:,2:n2,:)]; %%\nabla_2^T y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% initinalization
ZERO = zeros(n1,n2,n3); 
u  = x0;    x1 = ZERO;  x2 = ZERO; y = ZERO; g1 = ZERO; g2 = ZERO; 
txtr = ZERO;
z1 = ZERO;  z2 = ZERO; lbd11 = ZERO; lbd12 = ZERO; lbd2 = ZERO; 
lbd31 = ZERO; lbd32 = ZERO;  
MDu = beta1*(abs(d1h).^2+abs(d2h).^2) + beta2;
MDy = beta2 + abs(D).^2;
Mg11= beta2*abs(d1h).^2 + beta3;
Mg22= beta2*abs(d2h).^2 + beta3;
Mg12= beta2*d1h.*conj(d2h);
Mg21= beta2*d2h.*conj(d1h);
Mdet= beta2*beta3*(abs(d1h).^2 + abs(d2h).^2) + beta3*beta3; 
HTx0= HT(x0);
Im  = norm(I(:))^2;
MaxIt = 70;
Tol = 1e-2;
COOR  = zeros(1,MaxIt);
SNR   = zeros(1,MaxIt);
tic
for k = 1:MaxIt
    
    %%%% \tilde u
    tep1 = beta1*x1 + lbd11 + beta2*g1;
    tep2 = beta1*x2 + lbd12 + beta2*g2;
    Temp = PTx(tep1) + PTy(tep2) + lbd2 + beta2*y;
    un   = real(ifft2(fft2(Temp)./MDu));
  
    %%%% \tilde z
    tep1= g1 + lbd31/beta3;
    tep2= g2 + lbd32/beta3;
    [poz1,poz2] = Project(tep1,tep2,mu/beta3,p);
    zn1 = tep1 - poz1;
    zn2 = tep2 - poz2;
 
    %%%% \tilde x
    dxun = Px(un);
    dyun = Py(un);
    sk1 = dxun - lbd11/beta1;
    sk2 = dyun - lbd12/beta1;
    nsk = sqrt(sk1.^2 + sk2.^2); nsk(nsk==0)=1;
    nsk = max(1 - (tau/beta1)./nsk,0);
    xn1 = nsk.*sk1;
    xn2 = nsk.*sk2;
    
    %%%% \tilde g
    Temp = beta2*(un-y) - lbd2;
    b1 = Px(Temp) + beta3*zn1 - lbd31;
    b2 = Py(Temp) + beta3*zn2 - lbd32;
    nume1 = Mg22.*fft2(b1) - Mg12.*fft2(b2);
    nume2 = Mg11.*fft2(b2) - Mg21.*fft2(b1);    
    gn1 = real( ifft2(nume1./Mdet) );
    gn2 = real( ifft2(nume2./Mdet) );
    
    %%%% \tilde y
    txtrpre=txtr;
    txtr= -PTx(gn1)-PTy(gn2);
    Temp= HTx0 + beta2*(un+txtr) - lbd2;
    yn  = real(ifft2(fft2(Temp)./MDy));
    
    %%%% update lagrange multipliers
    lbdn11= lbd11 - beta1*(dxun-xn1);
    lbdn12= lbd12 - beta1*(dyun-xn2);
    lbdn2 = lbd2  - beta2*(un+txtr-yn);
    lbdn31= lbd31 - beta3*(zn1-gn1);
    lbdn32= lbd32 - beta3*(zn2-gn2);
    
    
    stopic = max(norm(u(:)-un(:))/norm(u(:)),norm(txtr(:)-txtrpre(:))/norm(txtr(:)))
    if stopic<Tol
        break;
    end
        
    u = un; g1 = gn1; g2=gn2; x1=xn1; x2=xn2; y=yn; z1=zn1; z2=zn2;
    lbd11 = lbdn11; lbd12 = lbdn12; lbd2 = lbdn2; lbd31 = lbdn31; 
    lbd32 = lbdn32;
  
    COOR(k)=corr(u(:),txtr(:));
    SNR(k)=10*log10(Im/norm(u(:)+txtr(:)-I(:))^2);
    
        
end
toc

figure; 
subplot(221); imshow(I,[]);    title('original')
subplot(222); imshow(x0,[]);   title('blure+missing')
subplot(223); imshow(u,[]);    title('re-cartoon')
subplot(224); imshow(txtr+0.5,[]); title('re-texture')
% imshow(u,'border','tigh'),rectangle('Position', [180 10 70 70],'EdgeColor','w','lineWidth',3)
% imshow(u(10:10+70,180:180+70,:),'border','tigh')
% imshow(u,'border','tigh'),rectangle('Position', [360 5 100 128],'EdgeColor','r','lineWidth',3)
a=[u(5:5+128,360:360+100,:),ones(129,5,3),0.5+txtr(5:5+128,360:360+100,:)];
imshow(a,'border','tigh')