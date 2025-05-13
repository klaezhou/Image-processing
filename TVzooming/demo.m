clc; 
clear;
close all
I = im2double(imread('cameraman.tif'));
[n1,n2,n3] = size(I);

h = fspecial('gaussian',7,2);                       %%filter B
xh = imfilter(I,h,'circular');  
s  = 4;                         %% downsample facor
xh = xh(1:s:n1,1:s:n2,:);       %% low-resolution
b  = 0.01*randn(n1/s,n2/s,n3);  %% Noise
xb = xh + b;                    %% noisy low-resolution

opts.alpha = 1;
opts.beta  = 100;
opts.s     = s;
opts.MaxIt = 200;
x = TVzoom(xb,h,opts);

%%%%% display results
figure; 
subplot(221); imagesc(I);  title('原始图片')
subplot(222); imagesc(xb); title('低精度')
subplot(223); imagesc(x);  title('恢复结果'); colormap gray















% x = 1:2:20.5; y =0.1*sin(0.3*x);  %%%1-demension
% subplot(211);plot(x,y,'o','MarkerSize',8);  axis([1 20 -0.1 0.1])
% xx=1:1:20.5; yy=interp(y,2);
% subplot(212);plot(xx,yy,'ro',x,y,'o','MarkerSize',8); axis([1 20 -0.1 0.1])
% 
% 
% [x,y,z] = peaks(15);  figure;subplot(221); mesh(x,y,z);   %%SR
% [xi,yi] = meshgrid(-3:.1:3,-3:.1:3);
% zi = interp2(x,y,z,xi,yi,'nearest'); subplot(222);mesh(xi,yi,zi);
% zi = interp2(x,y,z,xi,yi,'linear');  subplot(223);mesh(xi,yi,zi);
% [x3,y3,z3] = peaks(60);subplot(224); mesh(x3,y3,z3);pause
 
% I = double(imread('lighthouse2.png')); [n1,n2,n3]=size(I);s=4;
% Iz= zeros(n1*s,n2*s,n3); Iz(1:s:end,1:s:end,:)=I;
% Iz=NNInterp(I,4); imshow(Iz,[]);    %%% interpolation method

% % 
% I = double(imread('house.png'))/256; [n1,n2,n3]=size(I);
% s1=4;s2=4;
% Iz=I(1:s1:end,1:s2:end,:);
% Iz=NNInterp(I,4); imshow(Iz,[]);    %%% interpolation method


% h = [1,2,1;
%      2,4,2;
%      1,2,1];   %%% convolution method
% h = h/sum(abs(h(:)));  
% siz = size(h); center = [fix(siz(1)/2+1),fix(siz(2)/2+1)];
% P = zeros(n1*s,n2*s); 
% P(1:siz(1),1:siz(2)) = h;
% D = fft2(circshift(P,1-center));
% Iz=ifft2(D.*fft2(Iz)); imshow(Iz,[])

% a  =1:15;          %%%% 1-dimension
% da =downsample(a,3);
% ua =upsample(da,3);

% a  =reshape(1:8*8,8,8); %%%% 2-dimension
% da =a(1:2:8,1:2:8);   
% ua =zeros(8); ua(1:2:8,1:2:8)=da;

% a  =reshape(1:8*8,8,8);s=2;%%%% 2-dimension (advance)
% D1 = kron(speye(8/s),ones(1,s));   
% D2 = kron(speye(8/s),ones(1,s));
% R = kron(D1,D2)/s^2;   
% da=reshape(R*a(:),4,4);
% % ua=reshape(s^2*R'*da(:),8,8)