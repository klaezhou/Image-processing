clc;  
clear;
close all;
addpath ./../Images  solvers  utilits
%%%%% target images
I = double(imread('barbarargb.png'))/255; I = I(10:10+255,420:420+255,:);
I = double(imread('brickrgb2.jpg'))/255; I = I(1:256,1:256,:);
I = double(imread('barbara.png'))/255; I = I(257:end,257:end); 


%%%% blurry
[n1,n2,n3] = size(I);
h  = fspecial('disk',3);
x0 = imfilter(I,h,'circular');

%%%% missing pixels
dtx= 64; dty = 64; S1 = rand(dtx,dty)>0.2; S2 = ones(n1/dtx,n2/dty);
S  = (kron(S1,S2)>0.5)&(imread('mask2.bmp')>200);  %%% mask 1
% S = floor(double(imread('mask3.bmp'))/255);        %%% mask2
A = zeros(n1,n2,n3); for i=1:n3; A(:,:,i) = double(S); end;   S =A;  clear A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0= S.*x0;
p = 2;
opts.MaxIt= 30;     opts.Tol = 5e-4;
opts.tau  = 0.001;  opts.mu  = 0.05; 
opts.beta1= 0.005; opts.beta2= opts.beta1; opts.beta3 = opts.beta1; 
[u,v,SNR]=InpaintBlurtau(x0,h,S,p,opts,I);


figure; 
subplot(221); imshow(I,[]);    title('clean image')
subplot(222); imshow(x0,[]);   title('observed')
subplot(223); imshow(u,[]);  title('cartoon')
subplot(224); imshow(v,[]);  title('texture')


