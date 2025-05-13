clc;  
clear; 
close all;
addpath ./../Images  solvers  utilits

%%%%%% target image
I = double(imread('barbara.png'))/255; I = I(257:end,257:end);

%%%%% mask
[n1,n2,n3] = size(I);
S = floor(double(imread('mask3.bmp'))/255);    %%% mask
A = zeros(n1,n2,n3); for i=1:n3; A(:,:,i) = S; end;   S =A;  clear A
x0= S.*I;

p = inf;  opts.MaxIt= 50; 
opts.tau  = 0.01;  opts.mu   = 0.01;
opts.beta1= 0.1;   opts.beta2= opts.beta1;   opts.beta3 = opts.beta1; 
[u_S,v_S,SNR_S]=Inpainttau(x0,S,p,opts,I);



figure; 
subplot(221); imshow(I,[]);    title('clean image')
subplot(222); imshow(x0,[]);   title('observed')
subplot(223); imshow(u_S,[]);  title('cartoon')
subplot(224); imshow(v_S,[]);  title('texture')
