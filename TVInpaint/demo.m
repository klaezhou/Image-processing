clc;  
% clear;
close all;
clc
clear 
close all
I = im2double(imread('cameraman.tif')); 
% I = im2double(imread('housergb.png')); 
% I = im2double(imread('barbara.png')); I = I(257:end,257:end);

[n1,n2,n3] = size(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = floor(double(imread('mask3.bmp'))/255);      %%%%% Mask 1
A = zeros(n1,n2,n3); for i=1:n3; A(:,:,i) = S; end;   S =A;  clear A
%%%%%% Blur kernel
h  = fspecial('aver',3);                     %%%%%% blur + noise
xh = imfilter(I,h,'circular');
xb = S.*xh;                                  %%%%%% missing pixels     

opts.MaxIt= 200; 
opts.tau  = 0.01;    
opts.beta = 1; 
[x,PSNR]  = TVinpaint(xb,h,S,opts,I); 

%%%%%%%%%%%%%%%%%%%%显示结果及曲线
figure; 
plot(0:opts.MaxIt-1,PSNR,'-o','LineWidth',2,'MarkerSize',8);
legend('ADMM');
xlabel('Iteration   No.', 'FontSize',16); 
ylabel('PSNR  (dB)', 'FontSize',16);

figure
subplot(221); imshow(I);  title('原始图片')
subplot(222); imshow(xb); title('观测图片')
subplot(223); imshow(x);title('恢复结果')
