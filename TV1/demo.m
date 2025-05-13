clc; 
clear;
close all;
clc
clear
close all
addpath ../Images                     %% 加入图片所在路径\文件夹

I  = im2double(imread('cameraman.tif'));  %%%读取图片
h  = fspecial('average',7);               %%%模糊的图片
xh = imfilter(I,h,'circular'); 
x0 = imnoise(xh,'salt & pepper',0.5);     %%%impulse噪声

opts.MaxIt = 100;
opts.alpha = sum(abs(x0(:)-xh(:))); 
opts.beta1 = 10;
opts.beta2 = 10;

[x,PSNR]=TVL1(x0,h,opts,I);

figure; plot(0:opts.MaxIt-1,PSNR,'-.+','LineWidth',2,'MarkerSize',8);
legend('ADMM'); 
xlabel('Iteration   No.','FontSize',14); 
ylabel('SNR  (dB)','FontSize',14); 
 

