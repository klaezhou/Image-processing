clc
clear
close all
addpath ../Images
I = im2double(imread('shape.jpg'));    %%%测试图片
I = im2double(imread('chart.tiff'));
 I = im2double(imread('housergb.png'));


[n1,n2,n3] = size(I);
h  = fspecial('aver',5);    %%%% filter
x0 = imfilter(I,h,'circular')+ 0.20*randn(n1,n2,n3);  %%%加blur+noise

%%%%%%%%%%%%%%%%%%% ADMM 方法求解 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.beta = 1; %%%罚参数
opts.mu   = 10; %%%模型参数
opts.MaxIt= 150;  %%%最大迭代次数
opts.Tol  = 1e-5; 
[u, PSNR, Time, Itr]  = TV_deblur(x0,h,opts,I);      %


%%%%%%%%%%%%%显示结果
figure; 
subplot(221); imshow(I,[]);  title('原始图片')
subplot(222); imshow(x0,[]); title('模糊图片')
subplot(223); imshow(u,[]); title('恢复后图片')

figure; a = 0:opts.MaxIt-1; b = a+1;  %%绘制PSNR曲线
plot(0:opts.MaxIt-1,PSNR,'LineWidth',4);
legend('ADM')
xlabel('Iteration No.','Fontsize',16);
ylabel('SNR (dB)','Fontsize',16)


