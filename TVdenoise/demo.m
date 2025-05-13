clc; 
clear;
close all;
addpath ../Images                     %% 加入图片所在路径\文件夹
%I = im2double(imread('shape.jpg'));   %% 读图片
I = im2double(imread('chart.tiff'));
%I = im2double(imread('housergb.png'));

[n1,n2,n3] = size(I);            %% 获取图片大小
x0 = I+0.2*randn(n1,n2,n3);      %% 加高斯噪声
%%%%%%%%%%%% ADMM 求解优化问题  %%%%%%%%%%%%%%%%%%%%%%%%%%
opts.beta=1;      %%% ADMM算法的罚参数
opts.mu=5;        %%% 模型中的参数
opts.MaxIt=100;   %%% 迭代次数
opts.Tol=1e-4;       %%% 停机准则
[u,PSNR,Time] = TV_denoise(x0,opts,I);  %%% 主程序


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;   %%%%% 显示求解结果
subplot(221); imshow(I,[]);   title('ideal image') 
subplot(222); imshow(x0,[]);  title('degraded image')
subplot(223); imshow(u,[]);   title('denoised image')

figure;  %%%绘制PSRN 随迭代变化曲线
plot(0:opts.MaxIt-1,PSNR,'-','LineWidth',4); 
legend('ADM')
xlabel('Iteration No.','Fontsize',16);
ylabel('PSNR (dB)','Fontsize',16)

