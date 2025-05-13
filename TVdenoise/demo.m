clc; 
clear;
close all;
addpath ../Images                     %% ����ͼƬ����·��\�ļ���
%I = im2double(imread('shape.jpg'));   %% ��ͼƬ
I = im2double(imread('chart.tiff'));
%I = im2double(imread('housergb.png'));

[n1,n2,n3] = size(I);            %% ��ȡͼƬ��С
x0 = I+0.2*randn(n1,n2,n3);      %% �Ӹ�˹����
%%%%%%%%%%%% ADMM ����Ż�����  %%%%%%%%%%%%%%%%%%%%%%%%%%
opts.beta=1;      %%% ADMM�㷨�ķ�����
opts.mu=5;        %%% ģ���еĲ���
opts.MaxIt=100;   %%% ��������
opts.Tol=1e-4;       %%% ͣ��׼��
[u,PSNR,Time] = TV_denoise(x0,opts,I);  %%% ������


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;   %%%%% ��ʾ�����
subplot(221); imshow(I,[]);   title('ideal image') 
subplot(222); imshow(x0,[]);  title('degraded image')
subplot(223); imshow(u,[]);   title('denoised image')

figure;  %%%����PSRN ������仯����
plot(0:opts.MaxIt-1,PSNR,'-','LineWidth',4); 
legend('ADM')
xlabel('Iteration No.','Fontsize',16);
ylabel('PSNR (dB)','Fontsize',16)

