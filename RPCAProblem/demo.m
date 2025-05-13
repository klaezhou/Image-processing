clc
clear
close all
randn('state',0);   
rand('twister',0);
load Hall_airport_1000_1496_497_144_176_gray.mat
b = images(:,1:200) ; clear images
[m,n] = size(b) ;
sr    = 0.8;                        %%% Incomplete info.
Omega = randperm(m*n); p = round(sr*m*n); 
Omega = Omega(1:p);    Omega = Omega';  b = b(Omega);            
sigma = 0.001;         b = b + sigma * randn(p,1); %%%% add noise
fprintf('Vedio_Missing data:m:%3d,n:%3d,Omega:%4.2f **\n',m,n,sr);
opts.tau = 1/sqrt(m);            %% model parameter
opts.mu  = 0.01;
opts.beta  = 0.005/mean(abs(b)); % penalty
opts.MaxIt = 50;

[M,Y,X] = EADM3([m,n],Omega,b,opts); 
%
for i=1:200
    figure(1);
    I=reshape(M(:,i),144,176); subplot(221); imshow(I,[]); title('原始');pause(0.01);
    I=reshape(Y(:,i),144,176); subplot(222); imshow(I,[]); title('运动');pause(0.01);
    I=reshape(X(:,i),144,176); subplot(223); imshow(I,[]); title('静止');pause(0.01);
end
close

%
% for i=1:200
%     figure(2);
%     I=reshape(b(:,i),144,176); imshow(I,[]); title('监控视频');pause(0.01);
% end


% %%
% aviobj = avifile('ext_video.avi');
% for i=1:200
%     h = imshow(reshape(M(:,i),144,176),[]);
% %     I=[reshape(X(:,i),144,176),ones(144,3)*255,reshape(Y(:,i),144,176)];
% %     h = imshow(I,[],'border','tigh');colormap(gray);
%     frame = getframe(gca);
%     aviobj= addframe(aviobj,frame);
% end
% aviobj = close(aviobj);

