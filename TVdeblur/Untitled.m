clc
clear all
close all
I=im2double(imread('circles.tif')); I=I(1:4:end,1:4:end);
h = fspecial('average',9); I=conv2(I,h);
[n1,n2]=size(I);



% Px = @(x) [x(2:n1,:,:)-x(1:n1-1,:,:); x(1,:,:)-x(n1,:,:)]; %%% partial x 
% Py = @(x) [x(:,2:n2,:)-x(:,1:n2-1,:), x(:,1,:)-x(:,n2,:)]; %%% partial y
Px = @(x) [x(2:n1,:)-x(1:n1-1,:); zeros(1,n2)]; %%% partial x 
Py = @(x) [x(:,2:n2)-x(:,1:n2-1), zeros(n1,1)]; %%% partial y
Ix=Px(I); Iy=Py(I);
M  = Ix.^2+Iy.^2;
figure; imagesc(I); hold on; colormap gray; axis equal
[X,Y]=meshgrid(1:n1,1:n2);
quiver(X,Y,Iy,Ix,'r')
