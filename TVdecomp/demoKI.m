clc;   
clear; 
close all;
addpath ./../Images  solvers  utilits
%I = double(imread('barbara.png'))/255;  
%I = double(imread('weave.jpg'))/255;
cart = double(imread('TomAndJerry.png'))/255; 
text = double(imread('wool.png'))/255; 
I = 0.7*cart+0.3*text;  
Tol        = 0.04;
opts.MaxIt = 200; 
opts.beta  = 1; 
opts.tau = 0.01;
opts.mu  = 1;
[u_G,v_G] = TVdecom(I,opts,Tol);
%[u_G,v_G] = TVdecomNew(I,opts,Tol);
