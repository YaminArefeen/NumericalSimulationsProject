% copyright Luca Daniel, MIT 2018
% modified by Zijing 2018/10/21
clear;
close all;
addpath funcs;
load('MRI_tumor_image.mat');
Nx=64;   % larger matrix size shows higher resolution
Ny=64;
%% down-sampling the image to lower resolution
T1_image = real(ifft2c(crop(fft2c(double(T1_image)),[Nx,Ny])));
Tumor_mask = real(ifft2c(crop(fft2c(double(Tumor_mask)),[Nx,Ny])));
brainmask = real(ifft2c(crop(fft2c(double(brainmask)),[Nx,Ny])));

Tumor_mask=Tumor_mask>5*mean(Tumor_mask(:));
brainmask=brainmask>1*mean(brainmask(:));
figure; 
subplot(1,3,1); imshow(T1_image,[]); title('MR Image');
subplot(1,3,2); imshow(Tumor_mask,[]); title('Tumor Mask');
subplot(1,3,3); imshow(brainmask,[]); title('Brain Mask');

mask=zeros(size(T1_image));
mask(Tumor_mask)=1;              % 1 indicates tumor, 0 is air
mask(logical(brainmask-Tumor_mask))=2;    % 2 indicates health tissues
location_tumor = [49,41];
%select input evaluation functions
eval_u = 'eval_u_step';
% 2D Heat Process

[p,x_start,t_start,t_stop,max_dt_FE] = getParam_2DHeating_Brain(1,1,Nx,Ny,mask,location_tumor);
eval_f = 'eval_f_LinearSystem';
figure; imagesc(p.A);
p.Nx=Nx;
p.Ny=Ny;
% test FE function
timestep=0.001;
t_stop=10;
max_temp=10;
visualize=1;
step_interval_show=50;   % the time interval of image update ( 1 updata/ N steps )
[X] = ForwardEuler2D_Brain_POD(eval_f,x_start,eval_u,p,t_start,t_stop,timestep,T1_image,brainmask,max_temp,step_interval_show,visualize);
%% SVD model reduction
X_ref=X;
p.c=ones(size(p.B,1),1);
p.c=diag(p.c);

N_sample=1000;
idx = randperm(size(X,2));
X0=X(:,idx(1:N_sample));

figure;
plot(X0, 'linewidth', 1);
axis('tight');
xlabel('# nodes');
ylabel('T');
title('Simulated signal space');
%%
[U, S, ~] = svd(X0, 'econ');
S=diag(S);
figure;
plot(S, 'linewidth', 2);
xlabel('# nodes');
ylabel('T');
title('Singular Values');

figure;
plot(S(1:10), 'linewidth', 2);
xlabel('# nodes');
ylabel('T');
title('Largest 10 Singular Values');
axis tight;
sum(S(1:5))/sum(S)
%% q=2
A=p.A;
b=p.B;
c=p.c;
q=5;
V = U(:,1:q);
A1=V'*A*V;
b1= V'*b;
c1= V'*c;
q=10;
V = U(:,1:q);
A2=V'*A*V;
b2= V'*b;
c2= V'*c;
%% q=5
visualize=0;
p.A=A1;
p.B=b1;
x_start=zeros(size(A1,1),1);
tic;
[X] = ForwardEuler2D_Brain_POD(eval_f,x_start,eval_u,p,t_start,t_stop,timestep,T1_image,brainmask,max_temp,step_interval_show,visualize);
toc;
y1=c1'*X;
norm(X_ref(:)-y1(:),1)
Temperature_im=zeros(Nx,Ny);
Temperature_im(brainmask)=y1(:,end);
Temperature_ref=zeros(Nx,Ny);
Temperature_ref(brainmask)=X_ref(:,end);
figure; imshow3(cat(3,Temperature_im,Temperature_ref,100*abs(Temperature_ref-Temperature_im)),[],[1 3]); colormap('hot');
figure; imshow3(cat(3,Temperature_im,Temperature_ref,1000*abs(Temperature_ref-Temperature_im)),[],[1 3]); colormap('hot');
%% q=10
visualize=0;
p.A=A2;
p.B=b2;
x_start=zeros(size(A2,1),1);
tic;
[X] = ForwardEuler2D_Brain_POD(eval_f,x_start,eval_u,p,t_start,t_stop,timestep,T1_image,brainmask,max_temp,step_interval_show,visualize);
toc;
y2=c2'*X;
norm(X_ref(:)-y2(:),1)
Temperature_im=zeros(Nx,Ny);
Temperature_im(brainmask)=y2(:,end);
Temperature_ref=zeros(Nx,Ny);
Temperature_ref(brainmask)=X_ref(:,end);
figure; imshow3(cat(3,Temperature_im,Temperature_ref,1000*abs(Temperature_ref-Temperature_im)),[],[1 3]); colormap('hot');
figure; imshow3(cat(3,Temperature_im,Temperature_ref,10000*abs(Temperature_ref-Temperature_im)),[],[1 3]); colormap('hot');
norm(X_ref(:)-y2(:),1)
