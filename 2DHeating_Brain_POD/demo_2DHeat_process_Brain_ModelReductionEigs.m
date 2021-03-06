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
visualize=0;
step_interval_show=50;   % the time interval of image update ( 1 updata/ N steps )
tic
[X] = ForwardEuler2D_Brain_POD(eval_f,x_start,eval_u,p,t_start,t_stop,timestep,T1_image,brainmask,max_temp,step_interval_show,visualize);
toc
X_ref = X;
%% Eigenvalue Mode Reduction

% Get the full system parameters
A = p.A;
b = p.B;
p.c=ones(size(p.B,1),1);
p.c=diag(p.c);
c = p.c;

% Select the reduction
% q=size(A);
q=102;
type_trunc=4; %choose from residue:0, combo:1, slow:2, fast:3, halfhalf:4

% Determine the reduced system parameters
[rA, rb, rc] = getReducedEigSystem(A,b,c,q,type_trunc);

%%

% Solve for the reduced parameters
visualize=0;
timestep=0.001;
p.A=rA;
p.B=rb;
x_start=zeros(size(rA,1),1);
tic;
[X] = ForwardEuler2D_Brain_POD(eval_f,x_start,eval_u,p,t_start,t_stop,timestep,T1_image,brainmask,max_temp,step_interval_show,visualize);
toc;
y1=rc'*X;
norm(X_ref(:,end)-y1(:,end),1)

%%
Temperature_im=zeros(Nx,Ny);
Temperature_im(brainmask)=y1(:,end);
Temperature_ref=zeros(Nx,Ny);
Temperature_ref(brainmask)=X_ref(:,end);
figure; imshow3(cat(3,Temperature_im,Temperature_ref,100*abs(Temperature_ref-Temperature_im)),[],[1 3]); colormap('hot');
figure; imshow3(cat(3,Temperature_im,Temperature_ref,1000*abs(Temperature_ref-Temperature_im)),[],[1 3]); colormap('hot');
