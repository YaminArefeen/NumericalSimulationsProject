%This script demonstrates the ability of a reduced system to speed up the solver for
%acceptable accuracy

%Note that this is ONLY for the linear system solve; not the nonlinear
%system!

%% 0) Initializions
% clear;
clear all;
close all;
addpath funcs;
load('MRI_tumor_image.mat');
Nx=64;   % larger matrix size shows higher resolution
Ny=64;

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
%select input evaluation functions
eval_u = 'eval_u_step';

figure; imshow(mask,[]);
load('mask.mat')

location_tumor = [47,41];
[p] = getParam_2DHeating_Brain_Jacobian(1,1,Nx,Ny,mask,location_tumor);

sigma = 1e-1; %scaling for our nonlinear term.
F_jac = @(T,t) evalSysandJac(T,p,sigma);
%forward function which evaluates our system and Jaobian
F_nojac = @(T,t) eval_f_LinearSystem(T,p,sigma);
%forward function which evaluates just our system.

%Computing reference steady state solution
tol = 1e-8;
maxiters = 50;
x_ref = NewtonSolve(F_jac,zeros(size(p.A,1),1),p,sigma,tol,maxiters);

%The initial condition
x_start = zeros(size(p.A,1),1); 

%% 2) Demonstrate the TRAP solver for the entire system
dt_trap = 1;
t_stop = 300;
tolerances_newt = [1e-10,1e-10,1e0];

tic
trap_x = trap(x_start,dt_trap,t_stop,F_jac,tolerances_newt,x_ref);
trap_time = toc;

accuracy_trap = norm(x_ref - trap_x(:,end))/norm(x_ref);

fprintf('Accuracy Trap, %f, achieved in %f seconds\n'...
, accuracy_trap,trap_time)

figure
plot(trap_x(:,end))
hold on
plot(x_ref)

%% 3) Demonstrate the TRAP solver for the Eig reduced system
dt_trap=0.5;

%Calculate the reduced system
p.c=ones(size(p.B,1),1);
p.c=diag(p.c);
q=1000;
type_trunc=0; %choose from residue:0, combo:1, slow:2, fast:3
[rA, rb, rc] = getReducedEigSystem(p.A, p.B, p.c, q, type_trunc);

m = p;
m.A = rA;
m.B = rb;
m.c = rc;

Fr_jac = @(T,t) evalRedSysandJac(T,m,sigma);
%forward function which evaluates our system and Jaobian

xr_start = zeros(size(m.A,1),1); 

tic;
trap_xr = trap(xr_start,dt_trap,t_stop,Fr_jac,tolerances_newt);
trapr_time = toc;
y1 = rc'*trap_xr;

norm(x_ref(:,end)-y1(:,end),1)

Temperature_im=zeros(Nx,Ny);
Temperature_im(brainmask)=y1(:,end);
Temperature_ref=zeros(Nx,Ny);
Temperature_ref(brainmask)=x_ref(:,end);
figure;
title('Reduced System: Slowest Modes')
h=imshow(Temperature_im);
max_temp = max(max(Temperature_im(:,:,1)));
set(h, 'AlphaData', Temperature_im(:,:,1)/max_temp);
drawnow;
figure;
title('Full System')
h=imshow(Temperature_ref);
max_temp = max(max(Temperature_ref(:,:,1)));
set(h, 'AlphaData', Temperature_ref(:,:,1)/max_temp);
drawnow;

%% 4) Demonstrate the TRAP solver for the Truncated Balanced Realization (TBR) reduced system
dt_trap=0.5;

%Calculate the reduced system
p.c=ones(size(p.B,1),1);
p.c=diag(p.c);
q=2;
tic;
syst = ss(p.A, p.B, p.c', 0);
systred = reduce(syst, q);

m = p;
m.A = systred.A;
m.B = systred.B;
m.c = systred.C;

Fr_jac = @(T,t) evalRedSysandJac(T,m,sigma);
%forward function which evaluates our system and Jaobian

xr_start = zeros(size(m.A,1),1); 

trap_xr = trap(xr_start,dt_trap,t_stop,Fr_jac,tolerances_newt);
y1 = systred.C*trap_xr;
toc;
norm(x_ref(:,end)-y1(:,end),1)

Temperature_im=zeros(Nx,Ny);
Temperature_im(brainmask)=y1(:,end);
Temperature_ref=zeros(Nx,Ny);
Temperature_ref(brainmask)=x_ref(:,end);
figure;
title('Reduced System: TBR')
h=imshow(Temperature_im);
max_temp = max(max(Temperature_im(:,:,1)));
set(h, 'AlphaData', Temperature_im(:,:,1)/max_temp);
drawnow;
figure;
title('Full System')
h=imshow(Temperature_ref);
max_temp = max(max(Temperature_ref(:,:,1)));
set(h, 'AlphaData', Temperature_ref(:,:,1)/max_temp);
drawnow;

%% 5) Demonstrate the TRAP solver for the Moment Matching reduced system
dt_trap=0.5;

%Calculate the reduced system
p.c=ones(size(p.B,1),1);
p.c=diag(p.c);
q=2;
tic;
V(:,1) = (p.A\p.B)/norm(p.A\p.B);
for i=2:q
    V(:,i) = p.A\V(:,i-1)/norm(p.A\V(:,i-1));
end

m = p;
m.A = V'*p.A*V;
m.B = V'*p.B;
m.c = V'*p.c;

Fr_jac = @(T,t) evalRedSysandJac(T,m,sigma);
%forward function which evaluates our system and Jaobian

xr_start = zeros(size(m.A,1),1); 

tic;
trap_xr = trap(xr_start,dt_trap,t_stop,Fr_jac,tolerances_newt);
trapr_time = toc;
y1 = m.c'*trap_xr;
toc;
norm(x_ref(:,end)-y1(:,end),1)

Temperature_im=zeros(Nx,Ny);
Temperature_im(brainmask)=y1(:,end);
Temperature_ref=zeros(Nx,Ny);
Temperature_ref(brainmask)=x_ref(:,end);
figure;
title('Reduced System: Moment Matching')
h=imshow(Temperature_im);
max_temp = max(max(Temperature_im(:,:,1)));
set(h, 'AlphaData', Temperature_im(:,:,1)/max_temp);
drawnow;
figure;
title('Full System')
h=imshow(Temperature_ref);
max_temp = max(max(Temperature_ref(:,:,1)));
set(h, 'AlphaData', Temperature_ref(:,:,1)/max_temp);
drawnow;

%% 6) Demonstrate the TRAP solver for a Nonlinear Reduced System
dt_trap=0.5;

%Calculate V for the reduced system
p.c=ones(size(p.B,1),1);
p.c=diag(p.c);
q=10;

tic;
mask_brain = p.mask(p.mask>0); 
%Portions of the brain over which we simulate, will be 
%N x 1.  Essentially we are picking out portion of mask
%which just corresponds to decision voxels and not air.

vessel_indices = find(mask_brain>2);
%Indices corresponding to vessel voxels in the temperature
%vector and the Jacobian Matrix

vessel_mask = zeros(size(p.A));
vessel_mask(vessel_indices) = 1; 
%mask all places where we have a vessel 

% We're going to linearize around 0, as our heat fluctuation around that
% point will tend to be trivial.
J = p.A;
for ii = 1:length(vessel_indices)
    vidx = vessel_indices(ii);
    J(vidx,vidx) = J(vidx,vidx) + p.dx2*sigma* 2;
end

V(:,1) = (J\p.B)/norm(J\p.B);
for i=2:q
    V(:,i) = J\V(:,i-1)/norm(J\V(:,i-1));
end

m = p;
m.A = V'*J*V;
m.B = V'*p.B;
m.c = V'*p.c;

%We're going to use the linear model everywhere but at the blood vessels
Fr_jac = @(T,t) evalRedSysandJac(T,m,sigma);
%forward function which evaluates our system and Jaobian

xr_start = zeros(size(m.A,1),1); 

tic;
trap_xr = trap(xr_start,dt_trap,t_stop,Fr_jac,tolerances_newt);
trapr_time = toc;
y1 = trap_xr;
y1 = m.c'*trap_xr;
toc;
norm(x_ref(:,end)-y1(:,end),1)

Temperature_im=zeros(Nx,Ny);
Temperature_im(brainmask)=y1(:,end);
Temperature_ref=zeros(Nx,Ny);
Temperature_ref(brainmask)=x_ref(:,end);
figure;
title('Reduced System: Nonlinear Approximation')
h=imshow(Temperature_im);
max_temp = max(max(Temperature_ref(:,:,1)));
set(h, 'AlphaData', Temperature_im(:,:,1)/max_temp);
drawnow;
figure;
title('Full System')
h=imshow(Temperature_ref);
max_temp = max(max(Temperature_ref(:,:,1)));
set(h, 'AlphaData', Temperature_ref(:,:,1)/max_temp);
drawnow;