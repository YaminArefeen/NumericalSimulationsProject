function [all_x,t_all] = dyn_trap_withVisualization(x_init,dt,T,f,tolerances_newt,...
        max_dt,T1_image,brainmask,max_temp,visualize)
%In this function, we utilize the trap rule to perform time integration
%x_init ~ N x 1, initial state for our system
%dt ~ 1 x 1, time step over which we wish to integrate
%T ~ 1 x 1, time for which we wish to integrate
%f ~ function handler which evaluates function and Jacobian of system
    %at a particular state
%tolerances_newt ~ 3 x 1, Tolerances for the newton algorithm
    %first component is tolerance for the function, second is 
    %tolerance for the change in x, and third is tolerance for 
    %relative change in x.
%x_ref ~ N x 1, Reference solution if it exists, for the sake
%of checking
if visualize
    T1_image=T1_image-min(T1_image(:));
    T1_image=T1_image./max(T1_image(:));
    figure(2)
    imshow(T1_image,[0 1],'InitialMag', 'fit'); colormap('gray');
    hold on
    Temperature_im=zeros([size(T1_image),3]);
end
des_acc = .1; %desired accuracy

x = x_init; %current time integrated value as initial condition
t = 0; %Start time at t = 0;

c_tol =  tolerances_newt(1); %define tolerances for newt
dx_tol =  tolerances_newt(2); %define tolerances for newt
r_tol =  tolerances_newt(3); %define tolerances for newt

num_iters = round(T/dt);
all_x = zeros(length(x_init),num_iters);
t_all = zeros(num_iters,1);

for ii = 1:num_iters %loop while t is less than total desired time
    %Part which we know in trap rule
    
    gamma = x + dt/2*f(x,t); 
    
    %Define a function which outputs correct Jacobian
    %and function eval for trap newton solve so that 
    %I can use my pre-existing NewtonNd code
    trap_fjpoisson = @(input,t) trap_Jacobian(f,input,gamma,dt,t);
    
    %Get next point in time integration by solving newton system
    %F(x) = 0
    solutions = newtonNd(trap_fjpoisson,x,c_tol,dx_tol,r_tol,t);
    num_newt_step = size(solutions,2);
    
    if(num_newt_step > 5)
        dt = dt/2;
    else
        dt = 2*dt;
    end

    dt = min(dt,max_dt); %set max dt to 50

%     Update time step and new value
    t = t + dt;
    x = solutions(:,end); %recall newtonNd returns all solutions
    all_x(:,ii) = x; %keep track of solutions over time    
    t_all(ii)=t;
   if visualize
       Temperature_im(brainmask)=x;
       imshow(T1_image,[0 1],'InitialMag', 'fit'); colormap('gray');
       hold on
       h=imshow(Temperature_im);
       set(h, 'AlphaData', Temperature_im(:,:,1)/max_temp); title(['Time at ',num2str(t)]);
       drawnow;       
       hold on
   end                
    %break if we have reached desired end time
    fprintf('------Time: %d -------\n',t)
    if(t > T)
        fprintf('Have reached total simulation time\n')
        break
    end
end
all_x = all_x(:,1:ii);
end

function [f_trap,J_trap] = trap_Jacobian(f,input,gamma,dt,t)
    %this function evaluates the funtion value and Jacobian
    %used for the newton interation inside of the trapezoidal
    %integration scheme.
    [f_eval,J_eval] = f(input,t); %Evaluate the function and Jacobian
    %at the current point.  Keep in mind that this is not the 
    %right Jacobian/function eval for the Newton Step

    f_trap = input - dt/2*f_eval - gamma;
    J_trap = eye(size(J_eval)) - dt/2*J_eval;   
end

