function all_x = trap(x_init,dt,T,f,tolerances_newt,x_ref)
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

des_acc = .1; %desired accuracy

acc_flag = 1;

if(nargin < 6) % no reference solution exists
    acc_flag = 0; %do not check accuracy
end

x = x_init; %current time integrated value as initial condition
t = 0; %Start time at t = 0;

c_tol =  tolerances_newt(1); %define tolerances for newt
dx_tol =  tolerances_newt(2); %define tolerances for newt
r_tol =  tolerances_newt(3); %define tolerances for newt

num_iters = round(T/dt);
all_x = zeros(length(x_init),num_iters);

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
    
%     Update time step and new value
    t = t + dt;
    x = solutions(:,end); %recall newtonNd returns all solutions
    all_x(:,ii) = x; %keep track of solutions over time
    
    fprintf('Time: %d',t)

    if(acc_flag) %if we have a reference, check accuracy
        accuracy = norm(x_ref-x)/norm(x_ref);
        fprintf('|| Accuracy: %f\n',accuracy)
        
        if(accuracy <= des_acc)
            fprintf('Desired Accuracy of %d reached\n'...
                ,des_acc)
            break
        end            
    end
    
    %break if we have reached desired end time
    if(t > T)
        fprintf('Have reached total simulation time\n')
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

