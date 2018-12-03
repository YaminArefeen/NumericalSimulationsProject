function f = eval_f_nonlinearSystem(T,p,sigma_non)
%This function evaluates heat flows, based on the input 
%temperatures.  In other words, we are simply evaluating
%our forward function, f(T).
%~~~~~~~
%INPUTS
%~~~~~~~
%T ~ (N x 1), temperatures at the current state
%p.A ~ N x N, Linear portion of our system 
%p.B ~ N x 1, Vecotrized input, where we apply constant input
%p.dx2 ~ 1 x 1, Resolution squared of our image

%~~~~~~~
%Outputs
%~~~~~~~
%f ~ N x 1, Heat flows based on the current state.

mask_brain = p.mask(p.mask>0); 
%Portions of the brain over which we simulate.

vessel_indices = find(mask_brain>2);
%Indices corresponding to vessel voxels in the temperature
%vector and the Jacobian Matrix

vessel_mask = zeros(size(T(:)));
vessel_mask(vessel_indices) = 1; 
%mask all places where we have a vesseel 

f = p.A*T(:)+ vessel_mask.*...
    (p.dx2*(exp(sigma_non.*T(:))-exp(-sigma_non.*T(:))))+p.B;
%Evaluating the linear and nonlienar portions of our system
%seperately. Only evaluate nonlinearity at vessel locations,
%based on the convention of nonlineary we've defined.
