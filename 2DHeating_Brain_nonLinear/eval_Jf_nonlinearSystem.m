function J = eval_Jf_nonlinearSystem(T,p,sigma)
%In this function, we evalaute the Jacobian of our nonlinear
%system.  Since our nonlinear system is relatively benign
%we evaluate the Jacobian by just taking the linear Nodal 
%matrix, and adding the appropriate nonlinearities to the 
%rows of the diagonal which correspond to vessel nodes.
%~~~~~~
%Inputs
%~~~~~~
%T ~ N x 1, Current state of temperatures
%p.mask_brain ~ Nx x Ny, Total brain mask.  Note, this includes
%both brain voxels and air voxels, so (Nx)*(Ny) > N, since 
%the decision variables only exist on the brain
%p.A ~ N x N, Linear portion of our system.  Again, note that
%it is N x N, since we are only considering decision variables
%on the brain. 
%p.dx2 ~ 1 x 1, Resolution squared of our image, used for 
%nonlinearity
%~~~~~~~
%Outputs
%~~~~~~~
%J ~ N x N, Jacobian of our system evaluated at current T.

mask_brain = p.mask(p.mask>0); 
%Portions of the brain over which we simulate, will be 
%N x 1.  Essentially we are picking out portion of mask
%which just corresponds to decision voxels and not air.

vessel_indices = find(mask_brain>2);
%Voxel indices which correspond to vessels, based on the 
%convention defined.

J = p.A;
for ii = 1:length(vessel_indices)
    vidx = vessel_indices(ii);
    J(vidx,vidx) = ...
        J(vidx,vidx) + p.dx2.*...
            (sigma*exp(sigma.*T(vidx)) ... 
                + sigma*exp(-sigma.*T(vidx)));
end
%Our Jacobian will simply be the original nodal matrix + 
%the nonlinear term at voxels corresponding to vessels.
%Since the nonlinearity only applies on a voxel by voxel
%basis and has no influence on connections between voxels,
%we only need to add the nonlinearity to the diagonal elements
%of the Jacobian which correspond to vessels.