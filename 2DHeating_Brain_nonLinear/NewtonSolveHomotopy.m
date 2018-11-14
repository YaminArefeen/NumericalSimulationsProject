function [x0] = NewtonSolveHomotopy(F,x0,p,sigma,tol,maxiters)
    q0 = 0;
    G = @(x,q) q*F(x) + (1-q)*x;
    for iter=1:maxiters
        f = G(x0,q0);
        J = q0*eval_Jf_nonlinearSystem(x0,p,sigma) + (1-q0)*eye(size(p.A));
        dx = -f.'/J;
        nf(iter)=norm(f,inf);
        ndx(iter)=norm(dx,inf);
        new_x(:,iter)=x0+dx.';
        x0=new_x(:,iter);
        
        if q0 < 1
            q0 = round(q0,1) + 0.1;
        else
            q0 = 1;
        end
        nf(iter)
        ndx(iter)
        if nf(iter) < tol && ndx(iter) < tol && q0 == 1
             % check for convergence
            fprintf('Converged in %d iterations\n',iter);
            break; 
        end
    end

if iter==maxiters % check for non-convergence
    fprintf('Non-Convergence after %d iterations!!!\n',iter); 
end

end