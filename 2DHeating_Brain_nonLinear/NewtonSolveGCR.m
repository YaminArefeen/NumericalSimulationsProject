function [x0] = NewtonSolveGCR(F,x0,p,sigma,tol,maxiters,e,d)
    for iter=1:maxiters
        iter
        f = F(x0);
        [dx,r_norms] = tgcr(F,x0,e,-f,d,maxiters);
        nf(iter)=norm(f,inf);
        ndx(iter)=norm(dx,inf);
        new_x(:,iter)=x0+dx;
        x0=new_x(:,iter);
        if nf(iter) < tol && ndx(iter) < tol
             % check for convergence
            fprintf('Converged in %d iterations\n',iter);
            break; 
        end
    end

if iter==maxiters % check for non-convergence
    fprintf('Non-Convergence after %d iterations!!!\n',iter); 
end

end
