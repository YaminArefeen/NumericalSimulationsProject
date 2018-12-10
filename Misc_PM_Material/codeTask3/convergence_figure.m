n=10;

mydir  = pwd;
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end)-1); %navigate to folder above
addpath([newdir '/codeTask1']); %add appropriate directory 
addpath Utils;


cheb = zeros(n,1);
eigcc = zeros(n,1);
gcr = zeros(n,1);
for i = 20:20+n
    [A,b] = generate_matrix(3,0,i,15);
    [te,cc] = estimate_bounds(A,1e-6);
    eigcc(i-19,1) = te;
    cheb(i-19,1) = cc;
    [x,r_norms,gcr(i-19,1)] = tgcr(A,b,tol,maxiters);
end

hold on;
plot(cheb)
plot(eigcc)
plot(gcr)
legend('Condest', 'True Eigenvalues', 'GCR')