function estimate_bounds(A,tol)

% start with the Gershgorin Circle Estimation
true_small = 100;
true_large = 0;
for i=1:length(A)
    center = A(i,i);
    radius = sum(abs(A(:,i)))-abs(center);
    est_small = (center-radius);
    est_large = (center+radius);
    if est_small < true_small
        true_small = est_small;
    end
    if est_large > true_large
        true_large = est_large;
    end
end

if abs(true_small) > abs(true_large)
    fr = (sqrt(true_small/true_large)-1)/(sqrt(true_small/true_large)+1);
    gc_bound = log(tol/2)/log(fr);
else
    fr = (sqrt(true_large/true_small)-1)/(sqrt(true_large/true_small)+1);
    gc_bound = log(tol/2)/log(fr);
end

disp('The GC Bound is: ');
disp(gc_bound);

% with the true eigenstuff
V = eigs(sparse(A),1,'largestabs')
V2 = eigs(sparse(A),1,'smallestabs');

fr = (sqrt(V/V2)-1)/(sqrt(V/V2)+1);
gc_bound = log(tol/2)/log(fr);

disp('The true eigenstuff Bound is: ');
disp(gc_bound);

% with condest
V = condest(sparse(A));

fr = (sqrt(V)-1)/(sqrt(V)+1);
gc_bound = log(tol/2)/log(fr);

disp('The condest Bound is: ');
disp(gc_bound);

end