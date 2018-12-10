function [p,x_start,t_start,t_stop,max_dt_FE] = getParam_2DHeating(dx,dy,Nx,Ny,tmsize)

% Luca Daniel, MIT 2018 heatbar
% modified by Zijing 2018/9/18

%tmsize ~ tumorsize, dimension of the square tumor, assume to be even.

gamma =0.1; %thermal capacitance per unit length of the bar
km    =0.1; %thermal conductance through metal per unit length of the bar
Rair  = 100; %thermal resistance to air 
RcTissue     = (1/km)*dx*dy; %resistance of tissue
RcTumor       = (1/km)*dx*dy*.1;
Cstore = gamma*dx*dy; %the longer the section the larger the storage
     
N = Nx*Ny;
p.A    = zeros(N,N);
% coupling resistors Rc between (x,y) and (x+1,y), (x,y+1)
for x = 1:Nx
    for y=1:Ny
        
        if(abs(x-ceil(Nx/2)) <= tmsize/2 && abs(y-ceil(Ny/2)) <= tmsize/2)
            Rc = RcTumor;
        else
            Rc = RcTissue;
        end
        
        idx=(y-1)*Nx+x;
        xx=idx+1;
        yy=idx+Nx;
        if x<Ny
            p.A(idx,idx) = p.A(idx,idx)+(+1/Rc);
            p.A(idx,xx) = p.A(idx,xx)+(-1/Rc);
            p.A(xx,idx) = p.A(xx,idx)+(-1/Rc);
            p.A(xx,xx) = p.A(xx,xx)+(+1/Rc);
        end
        if y<Ny
            p.A(idx,idx) = p.A(idx,idx)+(+1/Rc);
            p.A(idx,yy) = p.A(idx,yy)+(-1/Rc);
            p.A(yy,idx) = p.A(yy,idx)+(-1/Rc);
            p.A(yy,yy) = p.A(yy,yy)+(+1/Rc);
        end
    end
end
% leakage resistor Rair between edge points and air
for x = 1:Nx
    y=2;
    idx=(y-1)*Nx+x;
    p.A(idx,idx) = p.A(idx,idx) + 1/Rair;
    y=Ny-1;
    idx=(y-1)*Nx+x;
    p.A(idx,idx) = p.A(idx,idx) + 1/Rair;
end
for y = 1:Ny
    x=2;
    idx=(y-1)*Nx+x;
    p.A(idx,idx) = p.A(idx,idx) + 1/Rair;
    x=Nx-1;
    idx=(y-1)*Nx+x;
    p.A(idx,idx) = p.A(idx,idx) + 1/Rair;
end

p.sqd   = eye(N,1);
p.B     = zeros(N,1);
idx=(floor(Ny/2))*Nx+ceil(Nx/2);
p.B(idx,1)= 1;


p.A     = p.A/Cstore; % note this will give a 1/dz^2 in A
							  % also pay attention to the negative sign
p.B     = p.B/Cstore;  % note this is important to make sure results
                       % will not depend on the number of sections N
      

x_start = ones(Nx,Ny)*15;
x_start([1,Nx],:)=0;
x_start(:,[1,Ny])=0;
x_start = x_start(:);
x_start = x_start.*0;
t_start = 0;

% slowest_timeconstant = min(abs(eig(p.A)));
% fastest_timeconstant = max(abs(eig(p.A)));

% to see steady state need to wait until the slowest mode settles
t_stop = 2000;

% usually Forward Euler is unstable for timestep>2/fastest_timeconstant
max_dt_FE = 0.1;
