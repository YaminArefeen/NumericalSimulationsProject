function visualizeResults2D(t,X,n,Nx,Ny,plottype)

% copyright Luca Daniel, MIT 2018


   subplot(2,1,1)
   %figure(1)
   plot(t(n),X(:,n),plottype); 
   hold on
   xlabel('time')
   subplot(2,1,2)
   %figure(2)
   imagesc(reshape(X(:,n),[Nx,Ny])); colorbar;
   %pause(1); 
   drawnow;