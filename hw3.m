% moving interface hw3 main function

%  function hw3

clear all
close all
clc
global x y 
global h 
% initial configuration of phi
 phi = initial();

tf=60; 
dt = 0.5*h; 

for tcont = 1: tf
    % reinitialization phi
     phi = reinitial(phi, tcont);
    % velocity extension
%      F = extension(phi,signphi,Fold,tcont);
      F = extension(phi,dt);
    % phi evolution
     phi = evolution(phi,F, dt);
     
     if floor(tcont/10) > floor((tcont-1)/10)
         figure(1)
     contour(x,y,phi,[0,0]), hold on
     axis([-5 5 -5 5]); axis equal
     xlabel('x position');
     ylabel('y position');
     title('phi advancing each 10 steps');
     fprintf('tcont is %d\n', tcont);
     end
end
    
    
    
    
    
    
    
    