%moving interface hw3, phi evolution solver
% dphidt + F(ext) |gradphi| = 0
% simplest first order scheme(reference, evolution, implementation and
% application of level set and fast marching methods for advancing fronts)
% phi^n+1 = phi^n - dt(max(Fij,0)grad^+(phi) + min(Fij,0)grad^-(phi))
% grad^+(phi) = (max(-phi^+x, phi^-x,0)^2 + max(phi^-y, -phi^+y,0)^2)^1/2;
% grad^-(phi) = (max(phi^+x, -phi^-x,0)^2 + max(phi^+y, -phi^-y,0)^2)^1/2;
% phi^+x, forward, phi^-x, backward

function phi = evolution(phiold,F,dt)

global indexi indexj im ip
global  h

grad1 = zeros(indexi,indexj);
grad2 = zeros(indexi,indexj);
rhs = zeros(indexi,indexj);
phi = zeros(indexi, indexj);  

dmx = zeros(indexi,indexj);
dpx = zeros(indexi,indexj);
dmy = zeros(indexi,indexj);
dpy = zeros(indexi,indexj);


for i= 1: indexi
     for j=1:indexj         
               dmx(i,j) = (phiold(i,j) - phiold(im(i),j))/h ; % x backward
               dpx(i,j) = (phiold(ip(i),j) - phiold(i,j))/h; %x forward
               dmy(i,j) = (phiold(i,j) - phiold(i,im(j)))/h; % y backward
               dpy(i,j) = (phiold(i,ip(j)) - phiold(i,j))/h; % y forward
     end
end

for i=1:indexi
    for j=1:indexj
%      
%         % grad^+(phi)               
%         grad1(i,j) = (max(max(-dpx(i,j),0),max(dmx(i,j),0))^2 + ...
%                       max(max(dmy(i,j),0), max(-dpy(i,j), 0))^2)^0.5;
%         %grad^-(phi)
%         grad2(i,j)=(max(max(dpx(i,j),0),max(-dmx(i,j),0)).^2 + ...
%                     max(max(dpy(i,j),0), max(-dmy(i,j),0)).^2)^0.5;

        % grad^+(phi)               
        grad1(i,j) = (max(-dpx(i,j),0)^2+max(dmx(i,j),0)^2 + ...
                      max(dmy(i,j),0)^2+ max(-dpy(i,j), 0)^2)^0.5;
        %grad^-(phi)
        grad2(i,j)=(max(dpx(i,j),0)^2+max(-dmx(i,j),0)^2 + ...
                    max(dpy(i,j),0)^2+ max(-dmy(i,j),0)^2)^0.5;
   
        % rhs of evolution
        rhs(i,j) = max(F(i,j),0)*grad1(i,j) + min(F(i,j),0)*grad2(i,j);
        % evolution update
        phi(i,j) = phiold(i,j) - dt * rhs(i,j);
    end
end








            





