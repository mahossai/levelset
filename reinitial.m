%moving interface hw3, phi(re)initialization
% dphi/dt + sign(phi0)(|gradphi|-1) = 0

function phi = reinitial(phiold,tcont)

global indexi indexj ip im
global h 
global x y 

outstep = 10;
dt = 0.5*h; 
pn = zeros(indexi,indexj);
phi = phiold;


if tcont == 1  
    tot = ceil(10/dt);
    for j1 = 1:indexi
    for j2=1:indexj
         pn(j1,j2) = sign((x(j1)^2 + y(j2)^2)^0.5-1);
    end
    end
    
elseif floor(tcont/outstep)>floor((tcont-1)/outstep)
       fprintf('do reinitial\n');
       tot = ceil(10/dt);
        for i= 1:indexi
            for j=1:indexj
                if phiold(i,j)*phiold(ip(i),j) <= 0 || phiold(i,j)*phiold(im(i),j) <= 0 ...
                   || phiold(i,j)*phiold(i,ip(j))<= 0 || phiold(i,j)*phiold(i,im(j))<= 0
                 
                else
                   phi(i,j) = 0;
                end
            end
        end
        pn = sign(phiold);
else
    return
end
      

for tao = 0:tot
    for j1 =1:indexi
        for j2 = 1:indexj
            
            if phi(j1,j2)*phi(im(j1),j2)<0 || phi(j1,j2)*phi(ip(j1),j2)<0 || ...
               phi(j1,j2)*phi(j1,im(j2))<0 || phi(j1,j2)*phi(j1,ip(j2))<0 ;
%               dd = phi(j1,j2)/sqrt(2);
%               phi(j1,j2) = phi(j1,j2) - dt*(pn(j1,j2)*abs(phi(j1,j2))- dd)/h;  
              
            else
               dmx = (phi(j1,j2) - phi(im(j1),j2))/h ; % x backward
               dpx = (phi(ip(j1),j2) - phi(j1,j2))/h; %x forward
               dmy = (phi(j1,j2) - phi(j1,im(j2)))/h; % y backward
               dpy = (phi(j1,ip(j2)) - phi(j1,j2))/h; % y forward
 
                  ap = max(dmx, 0);
                  bm = min(dpx, 0);
                  cp = max(dmy,0);
                  dm = min(dpy, 0);
                  am = min(dmx,0);
                  bp = max(dpx, 0);
                  cm = min(dmy, 0);
                  dp = max(dpy, 0);
               hh = max(pn(j1,j2),0)*((max(ap^2, bm^2)+max(cp^2, dm^2))^0.5 - 1)+ ...
                   min(pn(j1,j2),0)*((max(am^2,bp^2) + max(cm^2,dp^2))^0.5 - 1);
               phi(j1,j2) = phi(j1,j2) - dt*hh;
            end
        end
    end
end






