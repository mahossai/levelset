%moving interface hw3 part2: velocity extension

%  function F = extension(phi,signphi,Fold,tcont)

function F = extension(phi,dt)

% function F = extension(phi)
global indexi indexj im ip
global h
global x y 

  
F = zeros(indexi,indexj);
gradphix = zeros(indexi,indexj);
gradphiy = zeros(indexi,indexj);
gradphivalue = zeros(indexi,indexj);
nx = zeros(indexi,indexj);
ny = zeros(indexi,indexj);


signphi = sign(phi);

for j1 = 1:indexi
    for j2 = 1: indexj
 % centeral difference gradphi        
        gradphix(j1,j2) = (phi(ip(j1),j2) - phi(im(j1),j2))/2/h;
        gradphiy(j1,j2) = (phi(j1,ip(j2)) - phi(j1,im(j2)))/2/h;
        gradphivalue(j1,j2) =(gradphix(j1,j2)^2 + gradphiy(j1,j2)^2)^0.5;
        nx(j1,j2) = gradphix(j1,j2) / gradphivalue(j1,j2);
        ny(j1,j2) = gradphiy(j1,j2) / gradphivalue(j1,j2);
    end
end

%cos(theta)=nx;
%sin(theta)=-ny;

for i = 1:indexi
    for j=1:indexj
        if phi(i,j)*phi(ip(i),j) <= 0 || phi(i,j)*phi(im(i),j) <= 0 || ...
                   phi(i,j)*phi(i,ip(j))<= 0 || phi(i,j)*phi(i,im(j))<= 0
                            
%           F(i,j) = 2 - 4*nx(i,j)^2*ny(i,j)^2;
           F(i,j) = 2 + (1+ nx(i,j))/2*(2*nx(i,j)-1)^2; 
         %  F(i,j) = 10+(1+ nx(i,j))/2*(2*nx(i,j)-1)^2; 
%          F(i,j) = 1+10*(1-4*nx(i,j)^2*ny(i,j)^2);
         
        end
    end
end

 
 Fxp = zeros(indexi, indexj);
 Fxb = zeros(indexi,indexj);
 Fyp = zeros(indexi, indexj);
 Fyb = zeros(indexi, indexj);
 
count = 300; 

for tt = 1: count
for j1 = 1: indexi
    for j2 = 1: indexj  
        Fxp(j1,j2) = (F(ip(j1),j2) - F(j1,j2))/h; % forward of x
        Fxb(j1,j2) = (F(j1,j2) - F(im(j1),j2))/h; % backward of x
        Fyp(j1,j2) = (F(j1,ip(j2)) - F(j1,j2))/h; % forward of y
        Fyb(j1,j2) = (F(j1,j2) - F(j1,im(j2)))/h; % backward of y      
    end
end

for j1 = 1: indexi
    for j2 = 1: indexj
    if phi(j1,j2)*phi(im(j1),j2)<=0 || phi(j1,j2)*phi(ip(j1),j2)<=0 ||...
       phi(j1,j2)*phi(j1,im(j2))<= 0 || phi(j1,j2)*phi(j1,ip(j2))<= 0 
   
    else
    F(j1,j2) = F(j1,j2) - dt * (max(signphi(j1,j2)*nx(j1,j2),0)*Fxb(j1,j2) + ...
                  min(signphi(j1,j2)*nx(j1,j2),0)*Fxp(j1,j2))- ...
                  dt*(max(signphi(j1,j2)*ny(j1,j2),0)*Fyb(j1,j2) + ...
                  min(signphi(j1,j2)*ny(j1,j2),0)*Fyp(j1,j2));
    end
    end
end
end




