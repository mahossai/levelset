

function phi = initial()

global indexi indexj ip im
global x y xx yy
global h N_interface
% global phi

%----------------inital configuration-------------------

% global domain descritization
N = 100 ;
h = 10/(N-1); % grid spacing
x = -5: h: 5;
y = x;
% ip = zeros(N,1);
% im = zeros(N,1);
for i=1:N
    ip(i) = i+1;
    im(i) = i-1;
end
im(1)=N;
ip(N)=1;

indexi = N;
indexj = N;
phi= zeros(indexi,indexj); % phi initial

% interface descritization
N_interface = 400; 
s = linspace(0,2*pi,N_interface+1);
s(N_interface+1) = [];
xx = cos(s);
yy = sin(s);

%-------------setting initial phi(only near interface) as SDF
        for i=1:N_interface         
            labelx = ceil((xx(i)+5)/h);
            labely = ceil((yy(i)+5)/h);          
            phi(labelx:labelx+1, labely:labely+1) = 1;
        end
        distance = 2*ones(N_interface,1);    
% ------------- SDF -----------------
    for j1= 1:indexi
        for j2 = 1:indexj
            if phi(j1,j2) == 1
                for i = 1: N_interface;   
                      distance(i) = sum(([x(j1),y(j2)] - [xx(i),yy(i)]).^2).^0.5;
                end    
                Dij = min(distance);
                if x(j1)^2+y(j2)^2-1<0   % in circle
                   Dij = - Dij;
                end                      
                      phi(j1,j2) = Dij;
            end
        end
    end                  
%     surf(x,y,phi);