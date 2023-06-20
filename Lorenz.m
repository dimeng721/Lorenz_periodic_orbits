global m
global h

m=4e3;      %Iteration length
h=m/(0.001);      %Number of iterations

y0=[0.1;0.1;0.1];   %Starting point of the system
tspan=[0,m];
[t1,y1]=rk4(@lorenz,tspan,y0,h);    %Fourth Order Runge-Kutta Method
[t2,y2]=ode45(@lorenz,tspan,y0);    %ODE
plot3(y1(:,1),y1(:,2),y1(:,3));     %Lorenz Attractor
% figure();plot3(y2(:,1),y2(:,2),y2(:,3));

function [t,y]=rk4(f,tspan,y0,n)
%Fourth Order Runge-Kutta Method
%f:dy/dt=f(t,y); 
%Solution intervaltspan:[t0,tf]; 
%Initial value: y0=y(t0)=[x0,y0,z0];
%Number of sampling points n;
h=(tspan(2)-tspan(1))/n;
t=tspan(1):h:tspan(2);
y=zeros(length(y0),length(t));
y(:,1)=y0;
for m=2:length(t)
    k1=f(t(m-1),y(:,m-1));
    k2=f(t(m-1)+h/2,y(:,m-1)+k1*h/2);
    k3=f(t(m-1)+h/2,y(:,m-1)+k2*h/2);
    k4=f(t(m-1)+h,y(:,m-1)+k3*h);
    y(:,m)=y(:,m-1)+(k1+2*k2+2*k3+k4)*h/6;
end
t=t';
y=y';
end

function dy_dt=lorenz(t,y)
%y=[x,y,z] is the matrix storing the three variables.
p=10;
r=28;
b=8/3;
dx_dt=p*(y(2)-y(1));
dy_dt=r*y(1)-y(2)-y(1)*y(3);
dz_dt=y(1)*y(2)-b*y(3);
dy_dt=[dx_dt;dy_dt;dz_dt];
end