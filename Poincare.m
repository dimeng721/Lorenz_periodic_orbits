%Poincare map
%% Lorenz equations
clear
clc
close all
%Generate trajectories.
global m
global h

m=0.63;      %Time length
h=m/(0.001);      %Iteration number
x0=0:0.001:m;

y0=[17.0812943860444;27.6480832190505;27];
tspan=[0,m];
[t1,y1]=rk4(@lorenz,tspan,y0,h);
[t2,y2]=ode45(@lorenz,tspan,y0);
% plot3(y1(:,1),y1(:,2),y1(:,3));
% figure();plot3(y2(:,1),y2(:,2),y2(:,3));
y1=y1';
Lx=y1(1,:);
Ly=y1(2,:);
Lz=y1(3,:);

Plane=[0;0;1;-27];%Plane z=27 (positive direction).

[tP_List,yP_List]=Solve_Poincare(x0,y1,Plane);

%% Plot
%1 Poincare section 
%The first few points were not yet stable and did not reflect the 
%characteristics of the system, so they were abandoned
figure()
plot(yP_List(1,:),yP_List(2,:),'.')
% xlim([-1,0.6])
% ylim([-0.8,0.2])


function [t,y]=rk4(f,tspan,y0,n)
%Fourth-order Runge-Kutta
%f:dy/dt=f(t,y);
%Solution intervalt:  span:[t0,tf];
%Initial value: y0=y(t0)=[x0,y0,z0];
%Number of sampling points: n.
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
p=10;
r=28;
b=8/3;
dx_dt=p*(y(2)-y(1));
dy_dt=r*y(1)-y(2)-y(1)*y(3);
dz_dt=y(1)*y(2)-b*y(3);
dy_dt=[dx_dt;dy_dt;dz_dt];
end

%Distance from point to plane
function Dis=DistancePlane(xk,Plane)
% xk£¬Coordinate points£¬In the case of 3-dimensional coordinates, the size is a 3*N matrix.
% Plane£¬A plane of the form Ax+By+Cz+D=0

N=size(xk,2);
xk2=[xk;ones(1,N)];
Dis=dot(xk2,Plane*ones(1,N),1)./norm(Plane(1:end-1));
end

function [tP_List,yP_List]=Solve_Poincare(t,y,Plane)
%Cross-sectional equations: z=27 (Changed to take account of actual circumstances)
% Plane=[0;0;1;27]; Generally a plane perpendicular to an axis
%Generally only recorded from negative to positive crossing¡£
%If you want to record in reverse as well, you can set Plane=-Plane.

%Interpolation to obtain the intersection of a line and a surface
yP_List=[];
tP_List=[];
Dis=DistancePlane(y,Plane);
N=size(y,2);
for k=1:N-1
    if Dis(k)<=0 && Dis(k+1)>0
        t0=t(k);t1=t(k+1);
        yP0=y(:,k);yP1=y(:,k+1);
        Dis0=Dis(k);Dis1=Dis(k+1);
        yP=yP0+(yP1-yP0)/(Dis1-Dis0)*(0-Dis0);
        tP=t0+(t1-t0)/(Dis1-Dis0)*(0-Dis0);
        yP_List=[yP_List,yP];
        tP_List=[tP_List,tP];
    end
end
end
