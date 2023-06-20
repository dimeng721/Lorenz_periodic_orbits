%% Observe the motion of a trajectory with its similarity signature curve.
%Slow running

% A_temp = A(:,100000:end);
% [ kappa_1, kappa_s_1, tauVal_1, tau_sVal1 ] = compsig( A_temp );
% C = B;
% figure(1)
% for i=1:size(C,2)
% subplot(1,2,1)
% plot3(C(1,i),C(2,i),C(3,i),'.','linewidth',2)
% hold on
% pause(0.00001)
% subplot(1,2,2)
% plot3(A_temp(1,i),A_temp(2,i),A_temp(3,i),'.','linewidth',2)
% hold on
% end

%% Generate quasi-periodic orbits

global temp     %Index of the starting point of the current window in the whole trajectory.
global Dist     %Distance
global piecewise
temp = [];
Dist=[];
piecewise=[];

j=0;
temp = [1];
% piecewise=[53282];
m = 635;    %Initial point

%In order to increase the speed of the calculation, the Lorenz flow can be 
%processed in segments. But it is easy to miss quasi-periodic orbits
while temp(end)+1450<length(A_temp(:,m:1800000))
    j=j+1;
    fenzu(A_temp(:,m:990000),j,C(:,m:990000))   
end
% end


function fenzu(A_temp,zoo,data)
global temp
global Dist
global kappa
global piecewise

dist=zeros(1,1450);     %Record distance
if zoo == 1
    temp(1) = 1;
end

%Calculate distance
for i=1:1450
    dist(i) = -sqrt((data(1,temp(zoo))-data(1,i+temp(zoo)))^2 + (data(2,temp(zoo))-data(2,i+temp(zoo)))^2 + (data(3,temp(zoo))-data(3,i+temp(zoo)))^2);
%     if dist(i)>100
end
[peak, loca] = findpeaks(dist);
[val, index] = max(peak);
loca = loca(index);
% subplot(2,2,1);
% subplot(1,3,1)
% plot(1:length(dist),-dist)
% hold on;
temp=[temp temp(end) + loca(1)];

dist1=norm(A_temp(:,temp(end))-A_temp(:,piecewise(end)));
if dist1<0.33
    piecewise=[piecewise temp(end)];
end

if length(temp)>=2
%     subplot(2,2,2);
    subplot(1,2,2);
    plot3(data(1,temp(end-1):temp(end)),data(2,temp(end-1):temp(end)),data(3,temp(end-1):temp(end)));
    grid on;
    xlabel('$\tilde{\kappa}$', 'interpreter', 'latex');
    ylabel('$\tilde{\kappa _s}$', 'interpreter', 'latex')
    zlabel('$\tilde{\tau}$', 'interpreter', 'latex')
    hold on;
%     subplot(2,2,3);
    subplot(1,2,1);
    plot3(A_temp(1,temp(end-1):temp(end)),A_temp(2,temp(end-1):temp(end)),A_temp(3,temp(end-1):temp(end)),'color', [0 0.4470 0.7410]);     %,'color' [0 0.4470 0.7410]
%     view(0,0);
%     xlabel('x')
%     zlabel('z')
    xlabel('x');
    ylabel('y')
    zlabel('z')
    grid on;
    hold on;
    plot3(A_temp(1,temp(end)),A_temp(2,temp(end)),A_temp(3,temp(end)),'x','color',[218/255 83/255 23/255]);
    hold on;
end
Dist=[Dist sum(dist)/1000];
% subplot(2,2,4)
% plot(1:length(Dist),Dist)
% hold on
if length(Dist)>1
    kappa=[kappa Dist(end-1)-Dist(end)];
end
% figure(2);
% plot()

end

%% Checking quasi-periodic orbits that satisfy the conditions
% piecewise={};
% for j = 1:length(temp)-1
%     piecewise(j)={[temp(j)]};
% end
% for j = 1:length(temp)-1
%     for i=j+1:length(temp)-1
%         if norm(A_temp(:,1997+piecewise{j}(end))-A_temp(:,1997+temp(i)))<0.1&&temp(i)-piecewise{j}(end)<17500
%            piecewise{j}=[piecewise{j} temp(i)];
%         end
%     end
% end

