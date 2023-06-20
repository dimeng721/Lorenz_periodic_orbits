% Input is a discrete curve.  Uses approximations to curvature, torsion and
% arclength derivatives from Boutin and returns these values.

function [ kappa_, kappa_s_, tauVal_, tau_sVal ] = compsig( pointVec )
% This function returns kappa/similar kappa, kappas/similar kappas,
% tau/similar tau, and taus/similar taus
%   Detailed explanation goes here

    function d=dist(point1, point2)
        % point1 and point2 are two 3-by-1 colum vectors
        d=sqrt((point1(1)-point2(1))^2+(point1(2)-point2(2))^2+(point1(3)-point2(3))^2);
    end
    
    % Compute area of a triangle using Heron's formula 
    function area = heron(a, b, c)
        s=(a+b+c)/2;
        area=sqrt(s*(s-a)*(s-b)*(s-c));
    end

    % Compute kappa
    function k = kap(point1, point2, point3)
        a=dist(point1, point2);
        b=dist(point2, point3);
        c=dist(point1, point3);
        area = heron(a, b, c);
        k=4*area/(a*b*c);
    end
   
    % Compute kappa_s
    function ks = kaps(point1, point2, point3, point4)
        a=dist(point1, point2);
        b=dist(point2, point3);
        d=dist(point3, point4);
        
        ks = 3*(kap(point2, point3, point4) - kap(point1, point2, point3))/(a+b+d);
    end

    % Compute tau
    function tau = tauF(point1, point2, point3, point4)
       a = dist(point1, point2);
       b = dist(point2, point3);
       c = dist(point1, point3);
       e = dist(point2, point4);
       d = dist(point3, point4);
       f = dist(point1, point4);
       
       % heron formula
       area_base = heron(a, b, c);
       
       % calculate H
       V = 1/6 * det([point1 - point4, point2 - point4, point3 - point4]);
       H = 3*V/area_base;
       
       tau = 6*H/(d*e*f*kap(point1, point2, point3));
       
    end

    % Compute tau_s
    function tau_s = tau_sF(point1, point2, point3, point4, point5, point6)
       
        a = dist(point2, point3);
        b = dist(point3, point4);
        d = dist(point4, point5);
        g = dist(point1, point2);
        
        h = kap(point2, point3, point4)*a*b/2;
        
        tauDiff = tauF(point3, point4, point5, point6) - tauF(point1, point2, point3, point4);
        tauKappa = tauF(point2, point3, point4, point5) * kaps(point2, point3, point4, point5)/(6*kap(point2, point3, point4));
        tau_s = 4*(tauDiff + (2*a + 2*b - 2*d - 3*h + g)*tauKappa)/(2*a + 2*b + 2*d + h + g);
    end

    function fai1s=fai1_s(point1, point2, point3, point4)
         a=dist(point1, point2);
         b=dist(point2, point3);
         d=dist(point3, point4);
         fai1s=[3*(kap(point2, point3, point4) - kap(point1, point2, point3))/(a+b+d)]/(kap(point2, point3, point4))^2;
        
    end

    function fais=fai_s(point1, point2, point3, point4,point5)
        g=dist(point2,point3);
        h=dist(point3,point4);
        i=dist(point4,point5);
        fais=3*(fai1_s(point2, point3, point4,point5)-fai1_s(point1, point2, point3, point4))/(g+h+i);
        fais=fais/kap(point2, point3, point4);
    end

    n=size(pointVec,2);

    %Curvature
    kappa=zeros(1,n);
    for t=2:(n-1)
        kappa(t)=kap(pointVec(:,t-1), pointVec(:,t), pointVec(:,t+1));
    end

    %Derivative of curvature with respect to arc length s
    kappa_s = zeros(1,n);
    for t = 3:(n-2)
        kappa_s(t) = kaps(pointVec(:, t-1), pointVec(:,t), pointVec(:, t+1), pointVec(:, t+2));
    end

    %Torsion
    tauVal = zeros(1, n);
    for t = 3:(n-2)
        tauVal(t) = tauF(pointVec(:, t-1), pointVec(:, t), pointVec(:, t+1), pointVec(:, t+2));
    end

    %Derivative of torsion with respect to arc length s
    tau_sVal = zeros(1,n);
    for t = 4:(n-3)
        tau_sVal(t) = tau_sF(pointVec(:, t-2),pointVec(:, t-1), pointVec(:, t), pointVec(:, t+1), pointVec(:, t+2), pointVec(:, t+3));
    end
    
    %Similar curvature
    kappa_=zeros(1,n);
    for t=2:(n-1)
        kappa_(t)=kappa_s(t)/(kappa(t)^2);
    end
    
    %Similar curvature 2
%     fi_s= zeros(1,n);
%     for t=4:(n-3)
%        fi_s(t)=fai_s(pointVec(:, t-1), pointVec(:,t), pointVec(:, t+1), pointVec(:, t+2),pointVec(:, t+3));   
%     end 
    
    % Derivative of similar curvature with respect to similar arc length s
    kappa_s_=zeros(1,n);
    kappa_ss=zeros(1,n);
    for tt=3:(n-3)
        kappa_ss(1,tt)=(kappa_s(tt)-kappa_s(tt+1))/dist(pointVec(:,tt),pointVec(:,tt+1));
    end
    for t=3:(n-2)
        kappa_s_(1,t)=(kappa_ss(1,t)*kappa(1,t)-2*(kappa_s(1,t))^2)/(kappa(1,t)^4);
    end
    
    
    %Similar torsion
    tauVal_=zeros(1,n);
    for t=3:(n-2)
       tauVal_(t)=tauVal(t)/kappa(t) ;
    end

end

