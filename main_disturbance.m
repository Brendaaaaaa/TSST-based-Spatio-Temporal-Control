%%
% numerical simulation of scenario with disturbance
%%
clc;
clear all;
%%
i=1;
dt=0.001;

t(i)=0;

global l1 l2 rho k1 k2 k3 tau1 tau2 p T_f
l1=0.5;l2=0.5;
rho=8;
k1=0.1;k2=0.1;k3=10;
tau1=0.3;tau2=tau1;
p=2;
T_f=1.76;

t_end=1.7;

x1(i)=7.5;x2(i)=5;x3(i)=-10;
hat_theta(i)=1;hat_D(i)=1;
global mux1 h1 kappa w b
mux1=1;h1=1;kappa=0.1;w=1;b=1;
z11=1000;z12=2000;z21=1000;z22=3000;
%%
while t(i)<=t_end
    % TIME-VARYING FUNCTION
    mu(i)=1/(T_f-t(i))^p;
    % TARGET
    xd(i)=3*sin(0.5*t(i)); 
    dxd(i)=1.5*cos(0.5*t(i));
    ddxd(i)=-1.5*0.5*sin(0.5*t(i));
    disturbance(i)=cos(0.3*t(i))+2*x1(i)-3*x2(i);
    % FULL-STATE CONSTRAINTS AND INPUT SATURATION
    L1(i)=20*exp(-t(i))+10;dL1(i)=-20*exp(-t(i));
    H1(i)=18*exp(-t(i))+10;dH1(i)=-18*exp(-t(i));
    L2(i)=30;dL2(i)=0;
    H2(i)=20;dH2(i)=0;
    L3(i)=180;dL3(i)=0;
    H3(i)=150;dH3(i)=0;

    zeta1(i)=l1*(H1(i)*x1(i))/(H1(i)-x1(i))+l2*(L1(i)*x1(i))/(L1(i)+x1(i));
    zeta2(i)=l1*(H2(i)*x2(i))/(H2(i)-x2(i))+l2*(L2(i)*x2(i))/(L2(i)+x2(i));
    zeta3(i)=l1*(H3(i)*x3(i))/(H3(i)-x3(i))+l2*(L3(i)*x3(i))/(L3(i)+x3(i));
    
    alpha0(i)=l1*(H1(i)*xd(i))/(H1(i)-xd(i))+l2*(L1(i)*xd(i))/(L1(i)+xd(i));
    
    dalpha0(i)=(l1*(H1(i)^2)/(H1(i)-xd(i))^2+l2*(L1(i)^2)/(L1(i)+xd(i))^2)*dxd(i)-l1*(xd(i)^2*dH1(i))/(H1(i)-xd(i))^2+l2*(xd(i)^2*dL1(i))/(L1(i)+xd(i))^2;
    z1(i)=zeta1(i)-alpha0(i);
    f1(i)=l1*(H1(i)^2)/(H1(i)-x1(i))^2+l2*(L1(i)^2)/(L1(i)+x1(i))^2;
    g1(i)=-l1*(x1(i)^2*dH1(i))/(H1(i)-x1(i))^2+l2*(x1(i)^2*dL1(i))/(L1(i)+x1(i))^2;
    alpha1(i)=-f1(i)/(4*rho)*(x2(i)-zeta2(i))^2*z1(i)-f1(i)*z1(i)-mu(i)/(4*f1(i)*rho)*(g1(i)-dalpha0(i))^2*z1(i)-mu(i)*k1/(f1(i))*z1(i);
    
    if i==1
        alpha1d(i)=alpha1(i);
    end
    z2(i)=zeta2(i)-alpha1d(i);
    f2(i)=l1*(H2(i)^2)/(H2(i)-x2(i))^2+l2*(L2(i)^2)/(L2(i)+x2(i))^2;
    g2(i)=-l1*(x2(i)^2*dH2(i))/(H2(i)-x2(i))^2+l2*(x2(i)^2*dL2(i))/(L2(i)+x2(i))^2;
    beta(i)=(exp(-((z2(i)-mux1)^2)/(h1^2)))^2;
    alpha_s(i)=rho*f2(i)*z2(i)*hat_theta(i)*beta(i)^2;
    alpha_c(i)=rho*f2(i)*z2(i)*hat_D(i)*tanh(z2(i)^2/kappa);
    [b1(i)]=swt2(z1(i),z2(i),z11,z12,z21,z22);
    alpha2(i)=-f2(i)/(4*rho)*(x3(i)-zeta3(i))^2*z2(i)-mu(i)/(4*f2(i)*rho)*g2(i)^2*z2(i)-f2(i)*z2(i)-1/(f2(i)*tau1^2)*z2(i)-f1(i)/f2(i)*z1(i)-mu(i)*k2/(f2(i))*z2(i)-b1(i)/2*mu(i)*f2(i)*z2(i)-b1(i)*alpha_s(i)-(1-b1(i))*alpha_c(i);
    
    if i==1
        alpha2d(i)=alpha2(i);
    end
    z3(i)=zeta3(i)-alpha2d(i);
    f3(i)=l1*(H3(i)^2)/(H3(i)-x3(i))^2+l2*(L3(i)^2)/(L3(i)+x3(i))^2;
    g3(i)=-l1*(x3(i)^2*dH3(i))/(H3(i)-x3(i))^2+l2*(x3(i)^2*dL3(i))/(L3(i)+x3(i))^2;
    u(i)=-mu(i)/(4*f3(i)*rho)*g3(i)^2*z3(i)-1/(f3(i)*tau2^2)*z3(i)-f2(i)/f3(i)*z2(i)-mu(i)*k3/(f3(i))*z3(i);
    
    i=i+1;
    [x1(i),x2(i),x3(i)]=integrator_disturbance(dt,x1(i-1),x2(i-1),x3(i-1),u(i-1),disturbance(i-1));
    [alpha1d(i)]=DSCFunction(dt,alpha1d(i-1),alpha1(i-1),tau1);
    [alpha2d(i)]=DSCFunction(dt,alpha2d(i-1),alpha2(i-1),tau2);
    [hat_theta(i),hat_D(i)]=adapt(dt,mu(i-1),hat_theta(i-1),hat_D(i-1),f2(i-1),z2(i-1),beta(i-1),b1(i-1));
    t(i)=t(i-1)+dt;
end
%%
function [m2]=swt2(s1,s2,r11,r12,r21,r22)
global w b;
if abs(s1)<r11
    B1=1;
elseif (r11<=abs(s1))&&(abs(s1)<=r12)
    B1=(r12^2-s1^2)/(r12^2-r11^2)*exp(-((s1^2-r11^2)/(w*(r12^2-r11^2)))^(2*b));
else 
    B1=0;
end
if abs(s2)<r21
    B2=1;
elseif (r21<=abs(s2))&&(abs(s2)<=r22)
    B2=1*(r22^2-s2^2)/(r22^2-r21^2)*exp(-((s2^2-r21^2)/(w*(r22^2-r21^2)))^(2*b));
else 
    B2=0;
end
m2=B1*B2;
end