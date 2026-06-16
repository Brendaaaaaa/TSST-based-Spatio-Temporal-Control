function [x1,x2,x3]=integrator_disturbance(dt,x1,x2,x3,u,d)
dx1=x2;
dx2=x3+d;
dx3=u;

x1=x1+dx1*dt;
x2=x2+dx2*dt;
x3=x3+dx3*dt;