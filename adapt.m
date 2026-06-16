function [hat_theta,hat_D]=adapt(dt,mu,hat_theta,hat_D,f2,z2,beta,b1)
global kappa rho
dhat_theta=b1*rho*f2^2*z2^2*beta^2/mu-hat_theta;
dhat_D=(1-b1)*rho*f2^2*z2^2*tanh(z2^2/kappa)-hat_D;

hat_theta=hat_theta+dt*dhat_theta;
hat_D=hat_D+dt*dhat_D;