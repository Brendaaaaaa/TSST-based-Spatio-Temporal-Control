function [alphad]=DSCFunction(dt,alphad,alpha,tau)
dalphad=(alpha-alphad)/tau;
alphad=alphad+dalphad*dt;