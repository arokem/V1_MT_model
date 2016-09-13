function r=vonMises(T,Tpref,sigma)

%function r=vonMises(T,Tpref,sigma,b)
%Evaluates the Von Mises distribution, given: 
%<T> From [0,2pi] - the point in x to be evaluated
%<Tpref> The mean of the distribution
%<sigma> The bandwidth of the distribution
%
%100907 ASR made it


kappa=1/sigma;

r=exp(kappa*cos(T-Tpref))/(2*pi*besseli(0,kappa));