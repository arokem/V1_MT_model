function logL=logLikeTheta(theta)


global fit_makeTpref;
global fit_sigmaConn;
global fit_activityMT;

logL=[];
kappa=1./fit_sigmaConn;

for i=1:length(theta)
    logL=[logL -1*sum(kappa.*fit_activityMT.*cos(fit_makeTpref-theta(i)))/(sum(kappa))];
end

