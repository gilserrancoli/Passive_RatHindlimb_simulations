function [err, FT, Fpe, FMltilde, FMvtilde, cos_alpha, lTtilde, fse_el, fse_dam] = ForceEquilibrium_lMtildeState_forDebug(a,lMtilde,vMtilde,lM_projected,lMT,vMT,params,lMo,lTs,kT,shift,fiber_damping,tendon_damping)
% Hill-type muscle model: equilibrium between muscle and tendon forces
% All muscle-tendon characteristics are fully described in the publication
% and its online supplement

if size(lMtilde,1)==38
    FMo = params(1,:)';
    vMtildemax = params(5,:)';
else
    FMo = params(1,:);
    vMtildemax = params(5,:);
end

% Hill-type muscle model: geometric relationships
lM = lMtilde.*lMo;
lT = lMT - lM_projected;
lTtilde = lT./lTs;
cos_alpha = (lMT-lT)./lM;

% Tendon force-length characteristic
fse_el = (exp(kT.*(lTtilde - 0.995)))/5-0.25+shift;

% get muscle force-length characteristic
[Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(lMtilde,vMtilde,vMtildemax,fiber_damping);

% Tendon damping force
vTtilde=vMT./lMo-vMtilde./cos_alpha;
fse_dam=tendon_damping*vTtilde;

% Tendon normalized force
fse = fse_el + fse_dam;

% Active muscle force
Fce = a.*FMltilde.*FMvtilde;

% Muscle force
FM = Fce+Fpe;
% Tendon force
FT = FMo.*fse;

% Equilibrium between muscle and tendon forces
% Fm*cos(alpha) = Ft
err =  FM.*cos_alpha-fse;

end