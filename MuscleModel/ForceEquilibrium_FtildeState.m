% This function derives the Hill-equilibrium.
% More details in De Groote et al. (2016): DOI: 10.1007/s10439-016-1591-9
%
% Author: Antoine Falisse, modified by Gil Serrancolí to include tendon
% damping
% Date: 30/07/2022
% 
function [err, FT, Fce, Fiso, vMmax, Fpetilde, lMtilde, lTtilde] = ...
    ForceEquilibrium_FtildeState(a,fse_el,dfse_el,lMT,vMT,FMo_m,lMo_m,lTs_m,...
    alphao_m,vMmax_m,Fvparam,Fpparam,Faparam,fiber_damping,tendon_damping)

% FMo = ones(size(a,1),1)*params(1,:);
FMo = ones(size(a,1),1)*FMo_m;
lMo = ones(size(a,1),1)*lMo_m;
lTs = ones(size(a,1),1)*lTs_m;
alphao = ones(size(a,1),1)*alphao_m;
vMmax = ones(size(a,1),1)*vMmax_m;
Atendonsc = 35;
Atendon = ones(size(a,1),1)*Atendonsc;
volM = FMo.*lMo;
% massM = volM.*(1059.7)./(tension*1e6);

% Inverse tendon force-length characteristic
lTtilde = log(5*(fse_el + 0.25))./Atendon + 0.995;

% Hill-type muscle model: geometric relationships
lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
lMtilde = lM./lMo;

% Active muscle force-length characteristic
b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);
b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtilde-b23;
den3 = b33+b43*lMtilde;
FMtilde3 = b13*exp(-0.5*num3.^2./den3.^2);
num1 = lMtilde-b21;
den1 = b31+b41*lMtilde;
FMtilde1 = b11*exp(-0.5*num1.^2./den1.^2);
num2 = lMtilde-b22;
den2 = b32+b42*lMtilde;
FMtilde2 = b12*exp(-0.5*num2.^2./den2.^2);
FMltilde = FMtilde1+FMtilde2+FMtilde3;
Fiso = FMltilde;
% Active muscle force-velocity characteristic
vT = lTs.*dfse_el./(7*exp(35*(lTtilde-0.995)));
cos_alpha = (lMT-lTs.*lTtilde)./lM;
vM = (vMT-vT).*cos_alpha;
vMtilde = vM./vMmax;
e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);
FMvtilde = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;
% Active muscle force
% d = 0.01; % damping coefficient
Fcetilde = a.*FMltilde.*FMvtilde +fiber_damping.*vMtilde;
Fce = FMo.*Fcetilde;

% Passive muscle force-length characteristic
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtilde - 0.10e1) / e0);
% Passive muscle force
Fpetilde = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);

% Passive damping tendon force
fse_dam=(tendon_damping).*vT;

% Tendon force
fse=fse_el+fse_dam;

% Muscle force (non-normalized)
% FM = FMo.*(Fcetilde+Fpetilde);

% Tendon force
FT = fse .* FMo;

% Equilibrium between muscle and tendon forces
% err =  FM.*cos_alpha-FT;
err =  (Fcetilde+Fpetilde).*cos_alpha-fse;

end
