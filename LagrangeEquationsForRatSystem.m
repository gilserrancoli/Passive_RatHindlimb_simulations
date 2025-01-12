syms thetah(t) thetat(t) thetaf(t) g mfem mtib mfoot Ifem Itib Ifoot...
    lfem ltib rfem rtib rfoot lambda1 lambda2 h ...
    thetah_aux thetat_aux thetaf_aux ...
    thetadoth_aux thetadott_aux thetadotf_aux ...
    thetad2doth_aux thetad2dott_aux thetad2dotf_aux Mh Mk Ma;

thetadoth=diff(thetah);
thetadott=diff(thetat);
thetadotf=diff(thetaf);
thetad2doth=diff(thetah,2);
thetad2dott=diff(thetat,2);
thetad2dotf=diff(thetaf,2);

vGfem=[rfem*sin(thetah)*thetadoth; -rfem*cos(thetah)*thetadoth];
vGtib=[lfem*sin(thetah)*thetadoth + rtib*sin(thetat)*thetadott;...
      -lfem*cos(thetah)*thetadoth - rtib*cos(thetat)*thetadott];
vGfoot=[lfem*sin(thetah)*thetadoth + ltib*sin(thetat)*thetadott + rfoot*sin(thetaf)*thetadotf;...
       -lfem*cos(thetah)*thetadoth - rtib*cos(thetat)*thetadott - rfoot*cos(thetaf)*thetadotf];

Tfem= simplify((1/2)*mfem*sum(vGfem.^2)+(1/2)*Ifem*(thetadoth.^2));
Ttib= simplify((1/2)*mtib*sum(vGtib.^2)+(1/2)*Itib*(thetadott.^2));
Tfoot=simplify((1/2)*mfoot*sum(vGfoot.^2)+(1/2)*Ifoot*(thetadotf.^2));

Ufem= -mfem*g*rfem*sin(thetah);
Utib=-mtib*g*(-lfem*sin(thetah)-rtib*sin(thetat));
Ufoot=-mfoot*g*(-lfem*sin(thetah)-ltib*sin(thetat)-rfoot*sin(thetaf)); %the derivative of this should be 0

T=Tfem+Ttib+Tfoot;
U=Ufem+Utib+Ufoot;

lambda=[lambda1; lambda2];

c=[lfem*sin(thetah)+ltib*sin(thetat)-h; thetaf];
L=T-U+lambda'*c;

symthetas={     thetad2doth,    thetad2dott,    thetad2dotf,    thetadoth,      thetadott,      thetadotf,      thetah,     thetat,     thetaf};
nosymthetas={   thetad2doth_aux,thetad2dott_aux,thetad2dotf_aux thetadoth_aux,  thetadott_aux,  thetadotf_aux,  thetah_aux, thetat_aux, thetaf_aux,};

L0=subs(L,symthetas,nosymthetas);

%% For thetah
dL0dthetadoth=simplify(diff(L0,thetadoth_aux));
dLdtthetadoth=subs(dL0dthetadoth,nosymthetas,symthetas);
dLdtdtthetadoth=simplify(diff(dLdtthetadoth,t));
dL0dtdtthetadoth=subs(dLdtdtthetadoth,symthetas,nosymthetas);

dL0dthetah=simplify(diff(L0,thetah_aux));

eq(1,1)=simplify(dL0dtdtthetadoth-dL0dthetah-Mh-Mk);

%% For thetat
dL0dthetadott=simplify(diff(L0,thetadott_aux));
dLdtthetadott=subs(dL0dthetadott,nosymthetas,symthetas);
dLdtdtthetadott=simplify(diff(dLdtthetadott,t));
dL0dtdtthetadott=subs(dLdtdtthetadott,symthetas,nosymthetas);

dL0dthetat=simplify(diff(L0,thetat_aux));

eq(2,1)=simplify(dL0dtdtthetadott-dL0dthetat+Mk+Ma);


%% For thetaa
dL0dthetadotf=simplify(diff(L0,thetadotf_aux));
dLdtthetadotf=subs(dL0dthetadotf,nosymthetas,symthetas);
dLdtdtthetadotf=simplify(diff(dLdtthetadotf,t));
dL0dtdtthetadotf=subs(dLdtdtthetadotf,symthetas,nosymthetas);

dL0dthetaf=simplify(diff(L0,thetaf_aux));

eq(3,1)=simplify(dL0dtdtthetadotf-dL0dthetaf-Ma);



Mh_steady=subs(sol.Mh,{thetad2doth_aux,thetad2dott_aux,thetad2dotf_aux,thetadoth_aux,thetadott_aux,thetadotf_aux,thetaf_aux},{zeros(1,7)});
