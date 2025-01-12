%% Foot
m_foot=1.93e-3;
%assuming cylinder of 15 mm long by 1.5 mm radius...
%Iy' and Iz' equal, being Iy' perpendicular to the "plane" of the foot amd
%Iz' perp to the sagittal plane
r_foot=1.5e-3;
l_foot=15e-3;
Iyp_foot=(1/4)*m_foot*(r_foot^2)+(1/12)*m_foot*(l_foot^2);
Izp_foot=Iyp_foot;
Ixp_foot=(1/2)*m_foot*(r_foot^2);
S=[cosd(30) -sind(30) 0; sind(30) cosd(30) 0; 0 0 1]; %aprox rotation according to Johnson model
IIp=[Ixp_foot 0 0; 0 Iyp_foot 0; 0 0 Izp_foot];
II=S*IIp*S';

%% Tibia
m_tibia=0.00538;
%assuming cylinder of 42.6371 mm long by 2 mm radius
r_tibia=2e-3;
l_tibia=42.6371e-3;
Iy_tibia=(1/2)*m_tibia*(r_tibia^2);
Ix_tibia=(1/4)*m_tibia*(r_tibia^2)+(1/12)*m_tibia*(l_tibia^2);
Iz_tibia=Ix_tibia;
II_tibia=[Ix_tibia 0 0; 0 Iy_tibia 0; 0 0 Iz_tibia];

%%Femur
m_femur=13.51e-3;
%assuming cylinder of 30.8677 mm long by 2 mm radius
r_femur=2e-3;
l_femur=30.8677e-3;
Iy_femur=(1/2)*m_femur*(r_femur^2);
Ix_femur=(1/4)*m_femur*(r_femur^2)+(1/12)*m_femur*(l_femur^2);
Iz_femur=Ix_femur;
II_femur=[Ix_femur 0 0; 0 Iy_femur 0; 0 0 Iz_femur];

