%% SPM Motor Design
%% 1.0 - Source: file Mathcad  SPM_design_6.mcdx prof. Di Gerlando
%% DATA


function [eta_NP,Volume]=function_SPM_motor(hm,bm,N_pc)



%Constants
k=1000;
mu_0=4*pi*10^(-7);

% Copper resistivity at 20 °C
ro_cu_20=0.0175*k^(-2); 
% Copper critical temperature
teta_c_cu=235;
% Reference resistivity temperature
teta_ref_st=200;

% PM reference temperature
teta_PM_ref=160;

% Copper resistivity at θ_ref_st 
ro_cu_rif_st=ro_cu_20*(teta_c_cu+teta_ref_st)/(teta_c_cu+20);

% Peak point torque
TP=130; %Nm 

% Peak point speed
NP=8815; %rpm

% Rated torque 
Tn=100;

% Rated torque speed 
NTn=7000;

% Rated speed
Nn=9545;

% DC bus voltage (min )
Vdc_min=420;

% Maximum speed
NM=17000;

% Stator Bore Diameter
D=0.130;


% Inverter modulation index at minimum Vdc, at peak point
ma_peak=0.95;

% Fundamental inverter Max voltage at peak point
VinvM=Vdc_min/(2*sqrt(2))*ma_peak; %Vrms

% Air-gap peak B
B1=0.86;

% Linear current density
DELTA=100*k;

% Reaction angle
gamma=0;

% Efficiency (just core and mechanical losses)
eta_m=0.98;

% Presumed winding factor, for single layer winding without pitch shortening:
k_w=0.966;

% PM material: N48UZ-SGR 
Br20C=1.37; 
mu_rec_pu=1.05;
kBr=-0.1; %/°C 
ro_PM=1.4*k^(-2);

% PM remanence Br at θPM 
Br=Br20C*(1+0.01*kBr*(teta_PM_ref-20));

% Linearly extrapolated coercitive force: ≔ Ho ⎛⎝θPM⎞⎠
Ho=Br*(mu_0*mu_rec_pu)^(-1);

% Lenght
L=TP/(eta_m*cos(gamma)*k_w*pi/(2*sqrt(2))*B1*DELTA*D^2);

% Poles
p=6;

% Pole pitch
tau=pi*D/p;

% Peak point and maximum frequency
f_p=p*NP/120;

f_m=p*NM/120;

% Non-magnetic bandage width 
w_b=3*k^(-1);

% Mechanical air-gap width
delta_m=1*k^(-1);

% Magnetic air-gap width
delta=w_b+delta_m;

% % Magnet height hm DA OTTIMIZZARE
% hm=0.006632; 
% 
% % bm DA OTTIMIZARE
% bm=0.056

% Volume of one PM
V1PM=hm*bm*L;


%Air-gap flux density
B_i=mu_0*Ho*hm/(delta+hm*(mu_rec_pu)^-1);

%Eta PM: Eta_BM0,K_BM
eta_BM0=0.993;
k_BM=2.955;
eta_PM=eta_BM0-k_BM*hm;

% Peak flux density B1 in the air-gap, due to PMs
B1mPM=eta_PM*4/pi*B_i*sin(bm/tau*pi/2);

% Average air-gap flux density of one pole
B_av=B_i*bm/tau;

% Flux density of one pole
Phi=B_av*tau*L;

% Fundamental harmonic pole flux (PM at teta_PMref)
Phi_1=B1mPM*2/pi*tau*L;

% Rotor yoke external diameter
D_ext_rotor_yoke=D-2*(delta+hm);

%%Winding data definition

% No of slots/(pole-phase)
q=2;

% Slot pitch
tau_s=tau/(3*q);

% No of stator slots Ns (single layer winding, full pitch)
Ns=3*p*q;

% Coil pitch
yc=3*q;

% Pitch Shortening
epsilon=3*q-yc;

% Winding factor
kw=1/(2*q*sin(pi/(6*q)));
omega_p=2*pi*f_p;

%No of parallel paths:
a=p/2;

% Single conductor emf
E_cp=pi/sqrt(2)*f_p*Phi_1;

% Presumed E/V ratio, to take into account voltage drop: (to be checked and updated)
ro_EV=0.72;

% Theoretical No of series connected conductors
U_th=ro_EV*VinvM/(kw*E_cp);


% No u of conductors/slot (= No of turns/coil)
u=round(3*U_th*a/Ns);
U=Ns*u/(3*a);

% Phase emf
E_p=E_cp*U*kw;

% PM flux linkage
Psi_PMT=kw*U*Phi_1/(2*sqrt(2));
Psi_PMT*omega_p;

% Peak current
I_pp=TP*NP*pi/30/(eta_m*3*E_p);

% Peak path current
I_p_path=I_pp/a;

% Stator lamination geometry definition

% stacking factor
k_st=0.95;

% L_fe
L_fe=L*k_st;

% Flux density in the tooth (at the minimum width):
B_d=1.75;

% Stator tooth minimum width
b_dm=Phi/(2/pi*B_d*Ns/p*L_fe);

% Slot width
bs=pi*D/Ns-b_dm;

% Slot opening width
b_as=2/k;

% Heigth of the tooth horn
h_as=2/k;

% Heigth of the horn at the base
h_td=3/k;

% Equivalent magnetic air-gap
g=hm/mu_rec_pu+delta;

% Slotting auxiliary function
sigma_s=2/pi*(atan(1/2*b_as/g)-g/b_as*log(1+(1/2*b_as/g)^2));

% Carter's factor
Kc=1/(1-sigma_s*b_as*tau_s^-1);

%% Definition of slot and conductors geometry:

% Slot filling factor
alfa_cu=0.78;

% Peak surface current density
S_p=15*k^2;

% Path conductor cross section:
A_co=I_p_path/S_p;

% Total copper cross section in slot
A_cu_slot=u*A_co;

% Slot cross section:
A_slot=A_cu_slot/alfa_cu;

% Tooth height (= slot height):
h_d=A_slot/bs;

% Plate enamel / liner width
w_e=0.05/k;
w_l=0.1/k;

% Theoretical value of the plate width
w_pth=bs-2*w_l-2*w_e-0.1/k;

% Theoretical total conductor height
h_totcth=A_co/w_pth;

%No of plates for each rectangular conductors (SE SERVE DA OTTIMIZZARE)
% N_pc

% w_p
w_p=round(w_pth,4);

% No of plates for each rectangular conductors FISSO A DUE
% N_pc=2;
h_p_1=1.12/k;
A_p_1=4.83/k^2;

% Penetration depth:
s_p=sqrt(ro_cu_rif_st/(pi*f_p*mu_0)*bs/w_p);

% Csi/phi/psi
Csi_f=h_p_1/s_p;
phi_f=Csi_f*(sinh(2*Csi_f)+sin(2*Csi_f))/(cosh(2*Csi_f)-cos(2*Csi_f));
psi_f=2*Csi_f*(sinh(Csi_f)-sin(Csi_f))/(cosh(Csi_f)+cos(Csi_f));


% Estimated endwinding length
Lt=pi*(D+2*(h_as+h_d))/Ns*(yc-1)*pi/2;

% ohmic phase resistance
R_1_ohm_v=ro_cu_rif_st*(U*(L+Lt))/(a*A_p_1*N_pc);

N_pr=u*N_pc;

A_sf=(N_pr^2-1)/3;
beta_cc=sqrt(L/(L+Lt)*(1+2*w_e/h_p_1));

% ka contributions:       skin effect,  proximity effect
ka_Npc=1+L/(L+Lt)*(phi_f+A_sf*psi_f-1);

% Recalculation of the slot heigth
h_d=u*N_pc*(h_p_1+2*w_e)+2*w_l+0.1/k;

% Flux density in the stator yoke:
B_y=1.6;

% Actual slot fill factor
alfa_cu=(u*N_pc*A_p_1)/(h_d*bs);

% Stator yoke width
h_y=Phi/(2*L_fe*B_y);

% Stator core external diameter
D_se=D+2*(h_td+h_d+h_y);


% Equivalent circuit parameters evaluation:

% Equivalent air-gap
g_e=g*Kc;

% Main field specific permeance:
lambda_p=3/pi^2*k_w^2*mu_0*tau/g_e;

% Reaction inductance
L_ro=U^2/p*L*lambda_p;

% Slot leakage specific permeance
lambda_s=mu_0*((h_d-2*w_l)/(3*bs)+w_l/bs+(h_td-h_as)/(0.5*(bs+b_as))+h_as/b_as);

% Harmonic leakage coefficient
sigma_a=(5*q^2+1+epsilon^2*(4*q)^-1-1.5*epsilon^2-epsilon*(4*q)^-1)/(6*(3/pi)^2*k_w^2*q^2)-1;

% Harmonic leakage specific permeance
lambda_a=sigma_a*lambda_p;

% Air-gap specific permeance
lambda_d=mu_0*g/(b_as+0.8*g);

% Endwinding specific permeance
lambda_ew=0.4/k^2;

% Total leakage specific permeance
lambda=lambda_s/q+lambda_a+lambda_d/q+lambda_ew*Lt/L;

% Leakage inductance
L_L=U^2/p*L*lambda;

% Synchronous inductance
L_s=L_ro+L_L;
X_p=omega_p*L_s;

% Reactive voltage drop
X_p*I_pp;
Z_p=sqrt(R_1_ohm_v^2+(omega_p*L_s)^2);
X_p/Z_p;

% Actual peak point  ro_EV ratio ( gamma = 0)
V_p=sqrt((E_p+R_1_ohm_v*I_pp)^2+(omega_p*L_s*I_pp)^2);

% Short circuit current
I_sc=Psi_PMT/L_s;



% Losses estimation (at teta_cu = teta_cur_rif   and teta_PM = teta_PM_rif )

% Stator copper losses:
%Stator current:
I_TP=TP/(eta_m*3*Psi_PMT*p/2);

% Stator Cu losses
P_L_Cu=3*R_1_ohm_v*I_TP^2;

%Core losses

% Data
% Magnetic material: M235-35A (lamination width 0.35 mm, with 2.35 W/kg at 50 Hz, 1.5 T) 
gamma_lam=7650;
B_LE=0.1:0.1:2.5;
H_LE=[0 24.7 32.6 38.1 43.1 48.2 53.9 60.7 68.8 79.3 93.7 115 156 260 690 1950 4410 7630 12000 18200 26900 38700 54300 74400 99700 130800];

% Modified Steinmetz method is adopted, for which the specific core losses equal:
alfa_St=2.07;
beta_St=0;
K_eddy_St=2.199*10^-6;
K_h_St=0.02533;

% Resultant flux density in the air-gap
B_a_NP=mu_0/g*3*sqrt(2)/pi*k_w*U/p*I_TP;

B_gap_tot=abs(B1mPM+B_a_NP*exp(i*(pi/2+gamma)));

Phi_load=2/pi*B_gap_tot*tau*L;

% Tooth width bd at 1/3 of the min width:

b_dM=pi*(D+2*h_d)/Ns-bs;
b_d=b_dm+1/3*(b_dM-b_dm);

% Tooth and stator yoke flux density:

B_d_av=Phi_load/(2/pi*Ns/p*b_d*L*k_st);
B_y_st=Phi_load/(2*L_fe*h_y);

% Stator core masses:
M_teeth=gamma_lam*(pi/4*((D+2*h_d)^2-D^2)-Ns*h_d*bs)*L_fe;

M_st_y=gamma_lam*h_y*pi*(D+2*h_d+h_y)*L_fe;

% Presumed additional core loss coefficient (punching effects, field non-uniformity):
k_add_fe=1.15;

f=f_p;

P_fe_sp_id_av=K_h_St*f*B_d_av^(alfa_St+beta_St*B_d_av)+2*pi^2*K_eddy_St*f^2*B_d_av^2;
P_fe_sp_id_st=K_h_St*f*B_y_st^(alfa_St+beta_St*B_y_st)+2*pi^2*K_eddy_St*f^2*B_y_st^2;

% Losses in teeth and stator yoke
P_teeth=k_add_fe*P_fe_sp_id_av*M_teeth;
P_st_yoke=k_add_fe*P_fe_sp_id_st*M_st_y;

P_L_st_core=P_teeth+P_st_yoke;

% LOSSES IN THE PM
%Two mechanisms are considered: 
%1) slot opening modulation effect on the PM flux density (field model by straight lines and arcs);
%2) armature harmonic mmfs, moving with respect to the rotor

Ho_PM=8.929*10^5;
b_th=tau_s-b_as;
ni=2/pi;

x1=0:0.0001:b_th/2;
Bm1=0*x1+mu_0*Ho_PM*hm/(delta+hm*mu_rec_pu^-1);
x2=b_th/2+1e-6:0.0001:b_th/2+b_as/2;
Bm2=mu_0*Ho_PM*hm./(delta+hm*mu_rec_pu^-1+ni*pi/2*(x2-b_th/2));
x3=b_th/2+1e-6+b_as/2:0.0001:b_th/2+b_as;
Bm3=mu_0*Ho_PM*hm./(delta+hm*mu_rec_pu^-1+ni*pi/2*(-x3+b_as+b_th/2));
x4=b_th/2+1e-6+b_as:0.0001:b_th+b_as;
Bm4=0*x4+mu_0*Ho_PM*hm/(delta+hm*mu_rec_pu^-1);

xd=[x1 x2 x3 x4];
BMd=[Bm1 Bm2 Bm3 Bm4];


x=0:0.002*tau_s:tau_s;
BM=interp1(xd,BMd,x);


B_mo=1/tau_s*trapz(x,BM);


for h_m=1:15
B_m_(h_m)=2/tau_s*trapz(x,BM.*cos(h_m.*x./tau_s*2*pi));
end

B_mF=B_mo;
for h_m=1:15
B_mF=B_mF+(B_m_(h_m)*cos(h_m*x/tau_s*2*pi));
end
B_mac=B_mF-B_mo;

% Rotor mechanical speed:
Omega_p=omega_p*2/p;

% Rotor position (measured along the stator bore, from a slot axis):

period_p=(1/2/pi*omega_p)^(-1);

t=0:0.0005*period_p:0.1*period_p;

z_hP=0.5*D*Omega_p*t;

% Average radius within PM
R_q=0.5*D-delta-hm;

% Speed vq of the point q in the PM:
v_q=Omega_p*R_q;

B_mac=interp1(x,B_mac,z_hP);

% Electric field induced by slotting (emf per unit length is axial direction):
K_ec_slot=B_mac*v_q;

% ec current density distribution in the PM
Sec_slot=K_ec_slot*ro_PM^(-1);

% eddy current losses per unit PM volume: 
Losses_PM_volume_slot_t=ro_PM*Sec_slot.^2;

% time average specific PM losses
Losses_PM_vol_slot_av=1/max(t)*trapz(t,Losses_PM_volume_slot_t);

% Total losses in the PMs due to slotting effect:
P_L_PM_slot=p*bm*hm*L*Losses_PM_vol_slot_av;

% armature harmonic mmfs, moving with respect to the rotor
k_hm=8;
k_h=-k_hm:k_hm;
h=1+k_h*3*q;

k_wh=abs(sin(h*pi/6).*(q*sin(h*pi/6/q)).^(-1));

% air-gap peak flux density for the harmonic order h:

B_ah_NP=B_a_NP.*k_wh.*1./(abs(h));
Omega_h=Omega_p*(1-1./h);
b_m_tau=bm/tau;

N_me=20;
b_PMe=bm/N_me;

ro_av=2/pi*1/b_m_tau*N_me*abs(sin(h*b_m_tau*1/N_me*pi/2)./h);

B_ach_NP=B_ah_NP.*(1-ro_av);
K_ach=R_q.*Omega_h.*B_ach_NP;
S_ec_h=K_ach.*ro_PM^(-1);

P_L_PM_h=p*bm*hm*L*sum(1/2*ro_PM*S_ec_h.^2);

% Mechanical losses:
% Presumed mechanical losses in rated conditions;
% linear dependence on speed
P_L_mecc_n=600;
P_L_mecc_NP=P_L_mecc_n*NP/Nn;

% Total Losses at NP
P_L_tot_NP=P_L_Cu+P_L_st_core+P_L_PM_h+P_L_mecc_NP;

% Output power
P_out_NP=TP*NP*pi/30;

% Efficiency
eta_NP=P_out_NP/(P_out_NP+P_L_tot_NP);

%Volume 
Volume=(L+Lt)*pi*D_se^2/4*1000000000;


