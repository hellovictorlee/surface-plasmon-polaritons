clc;
clear all;
close all;

nx=299;
ny=199;
ds=1e-8;  
ep=8.85e-12;
mu=pi*4e-7;
c=1/sqrt(ep*mu);
n=0.5;   %stablility
dt=(n*(2^(-1/2))*ds)/c;
k=sqrt(mu/ep);

i0=15;
i1=nx-15;
j0=15;
j1=ny-15;

Exi=zeros(ny,nx-1);   % incident field
Eyi=zeros(ny-1,nx);
Hzxi=zeros(ny-1,nx-1);
Hzyi=zeros(ny-1,nx-1);
Hzi=Hzxi+Hzyi;

Ex=zeros(ny,nx-1);   % TFSF
Ey=zeros(ny-1,nx);
Hzx=zeros(ny-1,nx-1);
Hzy=zeros(ny-1,nx-1);
Hz=Hzx+Hzy;

eparray=ep*ones(ny,nx);
muarray=mu*ones(ny-1,nx-1);
conductivityE=zeros(ny,nx);
conductivityH=zeros(ny-1,nx-1);



Cai=(1 - (dt*conductivityE ./ (2*eparray)))./(1 + ((dt*conductivityE) ./ (2*eparray)));
Dai=(1 - (dt*conductivityH ./ (2*muarray)))./(1 + ((dt*conductivityH) ./ (2*muarray)));
Cbi=(dt ./ eparray) ./ (1 + ((dt*conductivityE) ./ (2*eparray)));
Dbi=(dt ./ muarray) ./ (1 + ((dt*conductivityH) ./ (2*muarray)));



%--------------------------------------------------------CPML
d=10;

psi_HzLi=zeros(ny-1,d);   % Left
psi_EyLi=zeros(ny-1,d);

psi_HzRi=zeros(ny-1,d);   % Right
psi_EyRi=zeros(ny-1,d);

psi_HzDi=zeros(d,nx-1);   % Down
psi_ExDi=zeros(d,nx-1);

psi_HzUi=zeros(d,nx-1);   % Up
psi_ExUi=zeros(d,nx-1);

psi_HzL=zeros(ny-1,d);   % Left
psi_EyL=zeros(ny-1,d);

psi_HzR=zeros(ny-1,d);   % Right
psi_EyR=zeros(ny-1,d);

psi_HzD=zeros(d,nx-1);   % Down
psi_ExD=zeros(d,nx-1);

psi_HzU=zeros(d,nx-1);   % Up
psi_ExU=zeros(d,nx-1);


m=3.5;
cm=(0.8*(m+1))/(k*ds);  % conductivity maximum value
%cm=0.7*(m+1)/(150*pi*dx);
D_H=(d:-1:1)';
L_H=(d:-1:1);
D_E=(d-0.5:-1:0.5)';
L_E=(d-0.5:-1:0.5);

segmaE_L=cm*((L_E./d).^m);
segmaE_L=repmat(segmaE_L,ny-1,1);
segmaH_L=cm*((L_H./d).^m);
segmaH_L=repmat(segmaH_L,ny-1,1);

segmaE_D=cm*((D_E./d).^m);
segmaE_D=repmat(segmaE_D,1,nx-1);
segmaH_D=cm*((D_H./d).^m);
segmaH_D=repmat(segmaH_D,1,nx-1);



km=1.1;   % kapa maximun
kapaE_L=1+(km-1)*((L_E./d).^m);
kapaE_L=repmat(kapaE_L,ny-1,1);
kapaH_L=1+(km-1)*((L_H./d).^m);
kapaH_L=repmat(kapaH_L,ny-1,1);
kapaE_R=fliplr(kapaE_L);
kapaH_R=fliplr(kapaH_L);

kapaE_D=1+(km-1)*((D_E./d).^m);
kapaE_D=repmat(kapaE_D,1,nx-1);
kapaH_D=1+(km-1)*((D_H./d).^m);
kapaH_D=repmat(kapaH_D,1,nx-1);
kapaE_U=flipud(kapaE_D);
kapaH_U=flipud(kapaH_D);



am=0;   % a maximun
ma=1;
aE_L=am*((1-(L_E./d)).^ma);
aE_L=repmat(aE_L,ny-1,1);
aH_L=am*((1-(L_H./d)).^ma);
aH_L=repmat(aH_L,ny-1,1);

aE_D=am*((1-(D_E./d)).^ma);
aE_D=repmat(aE_D,1,nx-1);
aH_D=am*((1-(D_H./d)).^ma);
aH_D=repmat(aH_D,1,nx-1);



be_L=exp(-((segmaE_L./kapaE_L)+aE_L)*(dt/ep));
be_R=fliplr(be_L);
be_D=exp(-((segmaE_D./kapaE_D)+aE_D)*(dt/ep));
be_U=flipud(be_D);


bh_L=exp(-((segmaH_L./kapaH_L)+aH_L)*(dt/ep));
bh_R=fliplr(bh_L);
bh_D=exp(-((segmaH_D./kapaH_D)+aH_D)*(dt/ep));
bh_U=flipud(bh_D);


ce_L=((segmaE_L./kapaE_L)./(segmaE_L+(aE_L.*kapaE_L))).*(be_L-1);
ce_R=fliplr(ce_L);
ce_D=((segmaE_D./kapaE_D)./(segmaE_D+(aE_D.*kapaE_D))).*(be_D-1);
ce_U=flipud(ce_D);


ch_L=((segmaH_L./kapaH_L)./(segmaH_L+(aH_L.*kapaH_L))).*(bh_L-1);
ch_R=fliplr(ch_L);
ch_D=((segmaH_D./kapaH_D)./(segmaH_D+(aH_D.*kapaH_D))).*(bh_D-1);
ch_U=flipud(ch_D);
%--------------------------------------------------------CPML



%--------------------Gaussian beam
H=1;
wavelength=2.5e-7;
freq=c/wavelength;
omega=2*pi*freq;

y_position=100;
x_position=12;
y_wing=60;
H_source=zeros(2*y_wing+1,1);

for j=1:2*y_wing+1
    
    H_source(j,1) = H*exp(-(1/2)*(3*(j-y_wing-1)/y_wing)^2);
    
end
%--------------------Gaussian beam





%-----------------------Drude model
J_p_Ex=zeros(ny,nx-1);
J_p_Ey=zeros(ny-1,nx);
r_p=zeros(ny,nx);
omega_p=zeros(ny,nx);


%-----circle
r=10;   % radian
xcen_1=50;
xcen_2=75;
xcen_3=100;
xcen_4=125;
xcen_5=150;
xcen_6=175;
xcen_7=200;
xcen_8=225;
xcen_9=250;

ycen_1=100;


out=ccc_eparray(eparray,ycen_1,xcen_1,r);   % circle_1
eparray=out;
out=ccc_r_p(r_p,ycen_1,xcen_1,r);
r_p=out;
out=ccc_omega_p(omega_p,ycen_1,xcen_1,r);
omega_p=out;

out=ccc_eparray(eparray,ycen_1,xcen_2,r);   % circle_2
eparray=out;
out=ccc_r_p(r_p,ycen_1,xcen_2,r);
r_p=out;
out=ccc_omega_p(omega_p,ycen_1,xcen_2,r);
omega_p=out;

out=ccc_eparray(eparray,ycen_1,xcen_3,r);   % circle_3
eparray=out;
out=ccc_r_p(r_p,ycen_1,xcen_3,r);
r_p=out;
out=ccc_omega_p(omega_p,ycen_1,xcen_3,r);
omega_p=out;

out=ccc_eparray(eparray,ycen_1,xcen_4,r);   % circle_4
eparray=out;
out=ccc_r_p(r_p,ycen_1,xcen_4,r);
r_p=out;
out=ccc_omega_p(omega_p,ycen_1,xcen_4,r);
omega_p=out;

out=ccc_eparray(eparray,ycen_1,xcen_5,r);   % circle_5
eparray=out;
out=ccc_r_p(r_p,ycen_1,xcen_5,r);
r_p=out;
out=ccc_omega_p(omega_p,ycen_1,xcen_5,r);
omega_p=out;

out=ccc_eparray(eparray,ycen_1,xcen_6,r);   % circle_6
eparray=out;
out=ccc_r_p(r_p,ycen_1,xcen_6,r);
r_p=out;
out=ccc_omega_p(omega_p,ycen_1,xcen_6,r);
omega_p=out;

out=ccc_eparray(eparray,ycen_1,xcen_7,r);   % circle_7
eparray=out;
out=ccc_r_p(r_p,ycen_1,xcen_7,r);
r_p=out;
out=ccc_omega_p(omega_p,ycen_1,xcen_7,r);
omega_p=out;

out=ccc_eparray(eparray,ycen_1,xcen_8,r);   % circle_8
eparray=out;
out=ccc_r_p(r_p,ycen_1,xcen_8,r);
r_p=out;
out=ccc_omega_p(omega_p,ycen_1,xcen_8,r);
omega_p=out;

out=ccc_eparray(eparray,ycen_1,xcen_9,r);   % circle_9
eparray=out;
out=ccc_r_p(r_p,ycen_1,xcen_9,r);
r_p=out;
out=ccc_omega_p(omega_p,ycen_1,xcen_9,r);
omega_p=out;
%-----circle


k_p    = (1-r_p*dt/2)       ./ (1+r_p*dt/2);
beta_p = (omega_p.^2*ep*dt/2) ./ (1+r_p*dt/2);
%-----------------------Drude model



Ca=(2*eparray - dt*beta_p - conductivityE*dt) ./ (2*eparray + dt*beta_p + conductivityE*dt);
Cb=(2*dt)                                     ./ (2*eparray + dt*beta_p + conductivityE*dt);
Da=(1-(dt*conductivityH./(2*muarray)))        ./ (1+ dt*conductivityH./(2*muarray));
Db=(dt./(muarray))                            ./ (1+ dt*conductivityH./(2*muarray));




for t=1:5000
    
    
    
    %------------------------------------------------------Eyi
    psi_EyLi = be_L .* psi_EyLi + ce_L .* (Hzi(:,2:d+1)     - Hzi(:,1:d))         / ds;   % Left
    psi_EyRi = be_R .* psi_EyRi + ce_R .* (Hzi(:,nx-d:nx-1) - Hzi(:,nx-d-1:nx-2)) / ds;   % Right
    
    
    Eyi(:,d+2:nx-d-1) = Cai(1:ny-1,d+2:nx-d-1) .* Eyi(:,d+2:nx-d-1) - Cbi(1:ny-1,d+2:nx-d-1) .* ((Hzi(:,d+2:nx-d-1) - Hzi(:,d+1:nx-d-2))   / ds);   % middle
    Eyi(:,2:d+1)      = Cai(1:ny-1,2:d+1)      .* Eyi(:,2:d+1)      - Cbi(1:ny-1,2:d+1)      .* ((Hzi(:,2:d+1)      - Hzi(:,1:d))         ./ (kapaE_L * ds));   % Left
    Eyi(:,nx-d:nx-1)  = Cai(1:ny-1,nx-d:nx-1)  .* Eyi(:,nx-d:nx-1)  - Cbi(1:ny-1,nx-d:nx-1)  .* ((Hzi(:,nx-d:nx-1)  - Hzi(:,nx-d-1:nx-2)) ./ (kapaE_R * ds));   % Right
    
    
    Eyi(:,2:d+1)     = Eyi(:,2:d+1)     - Cbi(1:ny-1,2:d+1)     .* psi_EyLi;   % Left
    Eyi(:,nx-d:nx-1) = Eyi(:,nx-d:nx-1) - Cbi(1:ny-1,nx-d:nx-1) .* psi_EyRi;   % Right
    %------------------------------------------------------Eyi
    %------------------------------------------------------Exi
    psi_ExDi = be_D .* psi_ExDi + ce_D .* (Hzi(2:d+1,:)     - Hzi(1:d,:))         / ds;   % Down
    psi_ExUi = be_U .* psi_ExUi + ce_U .* (Hzi(ny-d:ny-1,:) - Hzi(ny-d-1:ny-2,:)) / ds;   % Up 
    
    
    Exi(d+2:ny-d-1,:) = Cai(d+2:ny-d-1,1:nx-1) .* Exi(d+2:ny-d-1,:) + Cbi(d+2:ny-d-1,1:nx-1) .* ((Hzi(d+2:ny-d-1,:) - Hzi(d+1:ny-d-2,:))   / ds);   % middle 
    Exi(2:d+1,:)      = Cai(2:d+1,1:nx-1)      .* Exi(2:d+1,:)      + Cbi(2:d+1,1:nx-1)      .* ((Hzi(2:d+1,:)      - Hzi(1:d,:))         ./ (kapaE_D * ds));   % Down
    Exi(ny-d:ny-1,:)  = Cai(ny-d:ny-1,1:nx-1)  .* Exi(ny-d:ny-1,:)  + Cbi(ny-d:ny-1,1:nx-1)  .* ((Hzi(ny-d:ny-1,:)  - Hzi(ny-d-1:ny-2,:)) ./ (kapaE_U * ds));   % Up
    
    
    Exi(2:d+1,:)     = Exi(2:d+1,:)     + Cbi(2:d+1,1:nx-1)     .* psi_ExDi;   % Down   
    Exi(ny-d:ny-1,:) = Exi(ny-d:ny-1,:) + Cbi(ny-d:ny-1,1:nx-1) .* psi_ExUi;   % Up
    %------------------------------------------------------Exi
    %------------------------------------------------------Hzi
    psi_HzLi = bh_L .* psi_HzLi + ch_L .* (Eyi(:,2:d+1)     - Eyi(:,1:d))       / ds;   % Left
    psi_HzRi = bh_R .* psi_HzRi + ch_R .* (Eyi(:,nx-d+1:nx) - Eyi(:,nx-d:nx-1)) / ds;   % Right
    psi_HzDi = bh_D .* psi_HzDi + ch_D .* (Exi(2:d+1,:)     - Exi(1:d,:))       / ds;   % Down
    psi_HzUi = bh_U .* psi_HzUi + ch_U .* (Exi(ny-d+1:ny,:) - Exi(ny-d:ny-1,:)) / ds;   % Up
    
    
    Hzi(:,d+1:nx-d-1) = Dai(:,d+1:nx-d-1) .* Hzi(:,d+1:nx-d-1) - Dbi(:,d+1:nx-d-1) .* ((Eyi(:,d+2:nx-d)  - Eyi(:,d+1:nx-d-1)) / ds);   % middle                        
    Hzi(:,1:d)        = Dai(:,1:d)        .* Hzi(:,1:d)        - Dbi(:,1:d)        .* ((Eyi(:,2:d+1)     - Eyi(:,1:d))       ./ (kapaH_L * ds));   % Left
    Hzi(:,nx-d:nx-1)  = Dai(:,nx-d:nx-1)  .* Hzi(:,nx-d:nx-1)  - Dbi(:,nx-d:nx-1)  .* ((Eyi(:,nx-d+1:nx) - Eyi(:,nx-d:nx-1)) ./ (kapaH_R * ds));   % Right
    Hzi(d+1:ny-d-1,:) = Dai(d+1:ny-d-1,:) .* Hzi(d+1:ny-d-1,:) + Dbi(d+1:ny-d-1,:) .* ((Exi(d+2:ny-d,:)  - Exi(d+1:ny-d-1,:)) / ds);   % middle
    Hzi(1:d,:)        = Dai(1:d,:)        .* Hzi(1:d,:)        + Dbi(1:d,:)        .* ((Exi(2:d+1,:)     - Exi(1:d,:))       ./ (kapaH_D * ds));   % Down
    Hzi(ny-d:ny-1,:)  = Dai(ny-d:ny-1,:)  .* Hzi(ny-d:ny-1,:)  + Dbi(ny-d:ny-1,:)  .* ((Exi(ny-d+1:ny,:) - Exi(ny-d:ny-1,:)) ./ (kapaH_U * ds));   % Up
    
    
    Hzi(:,1:d)       = Hzi(:,1:d)       - Dbi(:,1:d)       .* psi_HzLi;   % Left
    Hzi(:,nx-d:nx-1) = Hzi(:,nx-d:nx-1) - Dbi(:,nx-d:nx-1) .* psi_HzRi;   % Right
    Hzi(1:d,:)       = Hzi(1:d,:)       + Dbi(1:d,:)       .* psi_HzDi;   % Down
    Hzi(ny-d:ny-1,:) = Hzi(ny-d:ny-1,:) + Dbi(ny-d:ny-1,:) .* psi_HzUi;   % Up
    %------------------------------------------------------Hzi
    
    
    
    
    
    
    
    
    
    
    
    %------------------------------------------------------Ey    
    Ey_temp = Ey;
    
    
    psi_EyL = be_L .* psi_EyL + ce_L .* (Hz(:,2:d+1)     - Hz(:,1:d))         / ds;   % Left
    psi_EyR = be_R .* psi_EyR + ce_R .* (Hz(:,nx-d:nx-1) - Hz(:,nx-d-1:nx-2)) / ds;   % Right
    
    
    Ey(:,2:d+1)      = Ca(1:ny-1,2:d+1)      .* Ey(:,2:d+1)      - Cb(1:ny-1,2:d+1)      .* ((Hz(:,2:d+1)      - Hz(:,1:d))         ./ (kapaE_L * ds) + (1/2)*(1+k_p(1:ny-1,2:d+1))      .* J_p_Ey(:,2:d+1));   % Left
    Ey(:,d+2:nx-d-1) = Ca(1:ny-1,d+2:nx-d-1) .* Ey(:,d+2:nx-d-1) - Cb(1:ny-1,d+2:nx-d-1) .* ((Hz(:,d+2:nx-d-1) - Hz(:,d+1:nx-d-2))   / ds             + (1/2)*(1+k_p(1:ny-1,d+2:nx-d-1)) .* J_p_Ey(:,d+2:nx-d-1));   % middle
    Ey(:,nx-d:nx-1)  = Ca(1:ny-1,nx-d:nx-1)  .* Ey(:,nx-d:nx-1)  - Cb(1:ny-1,nx-d:nx-1)  .* ((Hz(:,nx-d:nx-1)  - Hz(:,nx-d-1:nx-2)) ./ (kapaE_R * ds) + (1/2)*(1+k_p(1:ny-1,nx-d:nx-1))  .* J_p_Ey(:,nx-d:nx-1));   % Right
    
    
    Ey(:,2:d+1)     = Ey(:,2:d+1)     - Cb(1:ny-1,2:d+1)     .* psi_EyL;   % Left
    Ey(:,nx-d:nx-1) = Ey(:,nx-d:nx-1) - Cb(1:ny-1,nx-d:nx-1) .* psi_EyR;   % Right
    
    
    J_p_Ey = k_p(1:ny-1,1:nx) .* J_p_Ey + beta_p(1:ny-1,1:nx) .* (Ey+Ey_temp);
    %------------------------------------------------------Ey  
    
    
    
    %------------------------------------------------------Ex
    Ex_temp = Ex;
    
    
    psi_ExD = be_D .* psi_ExD + ce_D .* (Hz(2:d+1,:)     - Hz(1:d,:))         / ds;   % Down
    psi_ExU = be_U .* psi_ExU + ce_U .* (Hz(ny-d:ny-1,:) - Hz(ny-d-1:ny-2,:)) / ds;   % Up  
    
    
    Ex(2:d+1,:)      = Ca(2:d+1,1:nx-1)      .* Ex(2:d+1,:)      + Cb(2:d+1,1:nx-1)      .* ((Hz(2:d+1,:)      - Hz(1:d,:))         ./ (kapaE_D * ds) - (1/2)*(1+k_p(2:d+1,1:nx-1))      .*J_p_Ex(2:d+1,:));    % Down
    Ex(d+2:ny-d-1,:) = Ca(d+2:ny-d-1,1:nx-1) .* Ex(d+2:ny-d-1,:) + Cb(d+2:ny-d-1,1:nx-1) .* ((Hz(d+2:ny-d-1,:) - Hz(d+1:ny-d-2,:))   / ds             - (1/2)*(1+k_p(d+2:ny-d-1,1:nx-1)) .*J_p_Ex(d+2:ny-d-1,:));    % middle 
    Ex(ny-d:ny-1,:)  = Ca(ny-d:ny-1,1:nx-1)  .* Ex(ny-d:ny-1,:)  + Cb(ny-d:ny-1,1:nx-1)  .* ((Hz(ny-d:ny-1,:)  - Hz(ny-d-1:ny-2,:)) ./ (kapaE_U * ds) - (1/2)*(1+k_p(ny-d:ny-1,1:nx-1))  .*J_p_Ex(ny-d:ny-1,:));    % Up
    
    
    Ex(2:d+1,:)     = Ex(2:d+1,:)     + Cb(2:d+1,1:nx-1)     .* psi_ExD;   % Down   
    Ex(ny-d:ny-1,:) = Ex(ny-d:ny-1,:) + Cb(ny-d:ny-1,1:nx-1) .* psi_ExU;   % Up
    
    
    J_p_Ex = k_p(1:ny,1:nx-1) .* J_p_Ex + beta_p(1:ny,1:nx-1) .* (Ex+Ex_temp);
    %------------------------------------------------------Ex
    
    
    
    %------------------------------------------------------Hz
    psi_HzL = bh_L .* psi_HzL + ch_L .* (Ey(:,2:d+1)     - Ey(:,1:d))       / ds;   % Left
    psi_HzR = bh_R .* psi_HzR + ch_R .* (Ey(:,nx-d+1:nx) - Ey(:,nx-d:nx-1)) / ds;   % Right
    psi_HzD = bh_D .* psi_HzD + ch_D .* (Ex(2:d+1,:)     - Ex(1:d,:))       / ds;   % Down
    psi_HzU = bh_U .* psi_HzU + ch_U .* (Ex(ny-d+1:ny,:) - Ex(ny-d:ny-1,:)) / ds;   % Up
    
    
    Hz(:,d+1:nx-d-1) = Da(:,d+1:nx-d-1) .* Hz(:,d+1:nx-d-1) - Db(:,d+1:nx-d-1) .* ((Ey(:,d+2:nx-d)  - Ey(:,d+1:nx-d-1)) / ds);   % middle                        
    Hz(:,1:d)        = Da(:,1:d)        .* Hz(:,1:d)        - Db(:,1:d)        .* ((Ey(:,2:d+1)     - Ey(:,1:d))       ./ (kapaH_L * ds));   % Left
    Hz(:,nx-d:nx-1)  = Da(:,nx-d:nx-1)  .* Hz(:,nx-d:nx-1)  - Db(:,nx-d:nx-1)  .* ((Ey(:,nx-d+1:nx) - Ey(:,nx-d:nx-1)) ./ (kapaH_R * ds));   % Right
    Hz(d+1:ny-d-1,:) = Da(d+1:ny-d-1,:) .* Hz(d+1:ny-d-1,:) + Db(d+1:ny-d-1,:) .* ((Ex(d+2:ny-d,:)  - Ex(d+1:ny-d-1,:)) / ds);   % middle
    Hz(1:d,:)        = Da(1:d,:)        .* Hz(1:d,:)        + Db(1:d,:)        .* ((Ex(2:d+1,:)     - Ex(1:d,:))       ./ (kapaH_D * ds));   % Down
    Hz(ny-d:ny-1,:)  = Da(ny-d:ny-1,:)  .* Hz(ny-d:ny-1,:)  + Db(ny-d:ny-1,:)  .* ((Ex(ny-d+1:ny,:) - Ex(ny-d:ny-1,:)) ./ (kapaH_U * ds));   % Up
    
    
    Hz(:,1:d)       = Hz(:,1:d)       - Db(:,1:d)       .* psi_HzL;   % Left
    Hz(:,nx-d:nx-1) = Hz(:,nx-d:nx-1) - Db(:,nx-d:nx-1) .* psi_HzR;   % Right
    Hz(1:d,:)       = Hz(1:d,:)       + Db(1:d,:)       .* psi_HzD;   % Down
    Hz(ny-d:ny-1,:) = Hz(ny-d:ny-1,:) + Db(ny-d:ny-1,:) .* psi_HzU;   % Up
    %------------------------------------------------------Hz
    
    
    
    
    
    
    Hz(j0,i0:i1)   = Da(j0,i0:i1)   .* Hz(j0,i0:i1)   - Db(j0,i0:i1)   .* Exi(j0,i0:i1)   / ds;   % j=j0 Face
    Hz(j1,i0:i1)   = Da(j1,i0:i1)   .* Hz(j1,i0:i1)   + Db(j1,i0:i1)   .* Exi(j1+1,i0:i1) / ds;   % j=j1 Face
    Hz(j0:j1,i0)   = Da(j0:j1,i0)   .* Hz(j0:j1,i0)   + Db(j0:j1,i0)   .* Eyi(j0:j1,i0)   / ds;   % i=i0 Face
    Hz(j0:j1,i1)   = Da(j0:j1,i1)   .* Hz(j0:j1,i1)   - Db(j0:j1,i1)   .* Eyi(j0:j1,i1+1) / ds;   % i=i1 Face
    Ex(j0,i0:i1)   = Ca(j0,i0:i1)   .* Ex(j0,i0:i1)   - Cb(j0,i0:i1)   .* Hzi(j0,i0:i1)   / ds;   % j=j0-1/2 Face
    Ex(j1+1,i0:i1) = Ca(j1+1,i0:i1) .* Ex(j1+1,i0:i1) + Cb(j1+1,i0:i1) .* Hzi(j1,i0:i1)   / ds;   % j=j1+1/2 Face
    Ey(j0:j1,i0)   = Ca(j0:j1,i0)   .* Ey(j0:j1,i0)   + Cb(j0:j1,i0)   .* Hzi(j0:j1,i0)   / ds;   % i=i0-1/2 Face
    Ey(j0:j1,i1+1) = Ca(j0:j1,i1+1) .* Ey(j0:j1,i1+1) - Cb(j0:j1,i1+1) .* Hzi(j0:j1,i1)   / ds;   % i=i1+1/2 Face
    
    
    
    %--------------------Gaussian beam
    Hzi(y_position-y_wing:y_position+y_wing,x_position) = H_source.*sin(omega*t*dt);
    %--------------------Gaussian beam
    
    
    
    if mod(t,5)==0
        
        %
        surf(Hz(j0:j1,i0:i1));
        shading flat;
        colorbar;
        caxis([-3 3])
        axis([i0 i1 j0 j1 -3 3]);
        xlabel('x axis (10 nm)'); 
        ylabel('y axis (10 nm)'); 
        zlabel('Amplitude (A/m)');
        title(['Hz t = ', num2str(t)]);
        %}
        %{
        pcolor(Hz(j0:j1,i0:i1));
        shading flat;
        colorbar;
        caxis([-3,3]);
        xlabel('x axis (10 nm)'); 
        ylabel('y axis (10 nm)'); 
        title(['t = ', num2str(t)]);
        %}
        
        drawnow;
        
        %frame(t/5) = getframe(gcf); % get the frame
        
    end
    
end


%writegif('Hz_circle_array.gif',frame,0.05);


