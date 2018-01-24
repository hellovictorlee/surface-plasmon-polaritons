clc;
clear all;
close all;

nx=201;
dx=1.66e-9;
ep=8.854187817*1e-12;
mu=4*pi*1e-7;
c=1/sqrt(ep*mu);
S=1;   % Courant stability
dt=S*dx/c;
k=sqrt(mu/ep);
iL=15;
iR=nx-15;

Eyi=zeros(nx,1);
Hzi=zeros(nx-1,1);
Ey=zeros(nx,1);
Hz=zeros(nx-1,1);
eparray=ep*ones(nx,1);
muarray=mu*ones(nx-1,1);
conductivityE=zeros(nx,1);
conductivityH=zeros(nx-1,1);

Cai=(1 - (dt*conductivityE./(2*eparray))) ./(1 + dt*conductivityE./(2*eparray));
Cbi=(dt ./ (eparray))                     ./(1 + dt*conductivityE./(2*eparray));
Dai=(1 - (dt*conductivityH./(2*muarray))) ./(1 + dt*conductivityH./(2*muarray));
Dbi=(dt ./ (muarray))                     ./(1 + dt*conductivityH./(2*muarray));

%--------------------------------------------------------CPML
d=10;
psi_EyiL=zeros(d,1);
psi_HziL=zeros(d,1);
psi_EyiR=zeros(d,1);
psi_HziR=zeros(d,1);

psi_EyL=zeros(d,1);
psi_HzL=zeros(d,1);
psi_EyR=zeros(d,1);
psi_HzR=zeros(d,1);

m=3.5;
cm=(0.8*(m+1))/(k*dx);  % conductivity maximum value
%cm=0.7*(m+1)/(150*pi*dx);
mL_H_=(d:-1:1)';
mL_E_=(d-0.5:-1:0.5)';
segmaH_L=cm*((mL_H_./d).^m);
segmaE_L=cm*((mL_E_./d).^m);


km=1.1;   % kapa maximun
kapaH_L=1+(km-1)*((mL_H_./d).^m);
kapaE_L=1+(km-1)*((mL_E_./d).^m);

kapaH_R=wrev(kapaH_L);
kapaE_R=wrev(kapaE_L);


am=0.1;   % a maximun
ma=1;
aH_L=am*((1-(mL_H_./d)).^ma);
aE_L=am*((1-(mL_E_./d)).^ma);


bh_L=exp(-((segmaH_L./kapaH_L)+aH_L)*(dt/ep));
bh_R=wrev(bh_L);

be_L=exp(-((segmaE_L./kapaE_L)+aE_L)*(dt/ep));
be_R=wrev(be_L);


ch_L=((segmaH_L./kapaH_L)./(segmaH_L+(aH_L.*kapaH_L))).*(bh_L-1);
ch_R=wrev(ch_L);

ce_L=((segmaE_L./kapaE_L)./(segmaE_L+(aE_L.*kapaE_L))).*(be_L-1);
ce_R=wrev(ce_L);
%--------------------------------------------------------CPML


%-----------------------Drude model
freq = 0:1e12:4500e12;
omega = 2*pi*freq;

charge_e = 1.60217646*1e-19;         % 1ev = 1.60217646*1e-19 Joul
mass_e  = 9.10938188*1e-31;          % mass of electron

r_p=zeros(nx,1);
J_p=zeros(nx,1);
omega_p=zeros(nx,1);


r_p_value=2*pi*6.46e12;
omega_p_value=2*pi*2.183e15;
%omega_p_value=charge_e*sqrt(5.9e28/((0.642*mass_e)*ep));

r_p(nx/2:nx,1)=r_p_value;
omega_p(nx/2:nx,1)=omega_p_value;
ep_=11.575;

eparray(nx/2:nx,1)=ep_*ep;

k_p    = (1-r_p*dt/2)       ./ (1+r_p*dt/2);
beta_p = (omega_p.^2*ep*dt/2) ./ (1+r_p*dt/2);



%-----------------------Drude model
%-----------------------Drude model test
epr  = 1;   % incident

epr1 = ep_ - omega_p_value^2             ./         (omega.^2 + r_p_value^2) ;   % Drude model(real part)    epr1 = ep_ - omega_p^2 ./ (omega.^2 + r_p^2)
epr2 =       omega_p_value^2 * r_p_value ./ (omega.*(omega.^2 + r_p_value^2));   % Drude model(image part)   epr2 = omega_p^2 * r_p ./ (omega.*(omega.^2 + r_p^2))

n = real(sqrt(epr1 + 1i*epr2));
k = imag(sqrt(epr1 + 1i*epr2));
reflectance_true = abs( ((n+1i*k) - sqrt(epr)) ./ ((n+1i*k) + sqrt(epr)) ).^2;

nfreq = length(freq);   
Eyiphasor = zeros(1,nfreq);
Eyrphasor = zeros(1,nfreq);
%-----------------------Drude model test


Ca=(2*eparray - dt*beta_p - conductivityE*dt) ./ (2*eparray + dt*beta_p + conductivityE*dt);
Cb=(2*dt)                                     ./ (2*eparray + dt*beta_p + conductivityE*dt);
Da=(1-(dt*conductivityH./(2*muarray)))        ./ (1+ dt*conductivityH./(2*muarray));
Db=(dt./(muarray))                            ./ (1+ dt*conductivityH./(2*muarray));



E=1;   % amplitude of Gaussian pulse
f_wing=2250e12;
sigma=3/(2*pi*f_wing);
m=4*sigma;
f_center=2250e12;



for t=1:10000
    
    %-----------------------incident
    psi_HziL = bh_L .* psi_HziL + ch_L .* (Eyi(2:d+1)   - Eyi(1:d))     / dx;   % Left
    psi_HziR = bh_R .* psi_HziR + ch_R .* (Eyi(nx-d+1:nx) - Eyi(nx-d:nx-1)) / dx;   % Right
    
    Hzi(1:d)       = Dai(1:d)       .* Hzi(1:d)       - Dbi(1:d)       .* (Eyi(2:d+1)   - Eyi(1:d))      ./ (kapaH_L * dx);   % Left
    Hzi(d+1:nx-d-1) = Dai(d+1:nx-d-1) .* Hzi(d+1:nx-d-1) - Dbi(d+1:nx-d-1) .* (Eyi(d+2:nx-d) - Eyi(d+1:nx-d-1)) /dx;   % Middle
    Hzi(nx-d:nx-1)   = Dai(nx-d:nx-1)   .* Hzi(nx-d:nx-1)   - Dbi(nx-d:nx-1)   .* (Eyi(nx-d+1:nx) - Eyi(nx-d:nx-1))  ./ (kapaH_R * dx);   % Right
    
    Hzi(1:d)     = Hzi(1:d)     - Dbi(1:d)     .* psi_HziL;   % Left
    Hzi(nx-d:nx-1) = Hzi(nx-d:nx-1) - Dbi(nx-d:nx-1) .* psi_HziR;   % Right
    
    
    
    psi_EyiL = be_L .* psi_EyiL + ce_L .* (Hzi(2:d+1) - Hzi(1:d))         / dx;   % Left
    psi_EyiR = be_R .* psi_EyiR + ce_R .* (Hzi(nx-d:nx-1) - Hzi(nx-d-1:nx-2)) / dx;   % Right
    
    Eyi(2:d+1)     = Cai(2:d+1)     .* Eyi(2:d+1)     - Cbi(2:d+1)     .* (Hzi(2:d+1)     - Hzi(1:d))       ./ (kapaE_L * dx);   % Left
    Eyi(d+2:nx-d-1) = Cai(d+2:nx-d-1) .* Eyi(d+2:nx-d-1) - Cbi(d+2:nx-d-1) .* (Hzi(d+2:nx-d-1) - Hzi(d+1:nx-d-2))  /dx;   % Middle
    Eyi(nx-d:nx-1)   = Cai(nx-d:nx-1)   .* Eyi(nx-d:nx-1)   - Cbi(nx-d:nx-1)   .* (Hzi(nx-d:nx-1)   - Hzi(nx-d-1:nx-2)) ./ (kapaE_R * dx);   % Right
    
    Eyi(2:d+1)   = Eyi(2:d+1)   - Cbi(2:d+1)   .* psi_EyiL;   % Left
    Eyi(nx-d:nx-1) = Eyi(nx-d:nx-1) - Cbi(nx-d:nx-1) .* psi_EyiR;   % Right
    %-----------------------incident
    
    
    
    
    
    
    %-----------------------TF/SF
    psi_HzL = bh_L .* psi_HzL + ch_L .* (Ey(2:d+1)   - Ey(1:d))     / dx;   % Left
    psi_HzR = bh_R .* psi_HzR + ch_R .* (Ey(nx-d+1:nx) - Ey(nx-d:nx-1)) / dx;   % Right
    
    Hz(1:d)       = Da(1:d)       .* Hz(1:d)       - Db(1:d)       .* (Ey(2:d+1)   - Ey(1:d))      ./ (kapaH_L * dx);   % Left
    Hz(d+1:nx-d-1) = Da(d+1:nx-d-1) .* Hz(d+1:nx-d-1) - Db(d+1:nx-d-1) .* (Ey(d+2:nx-d) - Ey(d+1:nx-d-1)) / dx;   % Middle
    Hz(nx-d:nx-1)   = Da(nx-d:nx-1)   .* Hz(nx-d:nx-1)   - Db(nx-d:nx-1)   .* (Ey(nx-d+1:nx) - Ey(nx-d:nx-1))  ./ (kapaH_R * dx);   % Right
    
    Hz(1:d)     = Hz(1:d)     - Db(1:d)     .* psi_HzL;   % Left
    Hz(nx-d:nx-1) = Hz(nx-d:nx-1) - Db(nx-d:nx-1) .* psi_HzR;   % Right
    
    
    
    Ey_temp = Ey;
    
    psi_EyL = be_L .* psi_EyL + ce_L .* (Hz(2:d+1) - Hz(1:d))         / dx;   % Left
    psi_EyR = be_R .* psi_EyR + ce_R .* (Hz(nx-d:nx-1) - Hz(nx-d-1:nx-2)) / dx;   % Right
    
    Ey(2:d+1)     = Ca(2:d+1)     .* Ey(2:d+1)     - Cb(2:d+1)     .* ((Hz(2:d+1)     - Hz(1:d))       ./ (kapaE_L * dx) + (1/2)*(1+k_p(2:d+1))    .*J_p(2:d+1));   % Left
    Ey(d+2:nx-d-1) = Ca(d+2:nx-d-1) .* Ey(d+2:nx-d-1) - Cb(d+2:nx-d-1) .* ((Hz(d+2:nx-d-1) - Hz(d+1:nx-d-2))  / dx             + (1/2)*(1+k_p(d+2:nx-d-1)).*J_p(d+2:nx-d-1));   % Middle
    Ey(nx-d:nx-1)   = Ca(nx-d:nx-1)   .* Ey(nx-d:nx-1)   - Cb(nx-d:nx-1)   .* ((Hz(nx-d:nx-1)   - Hz(nx-d-1:nx-2)) ./ (kapaE_R * dx) + (1/2)*(1+k_p(nx-d:nx-1))  .*J_p(nx-d:nx-1));   % Right
    
    Ey(2:d+1)   = Ey(2:d+1)   - Cb(2:d+1)   .* psi_EyL;   % Left
    Ey(nx-d:nx-1) = Ey(nx-d:nx-1) - Cb(nx-d:nx-1) .* psi_EyR;   % Right
    
    J_p = k_p .* J_p + beta_p .* (Ey+Ey_temp);
    %-----------------------TF/SF
    
    
    
    
    
    
    Hz(iL-1) = Dai(iL-1) * Hz(iL-1) + Dbi(iL-1) * Eyi(iL)   / dx;   % left side of TF/SF interface
    Ey(iL)   = Cai(iL)   * Ey(iL)   + Cbi(iL)   * Hzi(iL-1) / dx;   % left side of TF/SF interface   
    
    
    
    
    Eyi(12) = E*exp(-(t*dt-m)^2/(2*(sigma)^2))*(cos(2*pi*f_center*(t*dt-m)));
    
    
    %-----------------------DFT
    kernelE = exp(-1j*omega*t*dt); 
    
    Eyiphasor = Eyiphasor + Eyi(12)*kernelE; 
    
    Eyrphasor = Eyrphasor + Ey(12)*kernelE;
    %-----------------------DFT
    
    
    
    if mod(t,100) == 0
        
        subplot(2,1,1) 
        plot(0:dx:(nx-1)*dx,Eyi);
        axis([0 (nx-1)*dx -1 1]); 
        xlabel('meter'); 
        ylabel('V/m'); 
        title(['Eyi at timestep n = ', num2str(t)]);      
        
        subplot(2,1,2)   
        plot(0:dx:(nx-1)*dx,Ey);
        axis([0 (nx-1)*dx -1 1]);   
        xlabel('meter'); 
        ylabel('V/m'); 
        title(['Eyt at timestep n = ', num2str(t)]); 
        
        drawnow;
        
    end
        
end


reflectance = (abs(Eyrphasor./Eyiphasor)).^2;


figure(2)
plot(freq,reflectance,'b',freq,reflectance_true,'r')
xlabel('frequency (THz)')
ylabel('Reflectance')
legend('reflectance','reflectance_true')
grid on

figure(3)
plot(c./freq,n,'r')
axis([0.45e-6 0.95e-6 0 0.13])
xlabel('wavelength')
%ylabel('')
legend('theory n','n')
grid on

figure(4)
plot(c./freq,k,'r')
axis([0.45e-6 0.95e-6 0 7])
xlabel('wavelength')
%ylabel('')
legend('theory k','k')
grid on
