%% Supplementary  Materials for Manuscript:
%% "On spectral changes of the seismic wave energy by a partially saturated crack due to the hysteresis of liquid bridges phenomenon"
%% Submitted to Geophysics
%% Instruction: 1) Two files: Main.m and MyFunctions.m needs to be saved in the same folder
%%              2) Run Main.m file and change input parameters of interest, described below
%%              3) Equations numbering corresponds to numbering in manuscript.
%% NOTE: Possible error messages!  
%% Comment: Attenuation Q-factor is defined only for sinusoidal waveform
%% ###########################################################
close all,
clear all,
warning('off');
MyFunctions;
%% Input parameters
E            = 30e9;                 % Young's modulus, Pa
nu           = 0.3;                  % Poisson's ratio, -
C            = 2*(1+nu)*(1-2*nu)/E;  % see equation #44
Kwe          = 2.25e9;               % Bulk modulus of liquid, Pa
Knw          = 4e6;                  % Bulk modulus of gas, Pa
gamma        = 73/1e3;               % Interface tension, Pa*m
theta_r      = 5/180*pi;             % Receeding contact angle, radians
theta_a      = 15/180*pi;            % Advancing contact angle, radians
n_cr         = 1e-3;                 % crack porosity
a            = 1e-3;                 % major semi-axis, m
b            = a/1e3;                % minor semi-axis, m
pwe          = 0e6;                  % Pressure in wetting phase (water), Pa
Stress       = -10e6;                % Total Confining Stress,  Pa (negative in compression)
dStress_amp  = 5e3;                  % Wave-stress amplitude, Pa
f            = 30;                   % frequency, [Hz]
Swe          = 0.5;                  % Saturation degree of the liquid, fraction
%%  PreCalculus ad Parameters        
theta        = (theta_a+theta_r)/2;  % initial contact angle (theta_r<theta<theta_a)
pcl          = b/a*E/2/(1-nu^2);     % Crack Closure pressure calculation
tmax         = 24/f;                 % Maximum time (number of periods)
time         = (0:1e-4/f:tmax);      % time, sec 
dStress      = dStress_amp*(sin(2*pi*f.*time)); % sinusoidal wave
% f1=30; f2=40; dStress      = dStress_amp*(sin(2*pi*f1.*time)+sin(2*pi*f2.*time))/2; % bi-sinusoidal wave
% dStress      = dStress_amp*((1-2*pi^2*f^2.*(time-0.08).^2).*exp(-pi^2*f^2.*(time-0.08).^2)); % Ricker waveform
if Stress+pwe>0
   error('Effective Confining Stress is negative') 
end
if pcl+Stress+pwe<0
   error('Effective confining stress is too hight and the crack is closed') 
end
if theta_a<=theta_r
     error('Choose theta_r<theta_a') 
end 
if theta>=theta_a || theta<=theta_r
     error('Choose theta_r<theta<theta_a') 
end 
if b>0.1*a
   error('Choose b<<a') 
end
if Swe>1
    error('Units for Swe is fraction, not %')
end
if theta_a      >= pi/2 
    error('Choose  theta_a<pi/2 %')
end
if theta_r      < 0 
    error('Choose  theta_r>=0 %')
end

%% Assign initial conditions
[beta,~]              = PredictBeta(Swe,a,b,pcl,Stress,pwe,gamma,theta); % Calculation of the contact line location for given Swe
[~,pcap,V_cr,Vnw,c]   = PcapSwe(a,b,pcl,Stress,pwe,beta,gamma,theta);    % Calculation of the initial pcap, V_cr and Vnw_cr
SOLUTION              = zeros(6,length(time));                           % assign solution array
SOLUTION(1,1)         = V_cr;                                            % Crack volume (m^2) per unit length perpendicular to the plane-strain
SOLUTION(2,1)         = Vnw;                                             % Gas volume (m^2) per unit length perpendicular to the plane-strain
SOLUTION(3,1)         = pwe;                                             % Pressure in the liquid phase, Pa
SOLUTION(4,1)         = pcap;                                            % Capillary pressure, Pa
SOLUTION(5,1)         = theta;                                           % Contact angles, radians
SOLUTION(6,1)         = beta;                                            % Contact line location (see eq #3)
%% Incremental (Iterative) Solution of equations 12-16:
for k=2:length(time)
    dstress_incr      = dStress(k)-dStress(k-1); % wave-Stress increment assignment
    [IncrSol]         = Pinning(a,b,pcl,Stress+dStress(k-1),dstress_incr,Kwe,Knw,gamma,SOLUTION(:,k-1)); % see equations #17&18 
    SOLUTION(:,k)     = SOLUTION(:,k-1) + IncrSol';
    if SOLUTION(5,k)  <= theta_r || SOLUTION(5,k) >= theta_a % Check the conditions if the contact line is slipping
        [IncrSol]     = Slipping(a,b,pcl,Stress+dStress(k-1),dstress_incr,Kwe,Knw,gamma,SOLUTION(:,k-1)); % see equations #19&20
        SOLUTION(:,k) = SOLUTION(:,k-1) + IncrSol';
    end
end
%% ASSIGN SOLUTION
SIGMA         = Stress+dStress;
PWE           = SOLUTION(3,:);
PCAP          = SOLUTION(4,:);
V_CR          = SOLUTION(1,:);
VNW           = SOLUTION(2,:);
BETA          = SOLUTION(6,:) ;
THETA         = SOLUTION(5,:) ;
%%  POST CALC
Strain_cr     = (V_CR-V_cr)/V_cr; % Crack-Strain, see equation #42 
Strain_REV    = C*dStress+n_cr*Strain_cr; % REV Strain, see equation #43
%% % the compressibility and Q are defined only for sinusoidal waveform
[Q_REV,C_REV] = Q_Comp(Strain_REV,dStress,f,time,dStress_amp) % see equations #45 & 48

%% FIGURES
figure(3) % Figure #3
subplot(311)
plot(Strain_cr(f*time<=5),dStress(f*time<=5),'k','Linewidth',2)
ylabel('\Delta\sigma, Pa','FontSize',14)
xlabel('\Delta\epsilon_c_r','FontSize',14)
grid on
axis tight
title('(a)','FontSize',14)
subplot(312)
plot(Strain_cr(f*time<=5),THETA(f*time<=5)*180/pi,'k','Linewidth',2)
xlabel('\Delta\epsilon_c_r','FontSize',14)
ylabel('\theta','FontSize',16)
grid on
axis tight
title('(b)','FontSize',14)
subplot(313)
plot(Strain_cr(f*time<=5),cos(BETA(f*time<=5)),'k','Linewidth',2)
xlabel('\Delta\epsilon_c_r','FontSize',14)
ylabel('c/a','FontSize',14),
axis tight
grid on
title('(c)','FontSize',14)

figure(4)% Figure #4
subplot(311)
plot(time,dStress,'k','Linewidth',2)
hold on
plot(time,PWE-PWE(1),'--k','Linewidth',2) 
plot(time,PCAP-PCAP(1),':k','Linewidth',2)
axis([0 10/f -dStress_amp*1.1 dStress_amp*1.1])
legend('\Delta\sigma','\Deltap_w_e','\Deltap_c_a_p')
xlabel('Time, sec','FontSize',14)
ylabel('Wave Stress, Pa','FontSize',14)
grid on
title('(a)','FontSize',14)
subplot(312)
plot(time,Strain_REV,'k','Linewidth',2)
axis([0 10/f -1.01*max(abs(Strain_REV)) 1.01*max(abs(Strain_REV))])
grid on
legend('\Delta\epsilon_R_E_V')
xlabel('Time, sec','FontSize',14)
ylabel('Wave Strain','FontSize',14)
title('(b)','FontSize',14)
subplot(313)
[~,Spectr_Strain_rate]=Spectrum(time(2:end),diff(Strain_REV)); % Spectral amplitudes of strain-rate
[Freq,Spectr_Stress_rate]=Spectrum(time(2:end),diff(dStress)); % Spectral amplitudes of stress-rate
semilogy(Freq,(Spectr_Stress_rate),'r','Linewidth',4)
hold on
semilogy(Freq,(Spectr_Strain_rate),'k','Linewidth',2)
grid on, hold on
axis([0 200 1e-3 1])
xlabel('Frequency, Hz','FontSize',14)
ylabel('Spectrum','FontSize',14)
legend('Spectrum of wave-stress-rate', 'Spectrum of wave-strain-rate')
title('(c)','FontSize',14)

