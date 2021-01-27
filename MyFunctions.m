%% Supplementary  Materials for Manuscript:
%% "On spectral changes of the seismic wave energy by a partially saturated crack due to the hysteresis of liquid bridges phenomenon"
%% Submitted to Geophysics
%% Instruction: 1) Two files: Main.m and MyFunctions.m needs to be saved in the same folder
%%              2) Run Main.m file and change input parameters of interest, described below
%%              3) Corresponding equations numbering corresponds to numbering in manuscript.
%% NOTE: Possible error messages!  
%% Comment: Attenuation Q-factor is defined only for sinusoidal waveform
%% ###########################################################
function MyFunctions
assignin('base','AssignSol',@AssignSol);
assignin('base','Spectrum',@Spectrum);
assignin('base','Pinning',@Pinning);
assignin('base','Slipping',@Slipping);
assignin('base','PcapSwe',@PcapSwe);
assignin('base','PartialDerivatives',@PartialDerivatives);
assignin('base','PredictBeta',@PredictBeta);
assignin('base','Q_Comp',@Q_Comp);
end
%%
function [Q,C]=Q_Comp(Strain,dStress,f,time,dStress_amp) % the compressibility and Q are defined only for sinusoidal waveform
C      = ((max(Strain(f*time>=1))-min(Strain(f*time>1)))/(max(dStress(f*time>1))-min(dStress(f*time>1)))); % see eq.# 45
dW     = trapz(Strain(f*time>=1&f*time<=2),dStress(f*time>=1&f*time<=2));                                  % area of hysteresis loop, see eq.# 47
W      = dStress_amp^2*C/2;                                                                             % see eq.# 46
Q      = 2*pi*W/dW;                                                                                        % see 1/eq.# 48
end
%%
function [Beta,Pcap]=PredictBeta(Swe,a,b,pcl,Stress,pwe,gamma,theta)
% Predict the contact line location eq(3) for a given Swe
beta             = 1e-5:1e-5:pi/2;
[swe,pcap,~,~,~] = PcapSwe(a,b,pcl,Stress,pwe,beta,gamma,theta);
Diff             = abs(Swe-swe);
Ind              = Diff==min(Diff);
Beta             = beta(Ind);
Pcap             = pcap(Ind);
end
%%
function [Freq,SpectrDensity]=Spectrum(time,signal)
% (Fourier) Spectral Amplitude calculation
Fs            = length(time)-1;
X             = fftshift(fft(signal));
X_mag         = abs(X);
df            = -Fs/2:1:Fs/2;
Fr            = df/max(time);
Amp           = X_mag/max(X_mag);
Freq          = Fr;
SpectrDensity = Amp;
end
%%
function [Vcr,Vnw,pwe,pcap,theta,beta]=AssignSol(Solution)
Vcr   = Solution(1);
Vnw   = Solution(2);
pwe   = Solution(3);
pcap  = Solution(4);
theta = Solution(5);
beta  = Solution(6);
end
%%
function [IncrSol]=Pinning(a,b,pcl,Stress,dsigma,Kwe,Knw,gamma,Solution)
%% See equations 17-18
[~,~,pwe,~,theta,beta]=AssignSol(Solution);
% Calculation of private derivatives
[~,~,dpcap_dsigma,dpcap_dpwe,dpcap_dtheta,~,Vcr,dVcr_dpcap,dVcr_dsigma,dVcr_dpwe,~,dVcr_dtheta,Vnw,dVnw_dpcap,dVnw_dsigma,dVnw_dpwe,dVnw_dtheta,~]=PartialDerivatives(a,b,pcl,Stress,pwe,beta,gamma,theta);
% Solution of linear algebra
Matr         = zeros(5,5);
Matr(1,:)    = [Kwe, -Kwe, Vcr - Vnw, 0, 0];
Matr(2,:)    = [0, Knw, Vnw, Vnw, 0];
Matr(3,:)    = [0, 0, dpcap_dpwe, -1, dpcap_dtheta];
Matr(4,:)    = [-1, 0, dVcr_dpwe, dVcr_dpcap, dVcr_dtheta];
Matr(5,:)    = [0, -1, dVnw_dpwe, dVnw_dpcap, dVnw_dtheta];
Vector       = -dsigma*[0, 0, dpcap_dsigma, dVcr_dsigma, dVnw_dsigma]';
Sol          = (Matr)\Vector;
dVtot = Sol(1); dVnw   = Sol(2); dpwe = Sol(3);
dpcap = Sol(4); dtheta = Sol(5); dbeta = 0;
IncrSol      = [dVtot dVnw dpwe dpcap dtheta dbeta];
end
%%
function [IncrSol]=Slipping(a,b,pcl,Stress,dsigma,Kwe,Knw,gamma,Solution)
%% See equations 19-20
% Calculation of private derivatives
[~,~,pwe,~,theta,beta]=AssignSol(Solution);
[~,~,dpcap_dsigma,dpcap_dpwe,~,dpcap_dbeta,Vcr,dVcr_dpcap,dVcr_dsigma,dVcr_dpwe,dVcr_dbeta,~,Vnw,dVnw_dpcap,dVnw_dsigma,dVnw_dpwe,~,dVnw_dbeta]=PartialDerivatives(a,b,pcl,Stress,pwe,beta,gamma,theta);
% Solution of linear algebra
Matr         = zeros(5,5);
Matr(1,:)    = [Kwe, -Kwe, Vcr - Vnw, 0, 0];
Matr(2,:)    = [0, Knw, Vnw, Vnw, 0];
Matr(3,:)    = [0, 0, dpcap_dpwe, -1, dpcap_dbeta];
Matr(4,:)    = [-1, 0, dVcr_dpwe, dVcr_dpcap, dVcr_dbeta];
Matr(5,:)    = [0, -1, dVnw_dpwe, dVnw_dpcap, dVnw_dbeta];
Vector       = -dsigma*[0, 0, dpcap_dsigma, dVcr_dsigma, dVnw_dsigma]';
Sol          = (Matr)\Vector;
dVtot = Sol(1); dVnw  = Sol(2); dpwe   = Sol(3);
dpcap = Sol(4); dbeta = Sol(5); dtheta = 0;
IncrSol      = [dVtot dVnw dpwe dpcap dtheta dbeta];
end
%%
function [Swe,pcap,dpcap_dsigma,dpcap_dpwe,dpcap_dtheta,dpcap_dbeta,Vcr,dVcr_dpcap,dVcr_dsigma,dVcr_dpwe,dVcr_dbeta,dVcr_dtheta,Vnw,dVnw_dpcap,dVnw_dsigma,dVnw_dpwe,dVnw_dtheta,dVnw_dbeta]=PartialDerivatives(a,b,pcl,Stress,pwe,beta,gamma,theta)
% Auxiliary parameters
X            = 1 + 8/pi*pcl*gamma/b .* cos(theta) ./ (pcl + Stress + pwe) .^ 2 ./ sin(beta) .* (beta + cot(beta) .* log(cos(beta))); % eq.21
Y            = (beta + cot(beta) .* log(cos(beta)));                                                                                 % eq.22
pcap         = pi/4 * (pcl + Stress + pwe) .* (sqrt(X) - 1)./Y;                                                                      % eq.5
dpcap_dsigma = pi*(sqrt(X)-1+(1 - X)/sqrt(X))/(4*Y);                                                                                 % eq.23
dpcap_dpwe   = dpcap_dsigma;                                                                                                         % eq.24
dpcap_dtheta = pi*(1-X)*(pcl+Stress+pwe)/8/Y/sqrt(X)*tan(theta);                                                                     % eg.25
dpcap_dbeta  = pi/8*(pcl+Stress+pwe)/sin(beta)/Y/sqrt(X)*((1-X)*cos(beta)+(sqrt(X)-1)^2*log(cos(beta))/sin(beta)/Y);                 % eq.26
% Vcr and derivatives
Vcr          = pi * a * b * (1 + (Stress + pwe + ((2 * beta) - sin((2 * beta))) / pi * pcap) / pcl);                                  % eq.7
dVcr_dsigma  = pi * a * b / pcl;                                                                                                     % eq.27
dVcr_dpwe    = dVcr_dsigma;                                                                                                          % eq.28
dVcr_dtheta  = 0;                                                                                                                    % eq.29
dVcr_dbeta   = a*b*(2 - 2*cos(2*beta))*pcap/pcl;                                                                                     % eq.30
dVcr_dpcap   = a*b*(2*beta - sin(2*beta))/pcl;                                                                                       % eq.31

% Vnw.1 and derivatives
Vnw1         = a * b / pcl .* (pcl + Stress + pwe - pcap .* (4 / pi .* (sin((2 * beta)) .* beta - (beta .^ 2) + 2 .* cos(beta) .^ 2 .* log(cos(beta))) ./ ((2 * beta) - sin((2 * beta))))) .* ((2 * beta) - sin((2 * beta))); % eq.9
dVnw1_dsigma = a*b*(2*beta - sin(2*beta))/pcl;                                                                                       % eq.32
dVnw1_dpwe   = dVnw1_dsigma;                                                                                                         % eq.33
dVnw1_dtheta = 0;                                                                                                                    % eq.34
dVnw1_dbeta  = 2*a*b*(pcl + Stress + pwe + 4*pcap*(beta + sin(2*beta)*log(cos(beta))/(1 - cos(2*beta)))/pi)*(1 - cos(2*beta))/pcl;   % eq.35
dVnw1_dpcap  = 4*a*b*(beta^2 - sin(2*beta)*beta - 2*cos(beta)^2*log(cos(beta)))/(pcl*pi);                                            % eq.36
% Vnw.2 and derivatives
wa           = b/pcl*(pcl+Stress+pwe+(beta+cot(beta)*log(cos(beta)))*pcap*2/pi);                                                     % eq.1
Vnw2         = (pi - 2*theta - sin(2*theta))*wa^2/cos(theta)^2;                                                                      % eq.10
dVnw2_dsigma = 2*b*Vnw2/(pcl*wa);                                                                                                    % eq.37
dVnw2_dpwe   = dVnw2_dsigma;                                                                                                         % eq.38
dVnw2_dtheta = 2*(-2 + tan(theta)*(pi - 2*theta))*wa^2/cos(theta)^2;                                                                 % eq.39
dVnw2_dbeta  = -4*b*pcap*log(cos(beta))*Vnw2/(pcl*pi*sin(beta)^2*wa);                                                                % eq.40
dVnw2_dpcap  = 4*b*(beta + cot(beta)*log(cos(beta)))*Vnw2/(pcl*pi*wa);                                                               % eq.41
% Sum of two terms
Vnw          = Vnw1         + Vnw2;                                                                                                  % eq.8
dVnw_dpcap   = dVnw1_dpcap  + dVnw2_dpcap;
dVnw_dbeta   = dVnw1_dbeta  + dVnw2_dbeta;
dVnw_dsigma  = dVnw1_dsigma + dVnw2_dsigma;
dVnw_dpwe    = dVnw1_dpwe   + dVnw2_dpwe;
dVnw_dtheta  = dVnw1_dtheta + dVnw2_dtheta;
Swe          = (Vcr-Vnw)/Vcr;                                                                                                       % eq.11
end
%%
function [Swe,pcap,Vcr,Vnw,c]=PcapSwe(a,b,pcl,Stress,pwe,beta,gamma,theta)
X            = 1 + 8/pi*pcl*gamma/b .* cos(theta) ./ (pcl + Stress + pwe) .^ 2 ./ sin(beta) .* (beta + cot(beta) .* log(cos(beta)));% eq.21
Y            = (beta + cot(beta) .* log(cos(beta)));                                                                                % eq.22
pcap         = pi/4 * (pcl + Stress + pwe) .* (sqrt(X) - 1)./Y;                                                                     % eq.5
Vcr          = pi * a * b * (1 + (Stress + pwe + ((2 * beta) - sin((2 * beta))) ./ pi .* pcap) ./ pcl);                              % eq.7
Vnw1         = a * b / pcl .* (pcl + Stress + pwe - pcap .* (4 / pi .* (sin((2 * beta)) .* beta - (beta .^ 2) + 2 .* cos(beta) .^ 2 .* log(cos(beta))) ./ ((2 * beta) - sin((2 * beta))))) .* ((2 * beta) - sin((2 * beta)));% eq.9
wa           = b/pcl*(pcl+Stress+pwe+(beta+cot(beta).*log(cos(beta))).*pcap*2/pi);                                                  % eq.1
Vnw2         = (pi - 2*theta - sin(2*theta))*wa.^2/cos(theta)^2;                                                                    % eq.10
Vnw          = Vnw1 + Vnw2;                                                                                                         % eq.8
Vwe          = Vcr-Vnw;
Swe          = Vwe./Vcr;                                                                                                            % eq.11
c            = a*cos(beta);                                                                                                         % eq.3
end
