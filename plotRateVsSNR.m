% Title : Plot rate versus SNR (Fig. 4) 
% File  : plotRateVsSNR.m
% -------------------------------------------------------------------------
% Description :
%   Plots acheivable rates for QPSK and 16QAM as a function of the SNR. The 
%   rate is computed both in the single-user (K=1) and multi-user (K=20) 
%   case. The number of pilots are optimized for each SNR.
% ------------------------------------------------------------------------- 
% Revisions   :
%   Date       Version  Author  Description
%   26-jan-15  1.0      svenja  created the file
% -------------------------------------------------------------------------                            
% Copyright (c) 2015, Sven Jacobsson
% All rights reserved.
% ------------------------------------------------------------------------- 

disp('plotting acheivable rate vs SNR...');

% PARAMETERS
%-------------------------------------------------------------------------%

snrdb = -30:2:10;       % SNR [dB]
N = 400;                % number of antennas
T = 1000;               % coherence interval [channel uses]
b = 1;                  % ADC resolution [bits]
K = 20;                 % number of users
rec = 'MRC';            % linear receiver (MRC/ZF)

% initialize variables

R16QAM_SU = nan(size(snrdb));
R16QAM_MU = nan(size(snrdb));

RQPSK_SU = nan(size(snrdb));
RQPSK_MU = nan(size(snrdb));

% SIMULATION
%-------------------------------------------------------------------------%

for i = 1:length(snrdb)
    
    fprintf('\nSNR = %i dB\n\n', snrdb(i));
    
    % 16QAM single-user rate
    R16QAM_SU(i)  = pilotOptSoft(@(Pk) simRate(snrdb(i), N, 1, T, Pk, b, 16, 'QAM', rec), T, 1);
    
    % 16QAM multi-user rate
    R16QAM_MU(i)  = pilotOptSoft(@(Pk) simRate(snrdb(i), N, K, T, Pk, b, 16, 'QAM', rec), T, K);
    
    % QPSK single-user rate
    RQPSK_SU(i)  = pilotOptSoft(@(Pk) simRate(snrdb(i), N, 1, T, Pk, b, 4, 'PSK', rec), T, 1);
    
    % QPSK multi-user rate
    RQPSK_MU(i)  = pilotOptSoft(@(Pk) simRate(snrdb(i), N, K, T, Pk, b, 4, 'PSK', rec), T, K);
 
end

% PLOT RESULTS
%-------------------------------------------------------------------------%

fig_RateVsN = figure(6); 
clf;
set(fig_RateVsN,'name','Rate vs N');
legh = {};
hold on;
plot(snrdb,R16QAM_SU,'r--x'); legh = [legh, '16QAM, 1 user'];
plot(snrdb,R16QAM_MU,'r-x'); legh = [legh, '16QAM, 20 users'];
plot(snrdb,RQPSK_SU,'b--o'); legh = [legh, 'QPSK, 1 user'];
plot(snrdb,RQPSK_MU,'b-o'); legh = [legh, 'QPSK, 20 users'];
grid on;
legh = legend(legh);
set(legh,'location','southeast','fontsize',14);
xlabel('SNR [dB]','fontsize',20);
ylabel('Rate per user [bits/channel use]','fontsize', 20); ylim([0 4]);


