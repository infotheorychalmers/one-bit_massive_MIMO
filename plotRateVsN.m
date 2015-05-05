% Title : Plot rate versus the number of antennas (Fig. 6) 
% File  : plotRateVsN.m
% -------------------------------------------------------------------------
% Description :
%   Plots acheivable rates for QPSK and 16QAM as a function of the number of
%   BS antennas (N). The rate is computed both in the single-user (K=1) and
%   multi-user (K=20) case. The number of pilots are optimized for all N.
% ------------------------------------------------------------------------- 
% Revisions   :
%   Date       Version  Author  Description
%   26-jan-15  1.0      svenja  created the file
% -------------------------------------------------------------------------                            
% Copyright (c) 2015, Sven Jacobsson
% All rights reserved. 
% ------------------------------------------------------------------------- 

disp('plotting acheivable rate vs N...');

% PARAMETERS
%-------------------------------------------------------------------------%

snrdb = -10;            % SNR [dB]
N = 50:25:500;          % number of antennas
T = 1000;               % coherence interval [channel uses]
b = 1;                  % ADC resolution [bits]
K = 20;                 % number of users
rec = 'MRC';            % linear receiver (MRC/ZF)

% initialize variables

R16QAM_SU = nan(size(N));
R16QAM_MU = nan(size(N));

RQPSK_SU = nan(size(N));
RQPSK_MU = nan(size(N));

% SIMULATION
%-------------------------------------------------------------------------%

for n = 1:length(N)
    
    fprintf('\nSNR = %i dB\n\n', N(n));
    
    % 16QAM single-user rate
    R16QAM_SU(n)  = pilotOptSoft(@(Pk) simRate(snrdb, N(n), 1, T, Pk, b, 16, 'QAM', rec), T, 1);
    
    % 16QAM multi-user rate
    R16QAM_MU(n)  = pilotOptSoft(@(Pk) simRate(snrdb, N(n), K, T, Pk, b, 16, 'QAM', rec), T, K);
    
    % QPSK single-user rate
    RQPSK_SU(n)  = pilotOptSoft(@(Pk) simRate(snrdb, N(n), 1, T, Pk, b, 4, 'PSK', rec), T, 1);
    
    % QPSK multi-user rate
    RQPSK_MU(n)  = pilotOptSoft(@(Pk) simRate(snrdb, N(n), K, T, Pk, b, 4, 'PSK', rec), T, K);
 
end

% PLOT RESULTS
%-------------------------------------------------------------------------%

fig_RateVsN = figure(6); 
clf;
set(fig_RateVsN,'name','Rate vs N');
legh = {};
hold on;
plot(N,R16QAM_SU,'r--x'); legh = [legh, '16QAM, 1 user'];
plot(N,R16QAM_MU,'r-x'); legh = [legh, '16QAM, 20 users'];
plot(N,RQPSK_SU,'b--o'); legh = [legh, 'QPSK, 1 user'];
plot(N,RQPSK_MU,'b-o'); legh = [legh, 'QPSK, 20 users'];
grid on;
legh = legend(legh);
set(legh,'location','southeast','fontsize',14);
xlabel('Number of antennas, N','fontsize',20);
ylabel('Rate per user [bits/channel use]','fontsize', 20); ylim([0 4]);


