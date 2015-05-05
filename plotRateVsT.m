% Title : Plot rate versus T (Fig. 5) 
% File  : plotRateVsT.m
% -------------------------------------------------------------------------
% Description :
%   Plots acheivable rates for QPSK and 16QAM as a function of the coherence
%   interval. The number of pilots are optimized for all T. Also computed
%   are the corresponding perfect CSIR rates.
% ------------------------------------------------------------------------- 
% Revisions   :
%   Date       Version  Author  Description
%   26-jan-15  1.0      svenja  created the file
% -------------------------------------------------------------------------                            
% Copyright (c) 2015, Sven Jacobsson
% All rights reserved.
% ------------------------------------------------------------------------- 

disp('plotting acheivable rate vs T...');

% PARAMETERS
%-------------------------------------------------------------------------%

snrdb = -10;            % SNR [dB]
N = 400;                % number of antennas
T = logspace(0,5,50);   % coherence interval [channel uses]
b = 1;                  % ADC resolution [bits]
K = 20;                 % number of users
rec = 'MRC';            % linear receiver (MRC/ZF)

% SIMULATION
%-------------------------------------------------------------------------%

% 16QAM rate with perfect CSIR
R16QAM_CSIR = simRate(snrdb, N, K, T, 0, b, 16, 'QAM', rec);

% 16QAM rate with no a priori CSI
R16QAM  = pilotOptSoft(@(Pk) simRate(snrdb, N, K, T, Pk, b, 16, 'QAM', rec), T, K);

% QPSK rate with perfect CSIR 
RQPSK_CSIR = simRate(snrdb, N, K, T, 0, b, 4, 'PSK', rec);

% QPSK rate with no a priori CSI
RQPSK  = pilotOptSoft(@(Pk) simRate(snrdb, N, K, T, Pk, b, 4, 'PSK', rec), T, K);

% PLOT RESULTS
%-------------------------------------------------------------------------%

fig_RateVsT = figure(5); 
clf;
set(fig_RateVsT,'name','Rate vs T');
legh = {};
semilogx(T,R16QAM_CSIR,'k--x'); legh = [legh, '16QAM (CSIR)'];
hold on;
semilogx(T,R16QAM,'r-x'); legh = [legh, '16QAM'];
semilogx(T,RQPSK_CSIR,'k--o'); legh = [legh, 'QPSK (CSIR)'];
semilogx(T,RQPSK,'b-o'); legh = [legh, 'QPSK'];
grid on;
legh = legend(legh);
set(legh,'location','southeast','fontsize',14);
xlabel('T','fontsize',20);
ylabel('Rate per user [bits/channel use]','fontsize', 20); ylim([0 4]);
