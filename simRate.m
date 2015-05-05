function [RSOFT, RHARD] = simRate(snrdb, N, K, T, Pk, b, M, sig, rec, plotout)
% Title : Acheivable rate computation: one-bit massive MIMO uplink
% File  : simRate.m
% -------------------------------------------------------------------------
% Description :
%   Computes the acheivable rate for PSK/QAM signalling in a one-bit
%   massive MIMO uplink system. The receiver employs MRC or ZF to map the
%   received signal over all antennas to a scalar output. When the channel 
%   fading coefficients that are required for MRC/ZF are unknown to the
%   receiver, LS channel estimation is used to estimate the channel based 
%   on pilot symbols.
%   The probabilities that are required to compute the acheivable rates are 
%   approximated by Monte Carlo analysis. The MRC/ZF output domain is 
%   discretized to a rectangular grid. For details see:
%   [S. Jacobsson, "Throughput analysis of massive MIMO uplink with one-bit
%   ADCs", Master's thesis, Chalmers University of Technology, Gothenburg, 
%   Sweden, Jan. 2015].
%   If the function is called without any inputs it will run the simulation
%   for some default parameters and display the receiver output.
%   
%   Output parameters:
%       RSOFT = Rate per user [b/s/Hz] with soft decoding
%       RHARD = Rate per user [b/s/Hz] with hard decoding
%
%   Input parameters:
%       snrdb   =   SNR in dB
%       N       =   number of antennas
%       K       =   number of users
%       T       =   cohrence interval measured in channel uses
%       Pk      =   number of pilots per user (0 = perfect CSIR)
%       b       =   number of quantization bits (0 = inf.prec., 1 = one-bit)
%       M       =   constellation size
%       sig     =   signalling scheme (PSK or QAM)
%       rec     =   receiver architechture (MRC or ZF)
%       plotout =   plot receiver output
%
%   Simulation parameters
%       numC    =   number of channel realizations (default = 100).
%       numU    =   number of noise/interference realizations per channel realization(default = 100).
%       numbin  =   number of output bins in discrete grid to compute soft probabilitie   
%
%   Display parameters:
%       showAxes    :   display x/y-axis (defauilt = 1)
%       showGrid    :   display MRC/ZF output discretization grid (default = 0)
%       showInpt    :   display input alphabet (default = 0)
%
% ------------------------------------------------------------------------- 
% Revisions   :
%   Date       Version  Author  Description
%   26-jan-15  1.0      svenja  created the file
%   29-jan-15  1.1      svenja  added description and comments
% -------------------------------------------------------------------------                              
% Copyright (c) 2015, Sven Jacobsson
% All rights reserved.
% =========================================================================


% INITIALIZATION PHASE
%-------------------------------------------------------------------------%

% include mtimesx directory for fast computation of 3-dim matrices.
addpath mtimesx;

% number of channel realisations and noise/interference realisations
numC = 200;
numU = 200;

% number of output bins per dimension for the discrete grid on which the
% soft probabilities are approximated
numbin = 100;

% display axes/grid/inputs
showAxes = 1;
showGrid = 0;
showInpt = 0;

% total number of inputs/outputs for MC simulation.
MAX = numC*numU;

tic;
fprintf('Running simRate.m for %i channel realisations and %i noise/interference realisations...\n ',numC,numU);

% set default parameters
if nargin == 0
    
    snrdb = 0;
    N = 400;
    K = 1;
    T = 1:1000;
    Pk = 20;
    b = 1;
    M = 16;
    sig = 'QAM';
    rec = 'MRC';
    
    % display receiver output signal
    plotout = 1; 
    
    % figure handle
    fig_out = figure(99); clf; set(fig_out,'name','Receiver output');
    
elseif nargin == 9
    
    plotout = 0; % do not display receiver output signal
    
end

% SNR in dB to linear units
snr = db2pow(snrdb);

% determine input modulation and cardinality (C = input alphabet)
if strcmp(sig,'PSK')
    C = pskmod(0:M-1,M,pi/M)';
elseif strcmp(sig,'QAM')
    C = qammod(0:M-1,M)';
    QAMnorm = sqrt(mean(abs(C).^2));
    C =C/QAMnorm;
else
    error('Only PSK or QAM is supported.');
end

% Rayleigh fading channel
H = (randn(N,K,numC) + 1i*randn(N,K,numC))/sqrt(2);

% TRAINING & CHANNEL ESTIMATION PHASE
%-------------------------------------------------------------------------%

% total number of pilots
Pt = K*Pk; 

% channel estimatiion: if P = 0 then assume perfect CSIR, otherwise perform
% LS channel estimation based on transmitted pilots.
if Pt == 0
    
    He = H;
    
else
    
    % pilots: each user sends at different time to ensure orthogonality. The
    % pilots are choosen at random from the input alphabet.
    Xp = zeros(K,Pt,numC);
    for k = 1:K
       Xp(k,1+(k-1)*Pk:k*Pk,:) =  ones(1,Pk,numC);
    end
    Xp = sqrt(snr) * sqrt(K) * C(randi([1,M],K,Pt,numC)) .* Xp;

    % additive noise acting on pilots
    Wp = (randn(N,Pt,numC) + 1i*randn(N,Pt,numC))/sqrt(2);

    % received pilot sequence (unquantised)
    Yp = mtimesx(H,Xp) + Wp;

    % quantised pilot sequence
    Rp = quantise(Yp,b);

    % LS estimation
    RXp = mtimesx(Rp,Xp,'c');
    XXp = mtimesx(Xp,Xp,'c');
    He = mtimesx(RXp,multinv(XXp));
    
end

% DATA TRANSMISSION & DETECTION PHASE
%-------------------------------------------------------------------------%

% soft estimate for the first user after MRC/ZF, used to approximate the
% necessary probabilities required to compute acheivable rate. 
xsoft = nan(M,MAX);

for m = 1:1:M
    
    % data symbols: all interfering users have their inputs choosen
    % randomly from the input alphabet.
    Xd = C(randi([1,M],K,numU,numC));
    Xd(1,:,:) = C(m);
    Xd = sqrt(snr)*Xd;

    % additive noise acting on data symbols
    Wd = (randn(N,numU,numC) + 1i*randn(N,numU,numC))/sqrt(2);
    
    % received data sequence (unquantised)
    Yd = mtimesx(H,Xd) + Wd;
    
    % quantised data sequence
    Rd = quantise(Yd,b);
    
    % soft decoding with MRC/ZF filter. The acquired channel estimate is
    % used to map the received signal vector to a scalar channel. 
    if strcmp(rec,'MRC')
        he = He(:,1,:);
        xsoft(m,:) = reshape(mtimesx(1./sum(conj(he).*he),mtimesx(he,'c',Rd)),1,[]);
    elseif strcmp(rec,'ZF')
        HeHe = mtimesx(He,'c',He);
        Xsoft = mtimesx(mtimesx(multinv(HeHe),He,'c'),Rd);
        xsoft(m,:) = reshape(Xsoft(1,:,:),1,[]);
    end
    
end

% hard decoding
if strcmp(sig,'PSK')
    xhard = pskdemod(xsoft+eps,M,pi/M);
elseif strcmp(sig,'QAM')
    xhard = qamdemod(QAMnorm/sqrt(snr)*xsoft+eps,M);
end

% PROBABILITIES & RATE COMPUTATION
%-------------------------------------------------------------------------%

% soft channel law and soft output probability. "edges" is used when 
% plotting the discretization grid.
[pSOFT, poutSOFT, edges] = softProb(xsoft,numbin);

% hard channel law and hard output probability.
[pHARD, poutHARD] = hardProb(xhard);

% acheivable rate with soft decoding.
RSOFT = rate(pSOFT,poutSOFT,Pt,T);

% acheivable rate with hard decoding.
RHARD = rate(pHARD,poutHARD,Pt,T);

%-------------------------------------------------------------------------%
% PLOT RECEIVER OUTPUT
%-------------------------------------------------------------------------%

% plot 1000 first receiver outputs for each input.
if plotout
    
    figure(99); hold all;

    % plot axes
    if showAxes
        line(sqrt(snr)*[-10, 10],[0, 0],'color','k','linewidth',2);
        line([0, 0],sqrt(snr)*[-10, 10],'color','k','linewidth',2);
    end

    % plot grid
    if showGrid
        for i = 1:length(edges)
           line([edges(i) edges(i)], [edges(1), edges(end)],'color', [0.7 0.7 0.7],'linewidth',0.1);
           line([edges(1), edges(end)],[edges(i) edges(i)],'color', [0.7 0.7 0.7],'linewidth',0.1);
        end
    end
    
    % plot MRC output
    if MAX <= 1000
        for m = 1:M
            plot(xsoft(m,:),'rx');
        end
    else
        for m = 1:M
            plot(xsoft(m,1:1000),'rx');
        end       
    end

    % plot input constellation
    if showInpt
        plot(sqrt(snr)*C,'ko','markersize',7);
    end

    axis equal; box on;
    axis(max(edges)*[-1, 1, -1, 1]);
    xlabel('In-phase','fontsize',20);
    ylabel('Quadrature','fontsize',20);
    title([num2str(M),sig,', ',rec,', SNR = ',num2str(snrdb),' dB, K = ',num2str(K),', Pk = ',num2str(Pt),', b = ',num2str(b)],'fontsize',16);
    
end

elapsedTime = toc;
fprintf('\tSimulation finished after %f seconds.\n',elapsedTime);

end

% FUNCTIONS
%-------------------------------------------------------------------------%

function R = quantise(Y,b)
% perform 1-bit quantisation if b = 1, otherwise return unquantised signal.
    if b == 0
        R = Y;
    elseif b == 1
        R = sign(real(Y)) + 1i*sign(imag(Y));
    end

end

function [pSOFT, poutSOFT,edges] = softProb(xsoft, numbin)
% approximate soft channel law and output probability by mapping MRC/ZF
% output on a discrete grid and counting the number of instances that the 
% MRC/ZF outputs occur in each grid region, for each input.

    % extract number of inputs and number of random
    % noise/interference/channel realisations
    [M,MAX] = size(xsoft);

    % find max/min-points of soft outputs to determine the the outer 
    % boundaries of the discretization grid.
    min_edges = min(min(min(real(xsoft))), min(min(imag(xsoft))));
    max_edges = max(max(max(real(xsoft))), max(max(imag(xsoft)))); 
    bnd_edges = 1.1*max(abs(min_edges),abs(max_edges));
    
    % edges that specify the grid boundaries per dimension  
    edges = linspace(-bnd_edges,bnd_edges+eps,numbin);

    % channel law: count number of occurences in each bin
    pSOFT = nan(M,length(edges)-1,length(edges)-1);
    for i = 1:M
        pSOFT(i,:,:) = histcn([real(xsoft(i,:))' ,imag(xsoft(i,:))'],edges,edges);
    end
    pSOFT = reshape(pSOFT,M,(length(edges)-1)^2)/MAX;

    % output probability: average the channel law over all inputs
    poutSOFT = mean(pSOFT);

end

function [pHARD, poutHARD] = hardProb(xhard)
% approximate hard MRC/ZF channel law and output probability by performing
% hard QAM/PSK decoiding and counting the number of instances that each
% hard output symbol occur, for each input symbol.

    % extract number of inputs and number of random
    % noise/interference/channel realisations
    [M,MAX] = size(xhard);

    % channel law
    pHARD = nan(M,M);
    for i = 1:1:M   
        for j = 1:1:M
            pHARD(i,j) = sum(xhard(i,:) == j-1);
        end
    end
    pHARD = pHARD/MAX;

    % output probability
    poutHARD = mean(pHARD);

end

function R = rate(p,pout,Pt,T)
% compute the acheivable rate by computing the mutual information and
% normaling with the number of transmitted pilots.

    % extract number of inputs
    M = size(p,1);

    % mutual information
    R = mean(sum(p .* log2(p ./ (repmat(pout,M,1)+eps) + eps),2));

    % normalisation
    R = (T-Pt)./T * R;
    R = max(R,0);

end

