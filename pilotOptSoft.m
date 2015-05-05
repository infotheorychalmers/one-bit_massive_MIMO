function R = pilotOptSoft(funR,T,K)
% Title : Optimize the number of pilots
% File  : pilotOptSoft.m
% -------------------------------------------------------------------------
% Description :
%   Optimizes the number of pilots transmitted in the one-bit massive MIMO
%   system. The acheivable rate as a function of the number of pilots is
%   concave (see [26] in paper). Therefore, to find the optimal number. 
%   The amount of pilots are incremented until the rate starts decreasing.
%
%   Output parameters:
%       R       =   Rate with soft decoding after pilot optimization
%
%   Input parameters:
%       funR    =   Function handle for the function "simRate.m" with the
%                   parameters "Pk" (pilots per user) as a variable and all
%                   other parameters as constants.
%       T       =   coherence interval.
%       K       =   number of users.
%
% ------------------------------------------------------------------------- 
% Revisions   :
%   Date       Version  Author  Description
%   26-jan-15  1.0      svenja  created the file
% -------------------------------------------------------------------------                            
% Copyright (c) 2015, Sven Jacobsson
% All rights reserved.
% ------------------------------------------------------------------------- 

% initialize variables
R_temp = zeros(length(T),length(T));
Rmax_new = 0;

% increment the number of pilots until the rate no longer increases.
for Pk = 1:1:floor(max(T)/K)
    
    tic; fprintf('P = %i, ',Pk);
    R_temp(Pk,:) = funR(Pk);
    
    if mod(Pk,5) == 0 
        
        Rmax_old = Rmax_new;
        Rmax_new = max(max(R_temp,[],1));
        
        if (Rmax_new == Rmax_old)
            break; 
        end
        
    end    
end

% optimise over pilots length
R = max(R_temp,[],1);

end