% To generate preambles suitable for TV white space MIMO communications.
% 
% Author: Sanjib Sur
% Institute: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 05/21/2014
% 
% Comments: Preamble generation for OFDM PHY layer in TV white space
% 
% 

function [Preamble] = generate_preamble()


%% Read global variables
tvwsGlobalVars;


if Txparams.numAntenna == 1
    Preamble = [STF LTF];
    
elseif Txparams.numAntenna == 2
    zero_vector = zeros(1, length(LTF));
    
    Preamble(1, :) = [STF LTF zero_vector];
    Preamble(2, :) = [STF zero_vector LTF];
    
elseif Txparams.numAntenna == 4
    
end

end
