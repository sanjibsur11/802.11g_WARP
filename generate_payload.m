% To generate payload suitable for TVWS OFDM communications
% 
% Author: Sanjib Sur
% Institute: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 05/22/2014
% 
% Comments: Payload generation for OFDM communications
% 
% 

function [Payload] = generate_payload()


%% Read global variables
tvwsGlobalVars;


%% Data generation

Payload = [];
% To visualize the constellation of payload
payload_constellation = [];

if Txparams.numAntenna == 1
    
    for ofdm_symb = 1:Txparams.numOFDMsymb
        
        data_symb = inbits((ofdm_symb - 1)*Timeparams.Nsd*symb2BitNum + 1 : ...
            (ofdm_symb - 1)*Timeparams.Nsd*symb2BitNum + Timeparams.Nsd*symb2BitNum);
        
        data_symb = reshape(data_symb, symb2BitNum, Timeparams.Nsd).';
        data_symb = fliplr(data_symb);
        data_symb = bi2de(data_symb);
        
        QAM_constellation = modulate(hmodem, data_symb);
        
        payload_constellation = [payload_constellation QAM_constellation'];
        
        % Distribute the data according the data subcarrier mapping
        freqDataBlock = zeros(1, Timeparams.Nst);
        for subc = 1:Timeparams.Nst
            if data_to_subcarrier_map(subc) ~= 0
                freqDataBlock(subc) = QAM_constellation(data_to_subcarrier_map(subc));
            end
        end
        % Add the pilot symbols
        for pilot_pos_idx = 1:length(pilot_position)
            freqDataBlock(pilot_position(pilot_pos_idx)) = 1;
        end
        
        % MATLAB IFFT assumes DC at index 0 and zero frequencies at the
        % middle
        freqDataBlock = fftshift(freqDataBlock);
        % IFFT to generate time domain data signal
        ifftData = ifft(freqDataBlock, Timeparams.Nst);
        % Add cyclic prefix: 25%
        Payload = [Payload ifftData(Timeparams.Nst - (Timeparams.Nst/4) + 1 : Timeparams.Nst) ifftData];
    end
    
    % Scale to span -1,1 range of DAC
    payloadScale = max([max(abs(real(Payload))), max(abs(imag(Payload)))]);
    Payload = (1/payloadScale) * Payload;
    
elseif Txparams.numAntenna == 2
    
    
end


if DEBUG_ON
    figure(901);
	plot(0, 0, 'ob', 'LineWidth', 2, 'MarkerSize', 10);
	hold on
	plot(real(payload_constellation), imag(payload_constellation),'+r', ...
        'LineWidth', 2, 'MarkerSize', 10);  
    hold off
    title('Tx constellation');
end

end