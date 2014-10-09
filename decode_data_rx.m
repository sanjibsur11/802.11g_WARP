% Data decoding logic for TVWS
% 
% Author: Sanjib Sur
% Institute: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 05/22/2014
% 
% Comments: This file contains the data decoding logic from OFDM symbols.
% Added pilot compensation which is derived from the original MUMIMO code.
% 
% 

function [outbits] = decode_data_rx(LTFpeakPos, dataPos, channel_vector, sigin)


%% Read global variables 
tvwsGlobalVars;


%% Data decoding starts here

outbits = [];

% To visualize constellation for received data
payload_constellation = [];
numSymbol = 0;
outsymb = [];
while numSymbol < Txparams.numOFDMsymb
    
    % TODO: Need to add data decoding for more than 1 Rx antennas
    for txantn = 1:Txparams.numAntenna
        for rxantn = 1:Rxparams.numAntenna

            dataOFDMsymb = sigin(dataPos + 1 : dataPos + Timeparams.Nst);

            % Freq compensation for data
            for k = 1 : Timeparams.Nst
                dataOFDMsymb(k) = dataOFDMsymb(k) * ...
                exp(sqrt(-1) * iFreq(txantn) * ((dataPos - LTFpeakPos) + k - 1));
            end
            
            % FFT to recover frequency domain data
            fftData = fft(dataOFDMsymb, Timeparams.Nst);         
            
            % Channel compensation
            for subc = 1:Timeparams.Nst
                fftData(subc) = fftData(subc) / channel_vector(txantn, subc);
            end
            
            % Need to add pilot compensation
            if usePilotComp
                if Timeparams.Nst == 64
                    th1 = angle(fftData(64-21+1));
                    th2 = angle(fftData(64-7+1));
                    th3 = angle(fftData(7+1));
                    th4 = angle(fftData(21+1));

                    % phase drift between two adjacent subcarriers         
                    dTheta = ((th3-th1)+(th4-th2))/2/(21+7);

                    % avgTheta is the phase corresponding to the middle subcarrier 0
                    % -21 + (14+28+42)/4 = 0 
                    avgTheta = th1 + (((th2-th1) + (th3-th1) + (th4-th1)) / 4);

                    % phase of leftmost subcarrier
                    th = avgTheta-26*dTheta; 

                    for k = (64-26+1):64
                        fftData(k) =  fftData(k) * exp(-1i*th);             
                        th = th + dTheta; 
                    end 

                    for k = 1:27
                        fftData(k) =  fftData(k) * exp(-1i*th);            
                        th = th + dTheta; 
                    end
                end
            end
            
            % Subcarrier to data mapping assumes guard bands at the left
            % and right end and DC at the middle
            freqDataBlock = fftshift(fftData);
            
            % Redistribute using mapping to extract the QAM constellation
            QAM_constellation = zeros(1, Timeparams.Nsd);
            for subc = 1:Timeparams.Nst
                if data_to_subcarrier_map(subc) ~= 0
                    QAM_constellation(data_to_subcarrier_map(subc)) = freqDataBlock(subc);
                end
            end
            
            payload_constellation = [payload_constellation QAM_constellation];
            data_symb = demodulate(hdemodem, QAM_constellation);
            outsymb = [outsymb data_symb];
                
        end
    end
    
    % Cyclic prefix is 25%
    dataPos = dataPos + Timeparams.Nst + (Timeparams.Nst/4);
    
numSymbol = numSymbol + 1;
end

outbits = de2bi(outsymb);
outbits = fliplr(outbits);
outbits = outbits';
outbits = outbits(:);

if DEBUG_ON
    figure(301);
	plot(0, 0, 'ob', 'LineWidth', 2, 'MarkerSize', 10);
	hold on
	plot(real(payload_constellation), imag(payload_constellation),'+r', ...
        'LineWidth', 2, 'MarkerSize', 10);  
    hold off
    title('Rx constellation');
end

end