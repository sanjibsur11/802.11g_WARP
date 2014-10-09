% Channel estimation for TVWS communication
% 
% Author: Sanjib Sur
% Institute: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 05/22/2014
% 
% Comments: Channel estimation for both SISO and MIMO OFDM communication
% TVWS. This function outputs the channel vector in frequency domain.
% 

function [channel_vector dataPos] = channel_estimation_rx(LTFpeakPos, ...
                                    sigin)
                                

%% Read global variables 
tvwsGlobalVars;


%% Channel estimation starts here

st = LTFpeakPos - Timeparams.Nst/2;

if VERBOSE2
	fprintf('LTFpeakPos = %d\n', LTFpeakPos);
end

%Estimate freq offset    
fftoffset = 0;
op = Timeparams.Nst/2 - fftoffset;
% Get the angular distortion
co = sum(sigin(st + op + 1 : st + op + Timeparams.Nst).*...
    conj(sigin(st + op + Timeparams.Nst + 1 : st + op + (2*Timeparams.Nst))));
iFreq = angle(co)/Timeparams.Nst;

% Fc is the carrier frequency
CFO_val = (Fc * iFreq)/(2*pi);
if VERBOSE1
    fprintf('CFO between Tx and Rx: %g kHz\n', CFO_val);
end

measured_CFO = [measured_CFO CFO_val];

% TODO: Need to add more than one Rx antenna case
% Estimate channel from each Tx antennas to each Rx antennas
for txantn = 1:Txparams.numAntenna
    for rxantn = 1:Rxparams.numAntenna
        
        % Freq. compensation
        delTheta = iFreq;
        for k = 1:Timeparams.Nst
            e0(k) = sigin(st + op + (txantn - 1) * length(LTF) + k);
            e0(k) = e0(k) * ...
                exp(sqrt(-1) * delTheta * ((txantn-1)* length(LTF) + k-1));
        end
        % estimate channel of each subcarrier
        LTF_FFT = fft(e0, Timeparams.Nst);
        for k = 1 : Timeparams.Nst %zero freq 0+1 does not carry data
            LTF_CH(txantn, k) = LTF_FFT(k)/LTF_shift(k);
        end
        % DC + left +right guard bands does not carry any data
        % DC
        LTF_CH(txantn, 1) = 0;
        % Left guard band
        LTF_CH(txantn, end/2 + 1 : end/2 + 1 + left_guard_len - 1) = zeros(1, left_guard_len);
        % Right guard band
        LTF_CH(txantn, end/2 - right_guard_len + 1 : end/2) = zeros(1, right_guard_len);
        
        % estimate signal power and noise floor from LTF1 and LTF2
        for k = 1 : 2 * Timeparams.Nst
            e1(k) = sigin(st + op + (txantn - 1) * length(LTF) + k);
            e1(k) = e1(k) * ...
                exp(sqrt(-1) * delTheta * ((txantn-1)* length(LTF) + k-1));
        end
        freqLTF1 = fft(e1(1 : Timeparams.Nst));
        freqLTF2 = fft(e1(Timeparams.Nst + 1 : (2 * Timeparams.Nst)));
        
        % Noise power is the variance of the LTF1 and LTF2
        % Noise and Signal power measured from Tx antennas to Rx antenna1
        noisePower(txantn) = sum(abs((freqLTF1-freqLTF2)*(freqLTF1-freqLTF2)'));
        noisePower(txantn) = noisePower/(Timeparams.Nst*(Timeparams.Nsd+Timeparams.Nsp));
        signalPower(txantn) = sum(abs(freqLTF1).^2) - noisePower(txantn);
        LTF_SINR(txantn) = 10*log10(signalPower(txantn)/noisePower(txantn));
        
        if VERBOSE2
            fprintf('RX antenna %d, LTF from Tx ant %d, signal power %g, noise power %g, Estimated SINR %g dB\n', rxantn, txantn, signalPower(txantn), noisePower(txantn), LTF_SINR(txantn));
        end
    end
end

% Find the joint SNR from all the Tx to all the Rx antennas
joint_SINR = 10*log10(sum(signalPower)/sum(noisePower));

channel_vector = LTF_CH;

if DEBUG_ON
    
    for txantn = 1:Txparams.numAntenna
        % Frequency domain channel response
        length_channel_vector = length(channel_vector(txantn,:));
        figure(306+2*(txantn-1));
        subplot(1, 2, 1);
        plot(abs([channel_vector(txantn, end/2 + 1 : end) 0 ...
            channel_vector(txantn, 2:end/2)]), ...
            '-+r', 'LineWidth', 2, 'MarkerSize', 10);
%         plot(abs(channel_vector), ...
%             '-+r', 'LineWidth', 2, 'MarkerSize', 10);        
        ylabel('LTF gain');
        subplot(1, 2, 2);
        plot(unwrap(angle([channel_vector(txantn, end/2 + 1 : end) 0 ...
            channel_vector(txantn, 2:end/2)])), ...
            '-or', 'LineWidth', 2, 'MarkerSize', 10);
%         plot(angle([channel_vector(txantn, end/2 + 1 : end) 0 ...
%             channel_vector(txantn, 2:end/2)]), ...
%             '-or', 'LineWidth', 2, 'MarkerSize', 10);
%         plot(unwrap(angle(channel_vector)), ...
%             '-or', 'LineWidth', 2, 'MarkerSize', 10); 
%         plot(angle(channel_vector), ...
%             '-or', 'LineWidth', 2, 'MarkerSize', 10);           
        ylabel('LTF phase');
        set(gcf,'NextPlot','add');
        axes;
        h = title('Frequency domain channel response');
        set(gca,'Visible','off');
        set(h,'Visible','on');
        
        % Time domain channel response
%         channel_vector_t = ifft(channel_vector(txantn, :));
%         figure(306+2*(txantn-1)+1);
%         plot(10*log10(abs(channel_vector_t).^2));
%         xlabel('Delay');
%         ylabel('Power');
%         title('Time-delay channel response');
    end
    
end


% Calculate the data starting position
dataPos = LTFpeakPos + (2 * Timeparams.Nst) + (Timeparams.Nst/4);

                                
end