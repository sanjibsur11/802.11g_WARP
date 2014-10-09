% TVWS received side code
% 
% Author: Sanjib Sur
% Institute: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 05/21/2014
% Comments: This file contains the receiver side code TV white space.


function [detect_pkt_flag, BER] = tvws_rx(RxData_board)

%% Read global variables
tvwsGlobalVars;


%% Plotting raw samples and psd

if DEBUG_ON
    
    % Raw samples
    length_rx = length(RxData_board);
    figure(200);
    plot(1:length_rx, real(RxData_board), 'b-', ...
         1:length_rx, imag(RxData_board), 'r-');
    xlabel('Sample count');
    ylabel('Rx samples');
    title('Received samples');
    
%     % Receiver FFT
%     sigdft = fft(RxData_board);
%     sigdft = sigdft(1:length_rx/2+1);
%     psdx = (1/(Fs*length_rx)).*abs(sigdft).^2;
%     psdx(2:end-1) = 2*psdx(2:end-1);
%     freq = 0:Fs/length(RxData_board):Fs/2;
%     figure(201);
%     plot(freq, 10*log10(psdx)); 
%     if USESIM
%         if PHASE_NOISE
%             title('PSD with phase noise for received samples');
%         else
%             title('PSD without phase noise for received samples');
%         end
%     else
%         title('PSD for received samples');
%     end       

end

% Downsample the received samples
sigin = decimate(RxData_board, osamp);


if DEBUG_ON
    
    % Raw samples
    length_rx = length(sigin);
    figure(201);
    plot(1:length_rx, real(sigin), 'b-', ...
         1:length_rx, imag(sigin), 'r-');
    xlabel('Sample count');
    ylabel('Rx samples');
    title('Received samples after downsampling');
    
end

%% Receiver logic starts here

% Process the preamble by first detecting the packet and then estimating
% the channel

%% Packet detection logic
[detect_pkt_flag LTFpeakPos] = detect_packet_rx(sigin);

% Only if packet is detected proceed
if ~detect_pkt_flag 
    BER = Inf;
	if VERBOSE1
        fprintf('No packet is detected!\n');
    end
    return;
end
    
if VERBOSE1
    fprintf('Packet is detected!\n');
end


%% Channel estimation logic 
% Get the channel estimation vector in frequency domain
[channel_vector dataPos] = channel_estimation_rx(LTFpeakPos, sigin);


%% Data decoding logic
outbits = decode_data_rx(LTFpeakPos, dataPos, channel_vector, sigin);


%% Visualize the data
mBER = sum(outbits ~= inbits)/length(inbits);    
BER = mBER;
    
if DEBUG_ON
        
    figure(202);
    stairs(inbits);
    axis([0 length(inbits) -0.5 1.5]);
    xlabel('Bit number');
    ylabel('Bit value');
    title('Input bits');
     
    figure(203);
    stairs(outbits);
    axis([0 length(outbits) -0.5 1.5]);
    xlabel('Bit number');
    ylabel('Bit value');
    title('Output bits');
     
    bit_error = xor(outbits, inbits);
    idx = (bit_error == 0);
    x = 1:length(bit_error);
    figure(204);
    plot(x(idx), bit_error(idx), 'g*', x(~idx), bit_error(~idx), 'ro');
    axis([0 length(bit_error) -0.5 1.5]);
    xlabel('Bit position');
    ylabel('Error');
    title('Error position');
     
end

end