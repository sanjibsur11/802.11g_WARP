% Implementation of TV white communication in UHF band using WURC + WARP
% boards
% 
% Author: Sanjib Sur
% Institute: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 05/22/2014
% 
% Comments: This file contains the SISO, STBC and Spatial Multiplexing PHY
% layer implmentation specifically for TV white space communication with
% using WURC daughterboard on WARP board. The file contains code for the
% transmitter side. For receiver side, check tvws_rx.m. OFDM is working
% fine in simulation. Need to debug WARP board. Change the WARP board
% oversampling rate to 2. Still the bits at the edges are corrupted.
% 
% % Version: 0.0.2
% Last modified: 05/23/2014
% 
% Comments: Fixed the issue of ifft and fft. Right now OFDM is working 
% perfectly in the WARP board. MATLAB ifft and fft functions assumes that 
% DC is at the index 1 and all the guard bands are in the middle of the 
% signal array. Oversampling by 2 or more is required. Tested with 1000
% packets in WARP, each packet containing 10 OFDM symbols and it is working
% fine.
% 

clc;
close all;
clear all;


%% Read global variables
tvwsGlobalVars;


%% Read configuration file
tvws_config;


%% Construct packets & samples
fg = fopen('databits.dat', 'w');
for k = 1:Txparams.numBitsTotal
	if (rand < 0.5)
        fprintf(fg, '1 ');
    else
        fprintf(fg, '0 ');
    end
    fprintf(fg, '\n');
end
fclose(fg);

% Input data bits
inbits = load('databits.dat');


%% Generate time domain preambles
Preamble = generate_preamble();


%% Generate time domain data signals
% Right now sending nothing in the payload
% Payload = [];
Payload = generate_payload();


%% Packet generation
if USESIM
    Packet = [Padding Preamble Payload];
else
    Packet = [Preamble Payload];
end


%% Transmitting samples

TxData_board = Packet;
% Upsample the complex samples according to the oversampling rate
for txantn = 1:Txparams.numAntenna
    new_TxData_board(txantn, :) = interp(TxData_board(txantn, :), osamp);
end
TxData_board = new_TxData_board;


%% Visualize data
if DEBUG_ON
    
    % Raw samples
    length_tx = length(TxData_board);
    figure(100);
    plot(1:length_tx, real(TxData_board), 'b-', ...
         1:length_tx, imag(TxData_board), 'r-');
    xlabel('Sample count');
    ylabel('Tx samples');
    title('Transmitted samples');
    
%     % Transmitter FFT
%     sigdft = fft(TxData_board);
%     sigdft = sigdft(1:length_tx/2+1);
%     psdx = (1/(Fs*length_tx)).*abs(sigdft).^2;
%     psdx(2:end-1) = 2*psdx(2:end-1);
%     freq = 0:Fs/length(TxData_board):Fs/2;
%     figure(101);
%     plot(freq, 10*log10(psdx)); 
%     title('PSD for transmitter samples');
    
    
end


%% ================= Start transmission rounds =====================

accBER = 0.0;
BER = 0.0;
correctDataPacket = 0;
header_lost = 0;
numDataPkt = 0;

while numDataPkt < Txparams.totalDataPkts
    
new_TxData_board = [];    
fprintf(1, '==== Transmission start for %d-th DATA packet ====\n',...
                                                    numDataPkt+1);
                                                
%% Transmit through channel 
RxData_board = tvws_channel(TxData_board);


%% Receiver side
preambleFound = 0;
[preambleFound, BER] = tvws_rx(RxData_board);

% Check whether the packet has correct preamble, otherwise discard the BER
% calculation
if preambleFound
    if BER < 0.05 % Consider packets which has BER less than 5%, otherwise,
        accBER = accBER + BER;
        correctDataPacket = correctDataPacket + 1;
    end
else
    header_lost = header_lost + 1;
end

numDataPkt = numDataPkt + 1;

fprintf('\n\n');


end     

%% Performance metrics

if correctDataPacket > 0
    fprintf ('Average BER = %0.6g%%\n', (accBER/correctDataPacket)*100);
    fprintf('Packet loss rate = %0.6g%%\n', ...
        ((numDataPkt-correctDataPacket)/numDataPkt)*100);
    fprintf('Packet header loss rate = %0.6g%%\n', (header_lost/numDataPkt)*100);
    fprintf('Packet loss due to BER = %0.6g%%\n', ...
        ((numDataPkt-correctDataPacket-header_lost)/numDataPkt)*100);
else
    fprintf('Nothing has been decoded!\n');
end


%% Visualize CFO
if DEBUG_ON
    figure(102);
	plot(measured_CFO, ...
	'-or', 'LineWidth', 2, 'MarkerSize', 10); 
    xlabel('Packet number');
    ylabel('kHz');
    title('CFO between Tx and Rx');
end


% Reset and disable the WARP board

if ~USESIM
    wl_basebandCmd(nodes,sum(WARPLab_RF_vector),'tx_rx_buff_dis');
    wl_interfaceCmd(nodes,sum(WARPLab_RF_vector),'tx_rx_dis');
end

fprintf('\n\n');