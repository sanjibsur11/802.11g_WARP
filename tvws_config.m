% A configuration file to modify various parameters for the transmitter and
% receiver in TV white space WARP communications
% 
% Author: Sanjib Sur
% Institute: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 05/23/2014
% 
% Comments: Various parameters handling file to have configuration across
% transmitter and receivers. 
% 

%% Read global variables
tvwsGlobalVars;


%% Basic mode control
% Simulation control
USESIM = 0; % whether to use simulation or WARP board
% Debugging control
% Show constellation or not
DEBUG_ON = 1;
% First level of print outs
VERBOSE1 = 1;
% Second level and detailed print outs
VERBOSE2 = 1;
% Suitable for debugging without adding noise to the signal
useFakeChannel = 0; % Use a fake channel to Debug
% Whether to use trace driven channel estimation results
TRACE_DRIVEN = 0;
% Whether to collect trace of the channels to later use for emulation
TRACE_COLLECT = 0;


%% Timing parameters
% Timing related parameters for both OFDM 
Timeparams = [];
% Number of Data subcarriers
Timeparams.Nsd = 48;                             
% Number of pilot subcarriers
Timeparams.Nsp = 4;                              
% Number of DC subcarriers
Timeparams.Ndc = 1;
% Number of guardbands
Timeparams.Ngc = 11;
% Total subcarriers
Timeparams.Nst = Timeparams.Nsd + Timeparams.Nsp + Timeparams.Ndc + ...
    Timeparams.Ngc;
% TODO: OFDM sample rate
Timeparams.fs = 1;                                
% OFDM sample Period
Timeparams.ts = 1/Timeparams.fs;
% FFT size: Combine the number of subcarriers
Timeparams.Sfft = Timeparams.Nst;
% DFT duration
Timeparams.Tdft = Timeparams.Sfft * Timeparams.ts;
% Gaurd Interval
Timeparams.Tg = Timeparams.Tdft/4;
% Total Symbol Duration
Timeparams.Tsym = Timeparams.Tdft + Timeparams.Tg;


%% Transmitter and Receiver Configuration
Txparams = [];
Rxparams = [];

% Carrier frequency offset between Tx and Rx
iFreq = 0;
% PHY transmission modes: SISO - 0, Diversity - 1, Spatial Multiplexing -
% 2 ...
Txparams.PHYmode = 0;

% How many antennas to use?
if Txparams.PHYmode == 0
    Txparams.numAntenna = 1;
    Rxparams.numAntenna = 1;
    
elseif Txparams.PHYmode == 1
    Txparams.numAntenna = 2;
    Rxparams.numAntenna = 1;
    
elseif Txparams.PHYmode == 2 % For spatial multiplexing, number of Rx
    % antennas has to be atleast equal or more than Tx antennas
    Txparams.numAntenna = 2;
    Rxparams.numAntenna = 2;
    
end

% Number of data packets
Txparams.totalDataPkts = 1;
% Number of OFDM symbols per packet
Txparams.numOFDMsymb = 10;
% Total number of data symbols = Number of data subcarriers * OFDM symbols
Txparams.numSymbTotal = Timeparams.Nsd * Txparams.numOFDMsymb;
% Modulation size: BPSK = 2, QAM = 4, 8-QAM = 8, 16-QAM = 16 etc.
Txparams.modulationM = 2;

% Number of bits per symbols
symb2BitNum = log2(Txparams.modulationM);
% Number of parallel stream to send through Tx antennas
if Txparams.PHYmode == 2 % Only spatial multiplexing will allow more than
    % 1 stream
    Txparams.streamNum = 1;
    
else
    Txparams.streamNum = 1;
    
end
% Number of bits send by Tx
Txparams.numBitsTotal = Txparams.numSymbTotal * symb2BitNum * ...
    Txparams.streamNum;


%% Simulation parameters
% Do we need any oversampling?
if USESIM
    osamp = 1;
else
    % Oversampling is required when transmitting through WARP board
    osamp = 2;
end

% Intended simulation SNR per antennas
SIMSNR = [25];


%% Preamble parameters
% STF
% Frequency domain STF: STF is using 64 sub-carriers. Do we need to change
% it?
if Timeparams.Nst == 64
STF_freq = [0 0 0 0 0 0 0 0 1+1i 0 0 0 -1+1i 0 0 0 -1-1i 0 0 0 1-1i 0 0 ...
    0 -1-1i 0 0 0 1-1i 0 0 0 0 0 0 0 1-1i 0 0 0 -1-1i 0 0 0 1-1i 0 0 0 ...
    -1-1i 0 0 0 -1+1i 0 0 0 1+1i 0 0 0 0 0 0 0];
end
% Before ifft, ensure the zero-frequency on index 1; minus-frequency begins
% from half (so the samples near the half-freq should be a sequence of 0s)
STF_time = ifft(fftshift(STF_freq));
STF_t_short = STF_time(1:Timeparams.Nst/4);
repeat_STF_cnt = 10;
% Repeat short STF sequence repeat_STF times
STF_repeat = repmat(STF_t_short, 1, repeat_STF_cnt);
% Scale to span -1,1 range of DAC
STFscale = max([max(abs(real(STF_repeat))), max(abs(imag(STF_repeat)))]);
% Time domain STF 
STF = (1/STFscale)*STF_repeat;
STF = STF;

% LTF
% Frequency domain LTF: LTF is using 64 sub-carriers. How do we estimate
% channels for lesser number of subcarriers?
if Timeparams.Nst == 64
LTF_freq_bot = [0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 ...
    -1 1 -1 1 1 1 1];
LTF_freq_top = [1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 ...
    -1 1 1 1 1 0 0 0 0 0];
end
LTF_freq = [LTF_freq_bot 0 LTF_freq_top];
LTF_shift = fftshift(LTF_freq);
LTF_time = ifft(LTF_shift);
% Scale to span -1,1 range of DAC
LTFscale = max([max(abs(real(LTF_time))), max(abs(imag(LTF_time)))]);
LTF_time = (1/LTFscale)*LTF_time;
% Concatenate two long training symbols and add cyclic prefix: 50% cyclic
% prefix is added for LTF
LTF = [LTF_time(end/2 + 1:end) repmat(LTF_time, 1, 2)];

% VLTF
% Frequency domain VLTF

% Mapping subcarrier format for precoding to Timeparams.Nsd data
% subcarriers ID
if Timeparams.Nst == 64
    data_to_subcarrier_map = [ 0     0     0     0     0     0    25 ...
        26    27    28    29     0    30    31    32    33    34    35 ...
        36    37    38    39    40    41    42     0    43    44    45 ...
        46    47    48     0     1     2     3     4     5     6     0 ...
        7     8      9    10    11    12    13    14    15    16    17 ...
        18    19     0    20    21    22    23    24     0     0     0 ...
        0     0];
        % Left and Right guard band length in terms of sub-carriers
        left_guard_len = 6;
        right_guard_len = 5;
        % Pilot position
        pilot_position = [12 26 40 54];
        
elseif Timeparams.Nst == 32
elseif Timeparams.Nst == 16
elseif Timeparams.Nst == 8
    
end


%% Packet detection parameters
STFcorrThrsh = 0.0005;


%% Initialization of WARP radio parameters
numTxNode = 1;
numRxNode = 1;

% We need to initialize separately if we are using WURC daughter board
if ~USESIM
    
    WARPLab_TxDelay = 1000;
    % WARPLab buffer size 32k
    WARPLab_TxLength = 2^15-WARPLab_TxDelay; % Length of transmission. In [0:2^15-1-TxDelay]
    WARPLab_CarrierChannel = 13; % Channel in the 2.4 GHz band. In [1:14] (avoid...
                         % 1 to 11); 5GHz in [15:37]
    WARPLab_TxGain_RF = 35; % Tx RF Gain. In [0:63] 
    WARPLab_TxGain_BB = 1; % Tx Baseband Gain. In [0:3]
    WARPLab_RxGain_BB = 10; % Rx Baseband Gain. In [0:31]
    WARPLab_RxGain_RF = 1; %2; % Rx RF Gain. In [1:3]
    
	USE_AGC = false;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the WARPLab experiment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Create a vector of node objects
    nodes = wl_initNodes(numRxNode + numTxNode);
    WARPLab_node_tx = nodes(1);
    WARPLab_node_rx = nodes(2:length(nodes));
    %{
    nodes = wl_initNodes(3);
    node_tx = nodes(2);
    node_rx = nodes(3);
    %}

    %Create a UDP broadcast trigger and tell each node to be ready for it
    WARPLab_eth_trig = wl_trigger_eth_udp_broadcast;
    wl_triggerManagerCmd(nodes,'add_ethernet_trigger',[WARPLab_eth_trig]);

    %Get IDs for the interfaces on the boards. Since this example assumes each
    %board has the same interface capabilities, we only need to get the IDs
    %from one of the boards
    [RFA,RFB] = wl_getInterfaceIDs(nodes(1));

    % Txparams.PHYmode: 0 SISO
    if Txparams.PHYmode == 0
        WARPLab_RF_vector = [RFA];
    else
    end

    %Set up the interface for the experiment
    wl_interfaceCmd(nodes,sum(WARPLab_RF_vector),'tx_gains',WARPLab_TxGain_BB,WARPLab_TxGain_RF);
    wl_interfaceCmd(nodes,sum(WARPLab_RF_vector),'channel',2.4,WARPLab_CarrierChannel);

    if(USE_AGC)
        wl_interfaceCmd(nodes,sum(WARPLab_RF_vector),'rx_gain_mode','automatic');
        wl_basebandCmd(nodes,'agc_target',-10);
        wl_basebandCmd(nodes,'agc_trig_delay', 500);
        wl_basebandCmd(nodes,'agc_dco', true);
    else
        wl_interfaceCmd(nodes,sum(WARPLab_RF_vector),'rx_gain_mode','manual');
        RxGainRF = 1; %Rx RF Gain in [1:3]
        RxGainBB = 15; %Rx Baseband Gain in [0:31]
        wl_interfaceCmd(nodes,sum(WARPLab_RF_vector),'rx_gains',WARPLab_RxGain_RF,WARPLab_RxGain_BB);
    end


    %We'll use the transmitter's I/Q buffer size to determine how long our
    %transmission can be
    %txLength = nodes(1).baseband.txIQLen;

    %Set up the baseband for the experiment
    wl_basebandCmd(nodes,'tx_delay',WARPLab_TxDelay);
    wl_basebandCmd(nodes,'tx_length',WARPLab_TxLength);
    
end


%% Padding
padding_size = 64;
for txantn = 1:Txparams.numAntenna
    Padding(txantn, :) = zeros(1, padding_size);
end


%% Modulation parameters
% M-QAM modulation using communication toolbox
hmodem = modem.qammod('M', Txparams.modulationM);
% M-QAM demodulation using communication toolbox
hdemodem = modem.qamdemod('M', Txparams.modulationM);


%% Channel estimation parameter
usePilotComp = 1;
% Considering carrier frequency of 2.4 GHz
Fc = 2400000;
measured_CFO = [];


%% Parameters check