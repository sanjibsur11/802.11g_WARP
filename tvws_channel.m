% Channel model for TVWS
% 
% Author: Sanjib Sur
% Institute: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 05/21/2014
% 
% Comments: This file contains the channel model for TV white space
% communications
% 
% 1. Fake channel: Extremely helpful for initial debugging and
% implementation of the PHY layer
% 
% 2. Additive Gaussian: Useful for simulations etc. and to calculate
% bit-error rate in simulation
% 
% 3. Trace-based channel: Very useful when we want to perform trace-based
% experiments. First, we can collect trace channel responses from WARP
% board and then apply the trace responses to do trace-based emulation
% 
% 4. WARP board: This is the main 60 GHz WARP transmit and receive channel
% 
% 
% 

function [RxData_board] = tvws_channel(TxData_board)

%% Read global variables 
tvwsGlobalVars;


%% Transmit through channel 

if USESIM % Use only simulation
    if useFakeChannel % Just bypass everything, easy to debug
        RxData_board = TxData_board;
        
    elseif TRACE_DRIVEN % If trace driven then apply the channel traces 
        % from the channel_trace_file
        
    else % Only add additive Gaussian noise
        
        for txantn = 1:Txparams.numAntenna
            RxData_board(txantn, :) = awgn(TxData_board(txantn, :), SIMSNR(txantn));
        end
        
    end
else % Use WARP board for transmission
    
    if length(TxData_board) < WARPLab_TxLength
        TxData_board = [TxData_board zeros(1, WARPLab_TxLength - ...
            length(TxData_board))];
    end
    
    % assign the data to transmit
    wl_basebandCmd(WARPLab_node_tx(1), WARPLab_RF_vector, 'write_IQ', TxData_board.');

    % start transmission and reception
    wl_interfaceCmd(WARPLab_node_tx, sum(WARPLab_RF_vector),'tx_en');
    wl_interfaceCmd(WARPLab_node_rx, sum(WARPLab_RF_vector),'rx_en');

    wl_basebandCmd(WARPLab_node_tx, sum(WARPLab_RF_vector),'tx_buff_en');
    wl_basebandCmd(WARPLab_node_rx, sum(WARPLab_RF_vector),'rx_buff_en');

    WARPLab_eth_trig.send();

    % get the received data
    RxData_board = wl_basebandCmd(WARPLab_node_rx(1), WARPLab_RF_vector, 'read_IQ', 0,...
                                                            WARPLab_TxLength+WARPLab_TxDelay);
    RxData_board = RxData_board.';
end

end