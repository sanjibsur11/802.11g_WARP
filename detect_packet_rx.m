% Detection of packets from STF im OFDM
% 
% Author: Sanjib Sur
% Institute: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 05/22/2014
% 
% Comments: This file contains the packet detection logic for OFDM PHY
% layer. Irrespective of the transmission mode being used, the packet
% detection logic remains same.
% 

function [detect_pkt_flag LTFpeakPos] = detect_packet_rx(sigin)


%% Read global variables 
tvwsGlobalVars;


%% Packet detection starts here

skipsamples = 50;
LTFpeakPos = 1;
ns = skipsamples;

% STF correlation from the sequence using the sliding window correlator
STF_correlator_output = [];
% Instanteneous signal energy
energy_output = [];
% To detect how many number of peaks went beyond the STFcorrThrsh value
STFpeakCnt = 0;
% Whether packet is detected or not
detect_pkt_flag = 0;


%% Packet detection from STF

% =========== Initialize parameters ===========
numPkt = 0;
smsEng_windowSize = (Timeparams.Nst/4);
smsEng = 0; %smoothed energy level
smsEngQ = zeros(Timeparams.Nst, 1);

SNRQsize = Timeparams.Nst;
SNRQ = zeros(SNRQsize, 1);
SNRthresh = 4;
SNRQthresh = 0.6;

selfcorrQsize = (Timeparams.Nst/4);
selfcorrQ = zeros(selfcorrQsize, 1);
selfcorrQthresh = 0.8; 

prevPktPos = 0;
corr2Engthresh = 0.9;
avgCrosscorr = zeros((Timeparams.Nst/4), 1);

if USESIM
    for k = 1:length(smsEngQ)
        smsEngQ(k) = 1e-8;
    end
end

while ns < length(sigin) - (Timeparams.Nst/4) % Corrlation length of 
    % (Timeparams.Nst/4)
    
    % Start locating LTF
    selfcorr = sum(sigin(ns + 1 : ns + (Timeparams.Nst/4)).* ...
        conj(sigin(ns - (Timeparams.Nst/4) + 1 : ns)));
    abs(selfcorr);
    energy = sum(abs(sigin(ns + 1 : ns + (Timeparams.Nst/4))).^2);
    
    STF_correlator_output = [STF_correlator_output selfcorr];
    energy_output = [energy_output energy];

    % --- Update smoothed energy level and SNR queue -----
    smsEng = smsEng*(1 - 1/smsEng_windowSize) + ...
        1/smsEng_windowSize*abs(sigin(ns + 1))^2;
    
    if (USESIM && abs(smsEng) <= 0)
        smsEng = 1e-8;
    end 
    
    %sliding window
    smsEngQ(end + 1) = smsEng;
    smsEngQ(1) = [];
    
    if (smsEngQ(1) <= 0)
        SNR = 0;
    else
        SNR = 10*log10(smsEng / smsEngQ(1));
    end
    
    if (SNR > SNRthresh && energy > STFcorrThrsh)
        SNRQ(end + 1) = 1;
    else
        SNRQ(end + 1) = 0;
    end
    
    SNRQ(1) = [];

    maxSTFcrossCorrPeak = 0;
    
    % ----- Update self-correlation queue -------
%     if ((ns - prevPktPos) > 2 * SNRQsize && ...
%         (abs(selfcorr) / energy > corr2Engthresh && ...
%         abs(selfcorr) / energy < 1 / corr2Engthresh))
%         selfcorrQ(end + 1) = 1;
    if ((abs(selfcorr) / energy > corr2Engthresh && ...
        abs(selfcorr) / energy < 1 / corr2Engthresh))
        selfcorrQ(end + 1) = 1;
    else
        selfcorrQ(end + 1) = 0;
    end
    selfcorrQ(1) = [];

    % ---- Decision: self-correlation and energy detection ----
    if ((sum(SNRQ(end - (Timeparams.Nst/4):end))) / ((Timeparams.Nst/4)) < ...
            SNRQthresh || sum(selfcorrQ) / length(selfcorrQ) < ...
            selfcorrQthresh)
        %keep sliding
        ns = ns + 1;
        continue;
    end

    % ---- Self-corr peak detected ---- 
	if (VERBOSE2)
    	fprintf('!!! Self-corr peak detected at sample %d level %f\n', ...
            ns, energy);
	end
    firstSTFpos = ns; %STF located
	selfcorrPeak = energy;

    % ---- Search for the first peak via cross-corr with STF ----
    af = 1; 
    peakPos = 0; 
    
    for k = 1:(Timeparams.Nst/4)
        crosscorrSTF = sum(sigin(ns + k : ns + k + (Timeparams.Nst/4)-1).*...
            conj(STF_time(1 + 1 : (Timeparams.Nst/4) + 1)));
        
        
        if k==1
            avgCrosscorr(1) = abs(crosscorrSTF);
            continue;
        end
        avgCrosscorr(k) = af*abs(crosscorrSTF) + (1-af)*avgCrosscorr(k-1);
    end
    
    for k=1 + 1 : (Timeparams.Nst/4) - 1
        % detect a peak in the avgCrosscorr curve
        if (avgCrosscorr(k) > maxSTFcrossCorrPeak &&...
               avgCrosscorr(k) > avgCrosscorr(k - 1) &&...
               (k<4||avgCrosscorr(k) > avgCrosscorr(k - 3)) &&...
               (k<6||avgCrosscorr(k) > avgCrosscorr(k - 5)) &&...
               avgCrosscorr(k) > avgCrosscorr(k + 1) &&...
               (k>(Timeparams.Nst/4) - 3 || avgCrosscorr(k) > avgCrosscorr(k + 3)) &&...
               (k>(Timeparams.Nst/4) - 5 || avgCrosscorr(k) > avgCrosscorr(k + 5)))
                
            maxSTFcrossCorrPeak = avgCrosscorr(k);
            peakPos = k;

            if VERBOSE2
                fprintf(1, 'peakPos %d, maxSTFcrossCorrPeak %f\n', ...
                       peakPos, maxSTFcrossCorrPeak);
            end

        end
    end
    
    ns = ns + peakPos;


    % ------- Detect other peaks ---------
	endSTF = 0;
    peakPos = 0;
    % Calculate avgCrosscorr of all the following sample positions
    numPeak = 1;
    psamples = ((Timeparams.Nst/4))*(repeat_STF_cnt - 1);
    
    for k = 1:psamples %max number of peaks is repeat_STF_cnt - 1
        crosscorrSTF = sum(sigin(ns+k:ns+k+((Timeparams.Nst/4))-1).*...
            conj(STF_time(1 + 1 : (Timeparams.Nst/4) + 1)));       
        if k==1
            avgCrosscorr(1) = abs(crosscorrSTF);
            continue;
        end
        avgCrosscorr(k) = af*abs(crosscorrSTF) + (1-af)*avgCrosscorr(k-1);
    end 
    
    for k=2:psamples 
        if (avgCrosscorr(k) > maxSTFcrossCorrPeak*0.8 &&...
           avgCrosscorr(k) > avgCrosscorr(k-1) && ... 
           (k<3||avgCrosscorr(k) > avgCrosscorr(k-2)) && ... 
           (k<4||avgCrosscorr(k) > avgCrosscorr(k-3)) &&...
           (k<6||avgCrosscorr(k) > avgCrosscorr(k-5)) &&...
           (k>psamples-1 || avgCrosscorr(k) > avgCrosscorr(k+1)) &&...
           (k>psamples-2 || avgCrosscorr(k) > avgCrosscorr(k+2)) &&...
           (k>psamples-3 || avgCrosscorr(k) > avgCrosscorr(k+3)) &&...
           (k>psamples-5 || avgCrosscorr(k) > avgCrosscorr(k+5)))
   
            numPeak = numPeak + 1;
            peakPos = k;
            if VERBOSE2
                fprintf(1, 'peak at %d, avgCrosscorr=%g\n', ...
                       ns+k, avgCrosscorr(k));
            end
        end

        if k-peakPos > 18 % no more peaks
            if VERBOSE2
                fprintf(1, 'No more peaks at %d\n', ns+k);
            end
            break;
        end
         
		
		% End of STF: self-corr ends (if there is a peak nearby, the distance will not be greater than 16)
        tempselfcorr = sum(sigin(ns+k+1:ns+k+(Timeparams.Nst/4)).*conj(sigin(ns+k-(Timeparams.Nst/4)+1:ns+k)));
        tempselfcorr = abs(tempselfcorr);
        tempEng = sum(abs(sigin(ns+k+1:ns+k+(Timeparams.Nst/4))).^2);
        if (tempselfcorr/tempEng>corr2Engthresh && ...
            tempselfcorr/tempEng<1/corr2Engthresh )
            selfcorrQ(end+1) = 1;
        else
            selfcorrQ(end+1) = 0;
        end
        selfcorrQ(1) = [];

        if (sum(selfcorrQ(end-(Timeparams.Nst/4)+1:end))/(Timeparams.Nst/4) < selfcorrQthresh)
            if VERBOSE2
               fprintf(1, 'case 1 End of STF at %d\n', ns + k);
            end
			endSTF = ns + k + (Timeparams.Nst/4);
			break; 
        end		
        
        if (ns+k-firstSTFpos > (2*Timeparams.Nst) + (Timeparams.Nst/2)) % 50% 
            % cycling prefix for LTF
            if VERBOSE2
                fprintf(1, 'case 2 End of STF at %d\n', ns+k);
            end
            endSTF = ns+k;
            break;
        end
        
    end %end for k=2:psamples
    
    if (numPeak < 2) 
            SNRQ = zeros(SNRQsize,1);
            selfcorrQ = zeros(selfcorrQsize,1); 
            fprintf(1, 'pkt lost, num of detected STF peaks is %d\n', numPeak);
            continue; 
    end
    
    detect_pkt_flag = 1;
    
    %STF ends, start point for LTF detection
    ns = endSTF;   
    
    % Get the exact position of LTF via crosscorr    
    for t = ns+10:ns+30        
        crosscorrLTF(t - ns - 9) = abs(sum(sigin(t+1:t+(Timeparams.Nst/4)).*...
            conj(LTF_time(1:(Timeparams.Nst/4)))));
    end   
    
% 	for t = ns:ns+20
%         crosscorrLTF(t - ns + 1) = abs(sum(sigin(t+1:t+(Timeparams.Nst/4)).*...
%             conj(LTF_time(1:(Timeparams.Nst/4)))));
%     end   
    
	[~, peakPos] = max(crosscorrLTF(1:end));

    LTFpeakPos = ns + 10 - 1 + peakPos;
%     LTFpeakPos = ns + peakPos;
    break;
    
end

if DEBUG_ON
    
    % Plotting STF correlation output
    length_STF = length(STF_correlator_output);
    figure(500);
    plot(1:length_STF, real(STF_correlator_output), 'b-', ...
        1:length_STF, imag(STF_correlator_output), 'r-');
    xlabel('Sample count');
    ylabel('STF correlation samples');
    title('STF correlation output');
    
    % Plotting instanteneous energy data
    
    figure(501);
    plot(energy_output);
    xlabel('Sample count');
    ylabel('Energy');
    title('Instanteneous energy');
    
end


end