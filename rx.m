function [rxbits conf] = rx(rxsignal,conf,k)
% Digital Receiver
%
%   [txsignal conf] = tx(txbits,conf,k) implements a complete causal
%   receiver in digital domain.
%
%   rxsignal    : received signal
%   conf        : configuration structure
%   k           : frame index
%
%   Outputs
%
%   rxbits      : received bits
%   conf        : configuration structure
%

f_c = conf.f_c;
f_s = conf.f_s;
nbits = conf.nbits;
modulation_order = conf.modulation_order;
os_factor = conf.os_factor;
npreamble = conf.npreamble;

% dummy 
% rxbits = zeros(conf.nbits,1);

% down conversion / de-multiplexing
carrier_seq = exp(-1i*2*pi*(f_c/f_s)*(1:length(rxsignal))).';
rxsyms_dc = rxsignal.*carrier_seq;
% low pass filter
rxsyms_lp = 2*lowpass(rxsyms_dc, conf);
% matched filter
rxsyms_mf = conv(rxsyms_lp, rrc(os_factor), 'same');
% frame synchronizationg
[data_idx, theta, magnitude] = frame_sync(rxsyms_mf, conf);

% % time estimation
% % Use preamble symbols to improve timing offset estimation
% cum_err = 0;
% for i=1:npreamble
%      start_ind = data_idx-npreamble+(i-1)*os_factor;
%      idx_range = start_ind:start_ind+os_factor-1;
%      segment = rxsyms_mf(idx_range);
%      cum_err = cum_err + exp(-1i*(1*pi)*(0:os_factor-1)*(1/os_factor))*(abs(segment).^2);
% end

% data_t = zeros(ceil(nbits/modulation_order), 1);
% for ii=1:ceil(nbits/modulation_order)
%     
%      start_idx  = data_idx+(ii-1)*os_factor;
%      idx_range  = start_idx:start_idx+os_factor-1;
%      segment    = rxsyms_mf(idx_range);
%      % Estimate timing error epsilon
%      [epsilon_tmp, cum_err] = epsilon_estimate(segment, cum_err, conf);
%      % epsilon(ii) = epsilon_tmp;
%      data_t(ii) = interpolator(rxsyms_mf, epsilon_tmp, start_idx, conf);
% end

data_t = rxsyms_mf(data_idx:os_factor:data_idx+os_factor*nbits/modulation_order-1);
% phase estimation
data_p = data_t*exp(-1i*theta)/magnitude;
norm_data_p = data_p/(sqrt((data_p')*data_p/1000));

rxbits = demodulator(norm_data_p, modulation_order);
end

function [epsilon_tmp, cum_err] = epsilon_estimate(segment, cum_err, conf)
os_factor = conf.os_factor;
tmp_seq = exp(-1i*(1*pi)*(0:os_factor-1)*(1/os_factor));
cum_err = cum_err + tmp_seq*segment;
epsilon_tmp = -angle(cum_err)/(2*pi);
end

function [beginning_of_data, phase_of_peak, magnitude_of_peak] = frame_sync(in_syms, conf)

npreamble = conf.npreamble;
os_factor = conf.os_factor;

preamble_syms = modulator(preamble_generate(npreamble), 1);
current_peak_value = 0;
samples_after_threshold = os_factor;
detection_threshold = 15;
% corVal = zeros(length(rxsyms_mf)-os_factor*npreamble, 1);

for i = os_factor*npreamble+1:length(in_syms)
    r = in_syms(i-os_factor*npreamble:os_factor:i-os_factor); 
    c = preamble_syms'*r;
    T = abs(c)^2/abs((r')*r);
    
    % corVal(i-os_factor*npreamble) = T; 
    
    if (T > detection_threshold || samples_after_threshold < os_factor)
        samples_after_threshold = samples_after_threshold - 1;
        if (T > current_peak_value)
            beginning_of_data = i;
            phase_of_peak = mod(angle(c),2*pi);
            magnitude_of_peak = abs(c)/npreamble;
            current_peak_value = T;
        end
    end
end
end

function [out_sym] = interpolator(in_syms, epsilon, start_idx, conf)
os_factor = conf.os_factor;
interpolator_type = conf.interpolator_type;
offset = floor(epsilon*os_factor);
remain = os_factor*epsilon - offset;
switch interpolator_type
    case 'none'
        out_sym = in_syms(start_idx+round(epsilon*os_factor));
    case 'linear'
        out_sym = (1-remain)*in_syms(start_idx+offset) ...
         + remain*in_syms(start_idx+offset+1);
    case 'cubic'
        idx_range = start_idx+offset-1:start_idx+offset+2;
        segment = in_syms(idx_range);
        A = [-1 1 -1 1; 0 0 0 1; 1 1 1 1; 8 4 2 1];
        v = A\segment;
        temp = [remain^3, remain^2, remain, 1].';
        out_sym = (v.')*temp;
    otherwise
        out_sym = insyms(start_idx+round(epsilon*os_factor));
end
end