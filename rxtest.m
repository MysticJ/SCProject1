function [rxbits conf] = rxtest(rxsignal, txsyms,conf,k)
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

rxsyms_mf = preprocess(rxsignal, conf);
% frame synchronizationg
[data_idx, theta] = frame_sync(rxsyms_mf, conf);

% time estimation
% Use preamble symbols to improve timing offset estimation
cum_err = 0;
for i=1:npreamble
     start_ind = data_idx-npreamble*os_factor+(i-1)*os_factor;
     idx_range = start_ind:start_ind+os_factor-1;
     segment = rxsyms_mf(idx_range);
     cum_err = cum_err + exp(-1i*(2*pi)*(0:os_factor-1)*(1/os_factor))*(abs(segment).^2);
end

data_time = zeros(ceil(nbits/modulation_order), 1);
data_phase = zeros(ceil(nbits/modulation_order), 1);
epsilon_arr = zeros(ceil(nbits/modulation_order), 1);
theta_hat = zeros(ceil(nbits/modulation_order)+1, 1); 
theta_hat(1) = theta;

for ii=1:ceil(nbits/modulation_order)
    
     start_idx  = data_idx+(ii-1)*os_factor;
     idx_range  = start_idx:start_idx+os_factor-1;
     segment    = rxsyms_mf(idx_range);
     % Estimate timing error epsilon
     [epsilon_tmp, cum_err] = epsilon_estimate(segment, cum_err, conf);
     epsilon_arr(ii) = epsilon_tmp;
     data_time(ii) = interpolator(rxsyms_mf, epsilon_tmp, start_idx, conf);
     
     deltaTheta = 1/4*angle(-data_time(ii)^4) + pi/2*(-1:4);
    
     [~, ind] = min(abs(deltaTheta - theta_hat(ii)));
     theta = deltaTheta(ind);

     % Lowpass filter phase
     theta_hat(ii+1) = mod(0.3*theta + 0.7*theta_hat(ii), 2*pi);
    
     % Phase correction and rotate the current symbol accordingly
     data_phase(ii) = data_time(ii) * exp(-1i * theta_hat(ii+1));
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the magnitude and phase of the estimated and original symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% txsyms_mf = preprocess(txsyms, conf);
% txsyms_sample = txsyms_mf(os_factor*npreamble+1:os_factor:length(txsyms_mf));
% figure(4); 
% subplot(2, 1, 1);
% plot(abs(txsyms_sample(1:100)), 'Linewidth', 1.5);hold on; 
% plot(abs(data_phase(1:100)), 'Linewidth', 1.5); hold off;
% ylabel('magnitude'); legend('txsyms', 'rxsyms after correction');
% subplot(2, 1, 2);
% plot(angle(txsyms_sample(1:100)), 'Linewidth', 1.5); hold on; 
% plot(angle(data_phase(1:100)), 'Linewidth', 1.5); hold off;
% ylabel('phase'); legend('txsyms', 'rxsyms after correction');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% phase estimation
data_p = data_phase; % *exp(-1i*theta)/magnitude;
norm_data_p = data_p/(sqrt((data_p')*data_p/(nbits/modulation_order)));

rxbits = demodulator(norm_data_p, modulation_order);
% rxbits = data_t;
end

function [out_syms] = preprocess(in_syms, conf)
f_c = conf.f_c; f_s = conf.f_s; os_factor = conf.os_factor;
carrier_seq = exp(-1i*2*pi*(f_c/f_s)*(1:length(in_syms))).';
in_syms_dc = in_syms.*carrier_seq;
in_syms_lp = 2*lowpass(in_syms_dc, conf);
out_syms = conv(in_syms_lp, rrc(os_factor), 'same');
end

function [epsilon_tmp, cum_err] = epsilon_estimate(segment, cum_err, conf)
os_factor = conf.os_factor;
tmp_seq = exp(-1i*(2*pi)*(0:os_factor-1)*(1/os_factor));
cum_err = cum_err + tmp_seq*(abs(segment).^2);
epsilon_tmp = -angle(cum_err)/(2*pi);
end

function [beginning_of_data, phase_of_peak] = frame_sync(in_syms, conf)

npreamble = conf.npreamble;
os_factor = conf.os_factor;

preamble_syms = modulator(preamble_generate(npreamble), 1);
current_peak_value = 0;
samples_after_threshold = os_factor;
detection_threshold = 15;
% corVal = zeros(length(in_syms)-os_factor*npreamble, 1);

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
            % magnitude_of_peak = abs(c)/npreamble;
            current_peak_value = T;
        end
    end
end
% figure(3);plot(corVal); ylabel('covariance with preamble')
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