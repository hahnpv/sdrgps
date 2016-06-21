function ca_code_value = ca_code_nco_v2(PRN_number, code_phase, strType)
% Code phase is specified as a number with 10 fractional bits: round(phi*2^10), where phi is a
% number from 0 (first chip) to 1022 (last chip).  The code phase will be
% interpreted modulo 1023*2^10.
% Uses linear interpolation between the exact chip values if strType ==
% 'linear', and square bits if strType == 'square'

if nargin == 0
    % no input argument means to run unit tests
    fprintf(1, 'ca_code_nco_v2(): Called with no arguments: running unit tests...\n');
    ca_code_nco_v2_unit_tests;
    return;
end
if nargin < 3
    strType = 'linear';
end

N_fractional = 10;  % code phase is scaled by 2^10
N_chips = 1023;     % Every code phase will be interpreted as modulo 1023
code_phase = code_phase(:); % Ensure that we get only column vectors

persistent ca_codes
if isempty(ca_codes)
    fprintf(1, 'generating prns...\n');
    tic;
    % Generate CA codes for all PRNs (1023 chips, PRNs 1-37):
    PRN_ar = (1:37).';
    ca_codes = zeros(N_chips, length(PRN_ar));
    for kPRN = 1:length(PRN_ar)
        ca_codes(:, kPRN) = cacode(PRN_ar(kPRN), 1).';
        ca_codes(:, kPRN) = ca_codes(:, kPRN)*2-1;  % convert to bipolar (+1/-1)
    end
    toc;
end


% tic;
% code_phase_axis = (0:1022).'*2^N_fractional;
% ca_code_value = interp1(code_phase_axis, ca_code, mod(code_phase, 2^10*1023), '*linear');
% toc;
% 

if strcmpi(strType, 'linear')
    % code with linear bit transitions (less aliasing)
%     tic;
    code_phase_before = mod(floor(code_phase/2^N_fractional), N_chips);
    code_phase_after = mod(code_phase_before+1, N_chips);
    fract_part = (code_phase - floor(code_phase/2^N_fractional)*2^N_fractional)/2^N_fractional;
    try
        ca_code_value = ca_codes(code_phase_before+1, PRN_number) .* (1-fract_part) + ca_codes(code_phase_after+1, PRN_number) .* (fract_part);
    catch E
        keyboard
    end
        
%     toc;
else
    % Square wave code
%     tic;
    code_phase_quant = mod(round(code_phase/2^N_fractional), N_chips);
    ca_code_value = ca_codes(code_phase_quant+1, PRN_number);
%     toc;
end




% figure;
% plot(code_phase_quant/N_chips);
% hold all;
% plot(2.1 + ca_code_value);


end

function ca_code_nco_v2_unit_tests()
    % generate the same PRN twice at slightly different rates and show the
    % resulting beating after filtering
    PRN_test = 22;
    N_oversampling = 100;
    N_sequences = 1;
    N_chips = 1023;
    nco_phase1 = (0:1/N_oversampling:N_chips*N_sequences-1/N_oversampling).';
    nco_phase2 = (0:1/N_oversampling:N_chips*N_sequences-1/N_oversampling).' + 0.25;   % this one is just phase advanced by 90 deg
    
    nco_output1 = ca_code_nco_v2(PRN_test, 2^10*nco_phase1, 'square');
    nco_output2 = ca_code_nco_v2(PRN_test, 2^10*nco_phase2, 'square');
    
    figure;
    plot(nco_output1, '.-');
    hold all;
    plot(nco_output2, '.-');
    title('Two NCOs vs time, same rate, green has a phase advance of 0.25 bit');
    
    % different rates
    PRN_test = 22;
    N_oversampling = 1;
    N_chips = 1023;
    N_sequences = 2*N_chips;
    nco_phase1 = (0:1/N_oversampling:N_chips*N_sequences-1/N_oversampling).';
    nco_phase2 = (0:1/N_oversampling:N_chips*N_sequences-1/N_oversampling).' * (N_sequences/(N_sequences-2));   % this one is just phase advanced by 90 deg
    
    nco_output1 = ca_code_nco_v2(PRN_test, 2^10*nco_phase1);
    nco_output2 = ca_code_nco_v2(PRN_test, 2^10*nco_phase2);
    
    beat = smooth(nco_output1 .* nco_output2, round(N_chips*N_oversampling));
    
    figure;
%     ax = subplot(211);
    plot(nco_output1, '-');
    hold all;
    plot(nco_output2, '-');
    plot(beat, '-');
    plot((nco_phase2-nco_phase1)/N_chips, '-');
    grid on;
    title('Two NCOs vs time and phase difference, green is faster');
    
    
end