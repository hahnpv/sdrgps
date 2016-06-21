% bInterleaved = true indicates that data_baseband is a vector containing
% interleaved real and imaginary parts instead of a complex vector
function [timestamp code_phase doppler vector_prompt vector_early vector_late carrier_phase_out phase_error vector_promptp1 vector_promptm1 vector_prompt0] = ...
    ca_code_tracking_v1(data_baseband, bInterleaved, PRN_number, fs, initial_doppler, initial_code_phase, start_offset_in_points, stop_offset_in_points, ...
    gain_t, gain_k)

if nargin == 0
    fprintf(1, 'ca_code_tracking_v1_unit_tests(): Called with no arguments: running unit tests...\n');
    ca_code_tracking_v1_unit_tests();
    return;
end
if nargin < 7
    start_offset_in_points = 0;
end
if nargin < 9
    gain_t = [];
    gain_k = [];
end

f_carrier = 1575.42e6;
f_chip = 1.023e6 * (1 + initial_doppler/f_carrier);
N_chips = 1023;
N_sequences = 10;
N_fract = 10;   % number of fractional bits on the code phase used in ca_code_nco_v2
N_per_block = round(fs/f_chip * N_chips*N_sequences);

delta_doppler = 0.1 * f_chip/(N_chips*N_sequences);

N_oversampling = 1;

if nargin < 8 || isempty(stop_offset_in_points)
    if bInterleaved
        stop_offset_in_points = length(data_baseband)/2-2*N_per_block;
    else
        stop_offset_in_points = length(data_baseband)-2*N_per_block;
    end
end


% Loop tuning parameters: 1 means "réponse pile", time constant is roughly
% 1/K_delay_locked_loop * N_per_block/fs
% Ki_delay_locked_loop = 0.1 * N_sequences/2;
% this works with N_sequences = 10, 10 times oversampling
Ki_delay_locked_loop = 0.1 * 2/N_sequences * 1;
% trying to make it works with N_sequences = 10, no oversampling
Ki_delay_locked_loop = 0.1 * 2/N_sequences * 5/N_oversampling;
f_unity_gain = Ki_delay_locked_loop/(2*pi); % unitless, normalized at the 'sampling rate' of the feedback system, with sampling rate fs/N_per_block
f_ii_corner = f_unity_gain/3;
Kii_delay_locked_loop = 2*pi*Ki_delay_locked_loop*f_ii_corner;

Ki_pll = 20 * N_sequences/2;
Kii_pll = 10 * N_sequences/2;
% this works with N_sequences = 10, 10 times oversampling
Ki_pll = 20*1 * 1;
Kii_pll = 10*1 * 1; % threshold of stability is around 30 for N_sequences = 5
% trying to make it works with N_sequences = 10, no oversampling
Ki_pll = 20*1 * 1 * 0.5 * (10/10);
Kii_pll = 10*1 * 1 * 0.5 * (10/10); % threshold of stability is around 30 for N_sequences = 5
% works well except for launch, N_oversampling = 1
Ki_pll = 20*1 * 1 * 0.5 * (10/10);
Kii_pll = 10*1 * 1 * 0.5 * (10/10); % threshold of stability is around 30 for N_sequences = 5
Kiii_pll = 10*1 * 1 * 0.5 * 0.1; % threshold of stability is around 30 for N_sequences = 5
% trying to make it work with N_oversampling = 3
Ki_pll = 20*1 * 1 * 0.5 * (10/10);
Kii_pll = 10*1 * 1 * 0.5 * (10/10); % threshold of stability is around 30 for N_sequences = 5
Kiii_pll = 10*1 * 1 * 0.5 * 0.1; % threshold of stability is around 30 for N_sequences = 5

phase_error_integral = 0;
phase_error_integral2 = 0;

current_gain = 1;   % for time-varying loop gains

% we loop over sub-sequences of length N_per_block

% if bInterleaved
%     ar_start = 1+start_offset_in_points:N_per_block:stop_offset_in_points-N_per_block;  % this assumes 0% overlap between blocks, we could improve things a bit by adding some overlap
% else
    ar_start = 1+start_offset_in_points:round(N_per_block/N_oversampling):stop_offset_in_points-N_per_block;  % this assumes 0% overlap between blocks, we could improve things a bit by adding some overlap
% end


predicted_code_phase = initial_code_phase;
predicted_doppler = initial_doppler;
predicted_carrier_phase_at_start_of_block = 0;

kOut = 0;
timestamp = NaN*zeros(length(ar_start), 1);
code_phase = NaN*zeros(length(ar_start), 1);
carrier_phase_out = NaN*zeros(length(ar_start), 1);
doppler = NaN*zeros(length(ar_start), 1);
vector_prompt = NaN*zeros(length(ar_start), 1);
vector_prompt0 = NaN*zeros(length(ar_start), 1);
vector_promptp1 = NaN*zeros(length(ar_start), 1);
vector_promptm1 = NaN*zeros(length(ar_start), 1);
vector_early = NaN*zeros(length(ar_start), 1);
vector_late = NaN*zeros(length(ar_start), 1);
code_phase_error_integral = 0;

freq_error = NaN*zeros(length(ar_start), 1);
phase_error = NaN*zeros(length(ar_start), 1);



% tic;
% ca_code_phase = predicted_code_phase + (0:N_per_block-1) * f_chip/fs;
% ca_signal_prompt = ca_code_nco_v2(PRN_number, 2^N_fract*ca_code_phase );
% toc;
% 
% kStart = 1;
% data_small = data_baseband(ar_start(kStart) + (0:N_per_block-1));
% data_small = data_small .* exp(-1j*2*pi*predicted_doppler/fs*(0:N_per_block-1).' - 1j*predicted_carrier_phase);
% tic;
% data_small = interp1((0:N_per_block-1).', data_small, (0:N_per_block-1).', '*linear');
% toc;

tic;
for kStart = 1:length(ar_start)
    ca_code_phase = predicted_code_phase + (0:N_per_block-1) * f_chip/fs;
    ca_signal_prompt = ca_code_nco_v2(PRN_number, 2^N_fract*ca_code_phase, 'linear' );
    ca_signal_early = ca_code_nco_v2(PRN_number, 2^N_fract*(ca_code_phase+0.5), 'linear' );
    ca_signal_late = ca_code_nco_v2(PRN_number, 2^N_fract*(ca_code_phase-0.5), 'linear' );
    % slice out data chunk
    if bInterleaved
        data_small = complex(double(data_baseband(1+2*(ar_start(kStart)-1) + (0:2:2*N_per_block-1))), double(data_baseband(1+2*(ar_start(kStart)-1) + (1:2:2*N_per_block-1))));
        data_small = data_small - mean(data_small); % in case the data is offset binary and has DC
    else
        data_small = data_baseband(ar_start(kStart) + (0:N_per_block-1));
    end
    % remove doppler shift, with correct initial phase!
    predicted_carrier_phase = predicted_carrier_phase_at_start_of_block + 2*pi*predicted_doppler/fs*(0:N_per_block-1).';
    data_small = data_small .* exp(-1j*predicted_carrier_phase);
    % :
    
    % We can now compute the outputs of the loop
    % de-spread the chip signal, integrate and dump
    kOut = kOut + 1;
    timestamp(kOut) = ar_start(kStart) + (N_per_block-1)/2;
    
    % vary loop gain vs time during dataset to help tracking at very low
    % SNRs:
    if ~isempty(gain_t)
        current_gain = interp1(gain_t, gain_k, timestamp(kOut)/fs, 'linear', NaN);
        % change extrapolation behavior: we hold first or last gain if out
        % of bound
        if isnan(current_gain)
            if timestamp(kOut)/fs > gain_t(end)
                current_gain = gain_k(end);
            end
            if timestamp(kOut)/fs < gain_t(1)
                current_gain = gain_k(1);
            end
        end
        
    end

    
    

    
    % We try three different doppler shifts and choose the one that gives the largest magnitude output
    data_small_despread = data_small .* ca_signal_prompt;
    data_small_despreadp1 = data_small_despread .* exp(-1j*2*pi*delta_doppler/fs*(0:N_per_block-1).');
    data_small_despreadm1 = data_small_despread .* exp(1j*2*pi*delta_doppler/fs*(0:N_per_block-1).');
    vector_prompt0(kOut) = sum(data_small_despread);
    vector_promptp1(kOut) = sum(data_small_despreadp1);
    vector_promptm1(kOut) = sum(data_small_despreadm1);
%     vector_prompt(kOut) = vector_prompt0(kOut);
    % Select the doppler shift which gave the largest correlation
    if abs(vector_promptp1(kOut)) > abs(vector_prompt0(kOut)) && current_gain > 0
        % we had better correlation when removing an extra positive doppler shift
        vector_prompt(kOut) = vector_promptp1(kOut);
        data_small = data_small .* exp(-1j*2*pi*delta_doppler/fs*(0:N_per_block-1).');
        predicted_doppler = predicted_doppler + delta_doppler;
        initial_doppler = initial_doppler + delta_doppler;
    elseif abs(vector_promptm1(kOut)) > abs(vector_prompt0(kOut)) && current_gain > 0
        % we had better correlation when removing an extra negative doppler shift
        vector_prompt(kOut) = vector_promptm1(kOut);
        data_small = data_small .* exp(1j*2*pi*delta_doppler/fs*(0:N_per_block-1).');
        predicted_doppler = predicted_doppler - delta_doppler;
        initial_doppler = initial_doppler - delta_doppler;
    else
        % this is the normal case
        vector_prompt(kOut) = vector_prompt0(kOut);
%         data_small
    end
    
%     vector_prompt(kOut) = sum(data_small .* ca_signal_prompt);
    vector_early(kOut) = sum(data_small .* ca_signal_early);
    vector_late(kOut) = sum(data_small .* ca_signal_late);
    
    
    
    % Compute code phase error
    code_phase_error = (abs(vector_early(kOut)) - abs(vector_late(kOut)))./(abs(vector_early(kOut)) + abs(vector_late(kOut)));   % this is in units of chips
    if isnan(code_phase_error)
        code_phase_error = 0;
    end

    code_phase_error_integral = code_phase_error_integral + current_gain*code_phase_error;
    code_phase(kOut) = predicted_code_phase + (N_per_block-1)/2/fs*f_chip + current_gain*code_phase_error;

    
    % compute carrier phase error
    if kStart == 1
        mean_prompt_ampl = abs(vector_prompt(kOut));
    end
    
    if real(vector_prompt(kOut)) > 0
        % Bit was +1
%         phase_error(kOut) = angle(vector_prompt(kOut));
        phase_error(kOut) = imag(vector_prompt(kOut))./mean_prompt_ampl;
        
    else
        % Bit was -1
%         phase_error(kOut) = angle(-vector_prompt(kOut));
        phase_error(kOut) = -imag(vector_prompt(kOut))./mean_prompt_ampl;
    end
    carrier_phase_out(kOut) = predicted_carrier_phase_at_start_of_block + phase_error(kOut) + 2*pi*(N_per_block-1)*predicted_doppler/fs;
    doppler(kOut) = predicted_doppler;
    
    % amplitude estimate used to scale phase error signal
    mean_prompt_ampl = mean_prompt_ampl * 0.999 + abs(vector_prompt(kOut)) * (1-0.999);
    
    % update values for next iteration
    if kStart < length(ar_start)
        % update the code phase prediction for the next iteration
        predicted_code_phase = predicted_code_phase + (ar_start(kStart+1)-ar_start(kStart)) * f_chip/fs ...
                            + Ki_delay_locked_loop * code_phase_error ...
                            + Kii_delay_locked_loop * code_phase_error_integral;
                        
                        
        % update predicted_doppler, and predicted_phase_offset
        % (implement the BPSK-tolerant phase-locked loop)
        phase_error_integral = phase_error_integral + current_gain*phase_error(kOut);
        phase_error_integral2 = phase_error_integral2 + phase_error_integral;
        predicted_doppler = initial_doppler + Ki_pll*current_gain*phase_error(kOut) + Kii_pll * phase_error_integral + Kiii_pll * phase_error_integral2;
        predicted_carrier_phase_at_start_of_block = predicted_carrier_phase_at_start_of_block + 2*pi*(ar_start(kStart+1)-ar_start(kStart)) * predicted_doppler/fs;
        
%         predicted_carrier_phase_at_start_of_block = predicted_carrier_phase_at_start_of_block + (ar_start(kStart+1)-ar_start(kStart)) * f_chip/fs;
       
       
    end
    
    
    if mod(kStart, round(length(ar_start)/10)) == 0
        toc;
    end
end

end

function ca_code_tracking_v1_unit_tests()
    close all;
    % generate a fake signal to track to demonstrate the behavior of the
    % function
    fs = 2.048e6;
    f_carrier = 1.57542e9;
    doppler_thruth = 1e3;
    initial_code_phase_thruth = 0.1;
    PRN_test = 22;
    f_chip = 1.023e6;
    N_oversampling = fs/f_chip;
    N_sequences = 5*300;
    N_chips = 1023;
    N_fract = 10;
    
    bInterleaved = 0;
    
    %%
    f_doppler = 100+25;
    doppler_thruth = 100;
    code_phase_thruth = initial_code_phase_thruth + (0:1/N_oversampling*(1+f_doppler/f_carrier):N_chips*N_sequences).';
    
    fake_signal = ca_code_nco_v2(PRN_test, 2^N_fract*code_phase_thruth);
    
    fake_signal = fake_signal .* exp(1j*2*pi*f_doppler/fs*(0:length(code_phase_thruth)-1).' + 1j*0.1);
    
    [timestamp code_phase doppler vector_prompt vector_early vector_late carrier_phase_out phase_error] = ca_code_tracking_v1(fake_signal, bInterleaved, PRN_test, fs, doppler_thruth, initial_code_phase_thruth);
    
    figure;
    plot(phase_error);
    title('Phase error when starting with a doppler error of 25 Hz');
%     return
    
    %%
    

    
    %% Test #1: a small frequency error on the code phase
    f_doppler = 15e3; 
    doppler_thruth = 15e3;
    code_phase_thruth = initial_code_phase_thruth + (0:1/N_oversampling*(1+f_doppler/f_carrier):N_chips*N_sequences).';
    code_phase_reference = initial_code_phase_thruth + (0:1/N_oversampling*1.0000:N_chips*N_sequences).';
    fake_signal = ca_code_nco_v2(PRN_test, 2^N_fract*code_phase_thruth);

    fake_signal = fake_signal .* exp(1j*2*pi*f_doppler/fs*(0:length(fake_signal)-1).');
    
    [timestamp code_phase doppler vector_prompt vector_early vector_late carrier_phase_out] = ca_code_tracking_v1(fake_signal, bInterleaved, PRN_test, fs, doppler_thruth, initial_code_phase_thruth);
    
    figure;
    subplot(211);
    plot(timestamp/fs, abs(vector_prompt), '.-');
    hold all;
    plot(timestamp/fs, abs(vector_early), '.-');
    plot(timestamp/fs, abs(vector_late), '.-');
    legend('prompt', 'early', 'late');
    title(sprintf('Loop response to a velocity step of f_{doppler} = %f kHz', f_doppler/1e3));
    
    offset_error = 1/f_chip * (abs(vector_early) - abs(vector_late))./(abs(vector_early) + abs(vector_late));
    
    subplot(212);
%     plot((0:length(code_phase_thruth)-1)/fs, (code_phase_thruth-code_phase_reference(1:length(code_phase_thruth)))/f_chip);
    plot((0:length(code_phase_thruth)-1)/fs, (interp1(timestamp/fs, code_phase, (1:length(code_phase_thruth)).'/fs, 'linear')-code_phase_thruth)/f_chip);
    hold all;
    plot(timestamp/fs, offset_error, '.-');
    legend('ground thruth error', 'in-loop error');
    ylabel('Code phase error [s]');
    
    
    
    %% Test #2: step response of the delay locked loop
    time_step = 0.1/f_chip;
    doppler_thruth = 0;
    code_phase_thruth = initial_code_phase_thruth + (0:1/N_oversampling*1.0000:N_chips*N_sequences).';
    code_phase_thruth(round(length(code_phase_thruth)/2):end) = code_phase_thruth(round(length(code_phase_thruth)/2):end) + time_step*f_chip;
    code_phase_reference = initial_code_phase_thruth + (0:1/N_oversampling*1.0000:N_chips*N_sequences).';
    fake_signal = ca_code_nco_v2(PRN_test, 2^N_fract*code_phase_thruth);

    fake_signal = fake_signal .* exp(1j*2*pi*doppler_thruth/fs*(0:length(fake_signal)-1).');
    
    [timestamp code_phase doppler vector_prompt vector_early vector_late] = ca_code_tracking_v1(fake_signal, bInterleaved, PRN_test, fs, doppler_thruth, initial_code_phase_thruth);
    
    figure;
    subplot(211);
    plot(timestamp/fs, abs(vector_prompt), '.-');
    hold all;
    plot(timestamp/fs, abs(vector_early), '.-');
    plot(timestamp/fs, abs(vector_late), '.-');
    legend('prompt', 'early', 'late');
    title(sprintf('Loop response to time step of %s ns', time_step/1e-9));
    
    offset_error = 1/f_chip * (abs(vector_early) - abs(vector_late))./(abs(vector_early) + abs(vector_late));
    
    subplot(212);
    plot((0:length(code_phase_thruth)-1)/fs, (code_phase_thruth-code_phase_reference(1:length(code_phase_thruth)))/f_chip);
    hold all;
    plot(timestamp/fs, offset_error, '.-');
    
    %% Test #3: mimick actual dataset's velocity profile
    f_doppler = 15e3; f_carrier = 1.575e9;
    doppler_t = [0; 0.5; 0.5+3.5; 0.5+3.5+60; 0.5+3.5+160];
    doppler_f = [-225; -225; 2*1525; -225; -225];
    f_doppler = interp1(doppler_t, doppler_f, (0:3e6).'/fs, 'linear');
    doppler_thruth = f_doppler(1);

    code_phase_thruth = cumsum(  (1 + f_doppler/f_carrier) .* f_chip/fs  );
    initial_code_phase_thruth = code_phase_thruth(1);
    
%     length(code_phase_thruth)
%     code_phase_reference = initial_code_phase_thruth + (0:1/N_oversampling*1.0000:N_chips*N_sequences).';
    fake_signal = ca_code_nco_v2(PRN_test, 2^N_fract*code_phase_thruth);

    fake_signal = fake_signal .* exp(1j*2*pi*cumsum(f_doppler/fs));
    
    [timestamp code_phase doppler vector_prompt vector_early vector_late carrier_phase_out phase_error vector_promptp1 vector_promptm1 vector_prompt0] = ...
        ca_code_tracking_v1(fake_signal, bInterleaved, PRN_test, fs, doppler_thruth, initial_code_phase_thruth);
    
    figure;
    ax = subplot(311);
    plot(timestamp/fs, abs(vector_prompt), '.-');
    hold all;
    plot(timestamp/fs, abs(vector_early), '.-');
    plot(timestamp/fs, abs(vector_late), '.-');
    legend('prompt', 'early', 'late');
    title(sprintf('Loop response to a velocity profile similar to the rocket'));
    
    offset_error = 1/f_chip * (abs(vector_early) - abs(vector_late))./(abs(vector_early) + abs(vector_late));
    
    ax(2) = subplot(312);
%     plot((0:length(code_phase_thruth)-1)/fs, (code_phase_thruth-code_phase_reference(1:length(code_phase_thruth)))/f_chip);
    plot((0:10:length(code_phase_thruth)-1)/fs, (interp1(timestamp/fs, code_phase, (1:10:length(code_phase_thruth)).'/fs, 'linear')-code_phase_thruth(1:10:end))/f_chip);
    hold all;
    plot(timestamp/fs, offset_error, '.-');
    legend('ground thruth error', 'in-loop error');
    ylabel('Code phase error [s]');
    
    ax(3) = subplot(313);
    plot(doppler_t, doppler_f);
    ylabel('Doppler shift [Hz]');
    title('Velocity profile scaled as Doppler shift');
    xlabel('time [s]');
    
    linkaxes(ax, 'x');
    
    figure;
    plot(timestamp/fs, phase_error);
    title('Phase error');
    
    figure;
    plot(timestamp/fs, abs(vector_prompt0));
    hold all;
    plot(timestamp/fs, abs(vector_promptp1));
    hold all;
    plot(timestamp/fs, abs(vector_promptm1));
    plot(timestamp/fs, abs(vector_prompt));
    
%     keyboard
    
    fprintf(1, 'Todo: handle bit transitions better (right now we just assume that we are seeded correctly wrt to the transitions\n');
end

