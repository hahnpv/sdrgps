clear
close all


%% 1. Load data files (place data as a global variable so it stays even when "clear" is called
N_max_bytes = inf; % should be an even number
szFile = 'C:\Users\phahn\Downloads\Feb6.bin';

global data_global    % declared as a global so that it stays in memory even when using "clear" ("clear all" would clear it)
if isempty(data_global)
    fprintf(1, 'Loading dataset...'); tic;
    
    f = fopen(szFile, 'rb');
    data_global = fread(f, N_max_bytes, 'uint8=>uint8');    % format is offset binary
    fclose(f);
    
    fprintf(1, 'Done (%.0f sec)\n', toc);
end

data = data_global; % note that this doesn't actually make a copy unless "data" is modified
N_pts_total = length(data)/2;



%% 2. Run the coarse/acquisition search

PRN_number_array = 22;  % for debugging, the strongest one
PRN_number_array = [1 3 10 14 22 31];  % PRNs found by GNSS-SDR
% PRN_number_array = [1 3 8 10 11 14 22 23 25 31];  % those gave reasonable tracks
% PRN_number_array = [1 3 10 11 14 22 25 31];  % those tracked through the launch
% PRN_number_array = (1:37);    % exhaustive search

% C/A code parameters
fs = 2048000;   % sampling rate
f_chip = 1.023e6;   % chip rate
N_chips = 1023; % number of chips in a sequence
N_sequences = 10;   % number of sequence integrated coherently
N_pts_per_code_period = N_chips*fs/f_chip;
N_pts_per_correlation = round(2*N_pts_per_code_period*N_sequences);


% Offsets in file at which to perform coarse search (technically, only one
% is needed, but its useful to take more since it requires no tracking loop to work)
ar_skip_pts = round(linspace(1e6, round(180*fs), 15).');
% ar_skip_pts = 1e6;


% Doppler search space parameters
max_doppler = 10e3;     % max doppler shift search space in Hz, whole search space is +/- max_doppler
doppler_delta = f_chip/(N_sequences*N_chips)/2/2;            % Doppler search increments in Hz
ar_doppler_search = -max_doppler:doppler_delta:max_doppler;
ar_doppler_search = -4e3:doppler_delta:4e3;
fprintf(1, 'todo: implement pre-doppler cancellation to reduce dyn range of doppler\n');



% Run the coarse/acquisition search:
bDisplay = 1;
[code_epoch_ar code_phase_ar doppler_ar correlation_amplitude_ar noise_amplitude_ar] = ...
    ca_code_acquisition_v3(data, fs, PRN_number_array, ar_skip_pts, ar_doppler_search, N_sequences, bDisplay);


fprintf(1, 'coarse/acquisition search computation speed: %f pts/sec\n', numel(code_epoch_ar)/toc);
fprintf(1, 'coarse/acquisition search computation time: %f sec/pts\n', toc/numel(code_epoch_ar));


%% 3. find location of bit transitions w.r.t. to code, 2nd try, this time with longer integration
% need to run for more than 20 chips, since there might be no bits
% transitions in the data stream
%
% approach is: integrate and dump by integer number of sequences
% with no tracking. then compute x(k).*conj(x(k-1)), and figure out where
% the bits flipped
N_sequences_per_bit = 20;
N_bits = 25;   % this should be high enough to make sure we catch a transition (I don't know if there is even a guarantee for this in the gps data stream... probably would need to be more than 1 word long at the very least)
N_fract = 10;
kPRN = 5;
N_skip_sequences_to_start_right_after_bit_transitions = 0*PRN_number_array;
for kPRN = 1:length(PRN_number_array)

    ca_code_phase = (0:f_chip/fs:N_chips*N_sequences_per_bit*N_bits-f_chip/fs).';
    N_pts_for_sync = length(ca_code_phase);

    ca_code = ca_code_nco_v2(PRN_number_array(kPRN), 2^N_fract*ca_code_phase);
    start_offset = code_epoch_ar(1, kPRN)-1;
    data_small = complex(double(data(1+2*(start_offset-1) + (0:2:2*N_pts_for_sync-1))), double(data(1+2*(start_offset-1) + (1:2:2*N_pts_for_sync-1))));
    data_small = data_small - mean(data_small); % in case the data is offset binary and has DC
    data_small = data_small .* ca_code;
    data_small = data_small .* exp(-1j*2*pi*doppler_ar(1, kPRN)/fs*(0:length(data_small)-1).');
    vector_prompt = reshape(data_small, round(N_chips*fs/f_chip), []);
    vector_prompt = sum(vector_prompt, 1);
    


    metric = 0;
    for k = 1:20
        % average into groups of 20 sequences for each bit, starting at sequence k
        vector_prompt_integrated = vector_prompt(k:end);
        vector_prompt_integrated = vector_prompt_integrated(1:floor(length(vector_prompt_integrated)/N_sequences_per_bit)*N_sequences_per_bit);  % ensure we have an integer number of bits for reshape() to work
        vector_prompt_integrated = sum(reshape(vector_prompt_integrated, N_sequences_per_bit, []), 1);
        vector_promptd = vector_prompt_integrated(1:end-1) .* conj(vector_prompt_integrated(1+1:end));
%         vector_promptd = vector_promptd(1:floor(length(vector_promptd)/N_sequences_per_bit)*N_sequences_per_bit);  % ensure we have an integer number of bits for reshape() to work
        metric(k) = mean(   abs(vector_promptd)  );
    end
    [valmmax posmax] = max(real(metric));

    N_skip_sequences_to_start_right_after_bit_transitions(kPRN) = mod(posmax-1, N_sequences_per_bit);

    figure;
    subplot(311);
    plot(abs(vector_prompt));
    title('Vector prompt amplitude vs chip number');
    subplot(312);
    plot(real(metric), '.-');
    hold all;
    plot(imag(metric), '.-');
    title(sprintf('Metric for bit transitions. PRN = %d', PRN_number_array(kPRN)));
    
    % for debugging: show that the code transitions are well aligned:
    ca_code_phase = (0:f_chip/fs:N_chips*N_sequences_per_bit*N_bits-f_chip/fs).';
    N_pts_for_sync = length(ca_code_phase);

    ca_code = ca_code_nco_v2(PRN_number_array(kPRN), 2^N_fract*ca_code_phase);
    start_offset = code_epoch_ar(1, kPRN)-1 + N_skip_sequences_to_start_right_after_bit_transitions(kPRN) * N_pts_per_code_period;
    data_small = complex(double(data(1+2*(start_offset-1) + (0:2:2*N_pts_for_sync-1))), double(data(1+2*(start_offset-1) + (1:2:2*N_pts_for_sync-1))));
    data_small = data_small - mean(data_small); % in case the data is offset binary and has DC
    data_small = data_small .* ca_code;
    data_small = data_small .* exp(-1j*2*pi*doppler_ar(1, kPRN)/fs*(0:length(data_small)-1).');
    vector_prompt = reshape(data_small, round(N_chips*fs/f_chip), []);
    vector_prompt = sum(vector_prompt, 1);
    
    subplot(313);
    plot(real(vector_prompt));
    hold all;
    plot(imag(vector_prompt));
    title(sprintf('Correlator outputs. PRN = %d', PRN_number_array(kPRN)));

end



%% run the code tracking
N_max = round(180*fs);
% N_max = round(50*fs);
bInterleaved = 1;   % data format is real, imag interleaved

% prepare plot to display results:
strLegend = {};
for kPRN = 1:length(PRN_number_array)
    strLegend{kPRN} = sprintf('PRN %d', PRN_number_array(kPRN));
end

color_order = get(0, 'DefaultAxesColorOrder');
figure;
ax = subplot(211);
plot(code_epoch_ar/fs, mod(-code_epoch_ar/fs, N_chips/f_chip), '.');
ylabel('Code phase - expected (modulo) [s]');
legend(strLegend);
ax(2) = subplot(212);
plot(code_epoch_ar/fs, doppler_ar, '.');
ylabel('Doppler [Hz]');


figure;
ax(3) = axes();
linkaxes(ax, 'x');
title('Correlator outputs');
xlabel('Time [s]');


% matrices to hold output results
timestamp_all = [];
code_phase_all = [];
doppler_all = [];

% Loop over each PRN number
for kPRN = 1:length(PRN_number_array)
    initial_doppler = doppler_ar(1, kPRN);
    PRN_number = PRN_number_array(kPRN);
    initial_code_phase = 0;

    % start the code tracking at the start of a code and also at the start
    % of a bit boundary.  this enables using integration lengths up to 20
    % chips (although we currently use 10)
    start_offset_in_points = round(code_epoch_ar(1, kPRN)-1  + N_skip_sequences_to_start_right_after_bit_transitions(kPRN) * N_pts_per_code_period );
    stop_offset_in_points = start_offset_in_points + N_max;
%     stop_offset_in_points = []; % this makes the function process the whole file

    % this is the function which does all the heavy lifting of code
    % tracking
    [timestamp code_phase doppler vector_prompt vector_early vector_late carrier_phase_out phase_error] = ca_code_tracking_v1(data, bInterleaved, PRN_number, fs, ...
                initial_doppler, initial_code_phase, start_offset_in_points, stop_offset_in_points);

    % save results in a matrix:
    timestamp_all(1:length(timestamp), kPRN) = timestamp;
    code_phase_all(1:length(code_phase), kPRN) = code_phase;
    doppler_all(1:length(doppler), kPRN) = doppler;

    % update display
    axes(ax(1));
    hold all;
    plot(timestamp/fs, mod(code_phase - timestamp*f_chip/fs, N_chips)/f_chip, 'color', color_order(mod(kPRN-1, length(color_order))+1, :));
    axes(ax(2));
    hold all;
    plot(timestamp/fs, doppler, 'color', color_order(mod(kPRN-1, length(color_order))+1, :));
    axes(ax(3));
    hold all;
    if kPRN == 1
        plot(timestamp/fs, real(vector_prompt), 'color', color_order(mod(kPRN-1, length(color_order))+1, :));
        last_max = max(real(vector_prompt));
    else
        plot(timestamp/fs, real(vector_prompt) + 1.1*max(real(vector_prompt)) + last_max, 'color', color_order(mod(kPRN-1, length(color_order))+1, :));
        last_max = max(real(vector_prompt)) + 1.1*max(real(vector_prompt)) + last_max;
    end
    legend(strLegend);
    
    drawnow;
end


%% one way of looking at the code phases vs time

% code_phase_all_minus_expected = (code_phase_all - timestamp_all*f_chip/fs);
% code_phase_all_minus_expected_common_offset = code_phase_all_minus_expected - repmat(code_phase_all_minus_expected(1, :), size(code_phase_all_minus_expected, 1), 1);
% 
% 
% figure;
% plot(timestamp_all/fs, -(code_phase_all_minus_expected_common_offset)/f_chip);


%% interpolate all phases to the same grid
% and generate pseudoranges estimates

% generate a common time grid
start_timestamp = max(min(timestamp_all, [], 1));
stop_timestamp = min(max(timestamp_all, [], 1));
common_timestamp = (start_timestamp:mean(mean(diff(timestamp_all))):stop_timestamp).';

% resample the phases on the same grid
code_phase_allr = [];
doppler_allr = [];
for kPRN = 1:size(timestamp_all, 2)
    code_phase_allr(:, kPRN) = interp1(timestamp_all(:, kPRN), code_phase_all(:, kPRN), common_timestamp, 'linear');
    doppler_allr(:, kPRN) = interp1(timestamp_all(:, kPRN), doppler_all(:, kPRN), common_timestamp, 'linear');
end

% subtract one common phase (kPRN = 5, PRN = 22) (because that's what
% GNSS-SDR did)
% dirty hack: we know that PRN 22 (kPRN = 5) was the closest sat in this dataset
kPRN_closest = 5;
figure;
plot(common_timestamp/fs, mod(-(code_phase_allr - repmat(code_phase_allr(:, kPRN_closest), 1, size(code_phase_allr, 2)))/f_chip, 20e-3)   );
title('Pseudoranges, subtracting all of PRN #22');
legend(strLegend);


code_phase_allr_aligned = -(code_phase_allr - repmat(common_timestamp/fs*f_chip, 1, size(code_phase_allr, 2)))/f_chip;
code_phase_allr_aligned = mod(code_phase_allr_aligned - min(code_phase_allr_aligned(:, kPRN_closest), [], 1)   , 20e-3);

figure;
plot(common_timestamp/fs,  code_phase_allr_aligned  );
title('Pseudoranges, subtracting only the smallest value');
legend(strLegend);


%%
% return

%% Overlay the observation data from the RINEX file and the observations we have just produced.

load_rinex_observation;

time_offset_rinex_file = 55.15; % the file seems to start at offset 55.15 sec into the dataset
doppler_offset_rinex = 115;     % the GNSS-computed Doppler shifts seem to have a 115 Hz offset (but I have no idea why)

%
color_order = get(0, 'DefaultAxesColorOrder');
for kPRN = 1:size(code_phase_allr_aligned, 2)
    axes(ax(1));
    hold on;
%     plot(common_timestamp/fs - time_offset_rinex_file,  code_phase_allr_aligned(:, kPRN), 'color', color_order(1+mod(kPRN-1, size(color_order, 1)), :)  );
    plot(common_timestamp/fs - time_offset_rinex_file,  mod(-(code_phase_allr(:, kPRN) - code_phase_allr(:, kPRN_closest))/f_chip, 20e-3), 'color', color_order(1+mod(kPRN-1, size(color_order, 1)), :)  );
    axes(ax(2));
    hold on;
    plot(common_timestamp/fs - time_offset_rinex_file,  doppler_allr(:, kPRN) - doppler_offset_rinex, 'color', color_order(1+mod(kPRN-1, size(color_order, 1)), :)  );
end
axes(ax(1));
legend(strLegendRINEX);
axes(ax(2));
legend(strLegend);

