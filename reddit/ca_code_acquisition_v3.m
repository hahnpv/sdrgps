% data_interleaved: complex data, centered at 0 frequency, real and imaginay parts interleaved, offset binary (uint8),
% zero is at +127, minimum at 0, maximum at 255
% fs: sampling rate, in Hz
% PRN_number_array: array of satellite vehicule PRN codes, possible values are from 1 to 37
% ar_start: offset in dataset at which to perform full doppler & code phase search
% ar_doppler_search: vector of doppler shifts to try, default is -max_doppler:f_chip/N_chip/4:max_doppler, where max_doppler is 10 kHz
% note that we return both an array of code epochs and one of code phase, but they are completely equivalent
% and are provided only for convenience. They are related by code_phase = code_freq * (t-code_epoch) =
% code_freq * t + code_phase, thus code_phase = -code_freq*code_epoch modulo N_pts_per_code_period, showing that code_phase is modulo N_pts_per_code_period, but code_epoch isn't
function [code_epoch_ar code_phase_ar doppler_ar correlation_amplitude_ar noise_amplitude_ar] = ca_code_acquisition_v3(data_interleaved, fs, PRN_number_array, ar_skip_pts, ar_doppler_search, N_sequences, bDisplay)
if nargin < 1
    fprintf(1, 'No argument passed. running unit tests...\n');
    ca_code_acquisition_v3_unittest();
    return
end

%% default parameters
if nargin < 3 || isempty(ar_skip_pts)
    ar_skip_pts = 0;
end
if nargin < 3
    ar_doppler_search = []; % we'll fill this up later in the function
end
if nargin < 5
    N_sequences = 1;
end
if nargin < 6
    bDisplay = 0;
end

%% physical parameters for the C/A code:
f_chip = 1.023e6;   % chip rate
N_chips = 1023; % number of chips in a sequence
% N_sequences = 10;   % number of sequence integrated coherently

N_pts_per_code_period = N_chips*fs/f_chip;
N_pts_per_correlation = round(2*N_pts_per_code_period*N_sequences);


%% Pre-compute PRN values:
tic;
code_phase_no_offset = (0:f_chip/fs:N_sequences*N_chips-f_chip/fs).';
PRN_sequences = NaN*zeros(length(code_phase_no_offset), length(PRN_number_array));

strLegend = {};
for kPRN = 1:length(PRN_number_array)
    ca_code_value = ca_code_nco_v2(PRN_number_array(kPRN), code_phase_no_offset*2^10, 'linear');
    PRN_sequences(:, kPRN) = flipud(ca_code_value);
    strLegend{kPRN} = sprintf('PRN %d', PRN_number_array(kPRN));
end

toc;

%% Doppler search space parameters
if isempty(ar_doppler_search)
    max_doppler = 10e3;     % max doppler shift search space in Hz, whole search space is +/- max_doppler
    doppler_delta = f_chip/(N_sequences*N_chips)/2/2;            % Doppler search increments in Hz
    ar_doppler_search = -max_doppler:doppler_delta:max_doppler;
    ar_doppler_search = -2.5e3:doppler_delta:2.5e3;
    fprintf(1, 'todo: implement pre-doppler cancellation to reduce dyn range of doppler\n');
end

% Pre-allocate output arrays
doppler_ar = NaN*zeros(length(ar_skip_pts), length(PRN_number_array));
correlation_amplitude_ar = NaN*zeros(length(ar_skip_pts), length(PRN_number_array));
code_epoch_ar = NaN*zeros(length(ar_skip_pts), length(PRN_number_array));
code_phase_ar = NaN*zeros(length(ar_skip_pts), length(PRN_number_array));
noise_amplitude_ar = NaN*zeros(length(ar_skip_pts), length(PRN_number_array));

%% Setup display figure
if bDisplay
    
    figure;
    ax = subplot(2, 1, 1);
    hLines1 = plot(ar_skip_pts/fs, doppler_ar, '.', 'markersize', 6);
    % xlabel('Time [s]');
    legend(strLegend);
    ylabel('Code epoch [s]');

    ax(2) = subplot(2, 1, 2);
    hLines2 = plot(ar_skip_pts/fs, doppler_ar, '.', 'markersize', 6);
    % title('Doppler');
    % xlabel('Time [s]');
    ylabel('Doppler [Hz]');
    xlabel('Time [s]');

    linkaxes(ax, 'x');
    set(gcf, 'position', [314.0000e+000   168.0000e+000   926.0000e+000   810.0000e+000]);
end

%% Actual computation loop
tic;
for kStart = 1:length(ar_skip_pts)
    for kPRN = 1:length(PRN_number_array)
        current_max = 0;
        current_code_offset = 0;
        current_doppler_ind = 0;
        for kDoppler = 1:length(ar_doppler_search)
            % slice out data and convert to complex, floating point numbers
            % for processing. we also perform conversion from offset binary
            % to zero-DC numbers here
            indices_small = 1 + 2*ar_skip_pts(kStart) + (0:2:2*N_pts_per_correlation-2).';
            data_small = complex(double(data_interleaved(indices_small)), double(data_interleaved(indices_small+1)));
            data_small = data_small - mean(data_small);
            % Remove Doppler shift and correlate data with PRN sequence
            data_shifted = data_small .* exp(-1j*2*pi*ar_doppler_search(kDoppler)/fs*(0:length(data_small)-1).');
            result = fftfilt((PRN_sequences(:, kPRN)), data_shifted);
            result = result(1+size(PRN_sequences, 1)-1:end);    % crop correlation startup
%             results_all(:, kDoppler) = result;  % save all results for debugging
            
            % compute a running max
            [valmax posmax] = max(abs(result));
            if valmax > current_max
                current_max = valmax;
                current_code_offset = posmax;
                current_doppler_ind = kDoppler;
            end
        end
        
        % save results for this PRN and starting offset:
        doppler_ar(kStart, kPRN) = ar_doppler_search(current_doppler_ind);
        correlation_amplitude_ar(kStart, kPRN) = current_max;
        code_epoch_ar(kStart, kPRN) = ar_skip_pts(kStart)+(current_code_offset);
        code_phase_ar(kStart, kPRN) = mod(-(current_code_offset-1) * f_chip/fs + N_chips/2, N_chips)-N_chips/2;
%         noise_amplitude_ar(kStart, kPRN) = noise_amplitude;


        % update display
        if bDisplay && (mod(kStart, 1) == 0)
            set(hLines1(kPRN), 'xdata', code_epoch_ar(:, kPRN)/fs, 'ydata', mod(code_epoch_ar(:, kPRN)/fs, N_chips/f_chip));
            set(hLines2(kPRN), 'xdata', code_epoch_ar(:, kPRN)/fs, 'ydata', doppler_ar(:, kPRN));
        %     set(hTime1, 'xdata', (1:length(data_small)), 'ydata', real(data_small));
        %     set(hTime2, 'xdata', (1:length(data_small)), 'ydata', abs(data_small));
            drawnow;
        end
    end
    toc;
end


end

function ca_code_acquisition_v3_unittest()
    close all
    
    %% parameters
    PRN_number_array = 22;  % for debugging
    % C/A code parameters
    f_chip = 1.023e6;   % chip rate
    fs = 10*f_chip; % sampling rate
    
    N_chips = 1023; % number of chips in a sequence
    N_sequences = 5;   % number of sequences integrated coherently
    N_pts_per_code_period = N_chips*fs/f_chip;
    N_pts_per_correlation = round(2*N_pts_per_code_period*N_sequences);

    %% inject a fake signal with a known doppler and code phase:
    % datar = 127 + 0*data(1:2:N_fake);
    % datai = 127 + 0*data(2:2:N_fake);
    doppler_ground_thruth = 1e3;
    code_phase_offset_ground_thruth = 50; % in units of chips
    code_phase_ground_thruth = code_phase_offset_ground_thruth + (0:f_chip/fs:4*N_sequences*N_chips-f_chip/fs).';   % we need at least 2 full sequences to do a correct search so we put a bit more

    ca_code_value = ca_code_nco_v2(PRN_number_array(1), code_phase_ground_thruth*2^10, 'linear');
    datar = 127+127*real(ca_code_value .* exp(1j*2*pi*doppler_ground_thruth/fs*(1:length(ca_code_value)).'));
    datai = 127+127*imag(ca_code_value .* exp(1j*2*pi*doppler_ground_thruth/fs*(1:length(ca_code_value)).'));
    % have to play games to match the format of the data files a bit
    data = [datar datai].';
    data = data(:);


    %% Test #1: code phase search only, no doppler, hardcoded code epoch ground thruth:
    ar_skip_pts = 0;
    ar_doppler_search = 1e3;

    bDisplay = 0;
    [code_epoch_ar code_phase_ar doppler_ar correlation_amplitude_ar noise_amplitude_ar] = ...
        ca_code_acquisition_v3(data, fs, PRN_number_array, ar_skip_pts, ar_doppler_search, N_sequences, bDisplay);
    if round(mod(code_epoch_ar+N_pts_per_code_period/2, N_pts_per_code_period)-N_pts_per_code_period/2) == -499 & round(100*code_phase_ar) == round(100*code_phase_ground_thruth(1)) & round(1*doppler_ar) == round(1*doppler_ground_thruth)
        fprintf(1, 'PASSED test #1\n');
    else
        fprintf(1, 'FAILED test #1\n');
        keyboard
    end
    
    figure;
    plot(ar_skip_pts+1, code_phase_ar, '.');
    hold all;
    plot((1:length(code_phase_ground_thruth)), mod(code_phase_ground_thruth, N_chips));
    plot(reshape([repmat(code_epoch_ar, 1, 2) NaN*zeros(length(code_epoch_ar), 1)].', [], 1), repmat([0 N_chips*1.1 NaN], length(code_epoch_ar)));
    

    
    %% Test #2: code phase and doppler search, hardcoded code epoch ground thruth:
    ar_skip_pts = 0;
    ar_doppler_search = [];

    bDisplay = 0;
    [code_epoch_ar code_phase_ar doppler_ar correlation_amplitude_ar noise_amplitude_ar] = ...
        ca_code_acquisition_v3(data, fs, PRN_number_array, ar_skip_pts, ar_doppler_search, N_sequences, bDisplay);
    if round(mod(code_epoch_ar+N_pts_per_code_period/2, N_pts_per_code_period)-N_pts_per_code_period/2) == -499 & round(100*code_phase_ar) == round(100*code_phase_ground_thruth(1)) & round(1*doppler_ar) == round(1*doppler_ground_thruth)
        fprintf(1, 'PASSED test #2\n');
    else
        fprintf(1, 'FAILED test #2\n');
    end
    

    

    %% Test #3: test for code epoch accuracy by searching for codes at known epochs:
    code_epoch_ground_thruth1 = 1e3;
    code_epoch_ground_thruth2 = 1.218e5;
    code_phase_ground_thruth = (0:f_chip/fs:4*N_sequences*N_chips-f_chip/fs).';   % we need at least 2 full sequences to do a correct search so we put a bit more
    ca_code_value = ca_code_nco_v2(PRN_number_array(1), code_phase_ground_thruth*2^10, 'linear');
    datar = 127*ones(code_epoch_ground_thruth2 + N_pts_per_correlation, 1);
    datai = 127*ones(code_epoch_ground_thruth2 + N_pts_per_correlation, 1);
    datar(code_epoch_ground_thruth1 + (0:N_pts_per_correlation/2-1)) = 127+127*real(ca_code_value(1:N_pts_per_correlation/2));
    datar(code_epoch_ground_thruth2 + (0:N_pts_per_correlation/2-1)) = 127+127*real(ca_code_value(1:N_pts_per_correlation/2));
%     datai(code_epoch_ground_thruth + (0:N_pts_per_code_period-1)) = 127+127*imag(ca_code_value(code_epoch_ground_thruth + (0:N_pts_per_code_period-1)));
    % have to play games to match the format of the data files a bit
    data = [datar datai].';
    data = data(:);
    
    ar_skip_pts = [100; 200; 950; 1.103e5];
    ar_doppler_search = 0;

    bDisplay = 0;
    [code_epoch_ar code_phase_ar doppler_ar correlation_amplitude_ar noise_amplitude_ar] = ...
        ca_code_acquisition_v3(data, fs, PRN_number_array, ar_skip_pts, ar_doppler_search, N_sequences, bDisplay);
    figure;
    plot(datar);
    hold all;
%     plot((1:length(code_phase_ground_thruth)), mod(code_phase_ground_thruth, N_chips));
    plot(reshape([repmat(code_epoch_ar, 1, 2) NaN*zeros(length(code_epoch_ar), 1)].', [], 1), repmat([0 N_chips*1.1 NaN], length(code_epoch_ar)));
    if all(code_epoch_ar(1:3) == code_epoch_ground_thruth1) && all(code_epoch_ar(4:end) == code_epoch_ground_thruth2)
        fprintf(1, 'PASSED test #3\n');
    else
        fprintf(1, 'FAILED test #3\n');
        keyboard
        
    end
    


    
    %% Test #4: code phase search only (for speedup), run a few times to check phase at a few starting offsets:
    code_phase_offset_ground_thruth = 50; % in units of chips
    code_phase_ground_thruth = code_phase_offset_ground_thruth + (0:f_chip/fs:4*N_sequences*N_chips-f_chip/fs).';   % we need at least 2 full sequences to do a correct search so we put a bit more
    ca_code_value = ca_code_nco_v2(PRN_number_array(1), code_phase_ground_thruth*2^10, 'linear');
    datar = 127+127*real(ca_code_value .* exp(1j*2*pi*doppler_ground_thruth/fs*(1:length(ca_code_value)).'));
    datai = 127+127*imag(ca_code_value .* exp(1j*2*pi*doppler_ground_thruth/fs*(1:length(ca_code_value)).'));
    % have to play games to match the format of the data files a bit
    data = [datar datai].';
    data = data(:);
    
    ar_skip_pts = [100; 1e3; N_pts_per_correlation];
    ar_doppler_search = 1e3;

    bDisplay = 0;
    [code_epoch_ar code_phase_ar doppler_ar correlation_amplitude_ar noise_amplitude_ar] = ...
        ca_code_acquisition_v3(data, fs, PRN_number_array, ar_skip_pts, ar_doppler_search, N_sequences, bDisplay);
    figure;
    plot(ar_skip_pts+1, code_phase_ar, '.');
    hold all;
    plot((1:length(code_phase_ground_thruth)), mod(code_phase_ground_thruth, N_chips));
    plot(reshape([repmat(code_epoch_ar, 1, 2) NaN*zeros(length(code_epoch_ar), 1)].', [], 1), repmat([0 N_chips*1.1 NaN], length(code_epoch_ar)));
    if round(100*code_phase_ar) == round(100*mod(interp1((1:length(code_phase_ground_thruth)).', code_phase_ground_thruth, ar_skip_pts+1, 'linear'),N_chips)) & ...
            round(1*doppler_ar) == round(1*doppler_ground_thruth)
        fprintf(1, 'PASSED test #4\n');
    else
        fprintf(1, 'FAILED test #4\n');
        keyboard
        
    end

end