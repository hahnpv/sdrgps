

% load and display RINEX observation data
strFile = 'd:\Downloads\GSDR086x11.16O';
f = fopen(strFile, 'r');
% skip the first header lines:
for kSkip = 1:16
    fgetl(f);
end

kOut = 0;
timestamps = [];
obs_code_phase = NaN*zeros(37, 1);
obs_doppler = NaN*zeros(37, 1);

while ~feof(f)
    strLine = fgetl(f);
%     disp(strLine);
    if strLine(1) == '>'
        % this is a timestamp line
        timestamp_vector = sscanf(strLine(2:end), '%f');
        kOut = kOut + 1;
        timestamps(kOut) = timestamp_vector(5)*60+timestamp_vector(6);
        obs_code_phase(:, kOut) = NaN;
        obs_doppler(:, kOut) = NaN;
        continue;
    end
    
    if strLine(1) == 'G'
        % this is an observation line
        obs_vector = sscanf(strLine(2:end), '%f');
%         sat_num = str2double(strLine(2:3))
        sat_num = obs_vector(1);
%         disp(obs_vector);
        obs_code_phase(sat_num, kOut) = obs_vector(2);
        obs_doppler(sat_num, kOut) = obs_vector(6);
%         obs_code_phase(sat_num, kOut) = 1;
    end
end

fclose(f);

%% remove sats that were never observed in the dataset:
ar_PRN_logical = true(size(obs_code_phase, 1), 1);
for k = 1:size(obs_code_phase, 1)
    if all(isnan(obs_code_phase(k, :)))
        ar_PRN_logical(k) = false;
    end
end

obs_code_phase = obs_code_phase(ar_PRN_logical, :);
obs_doppler = obs_doppler(ar_PRN_logical, :);

strLegendRINEX = {};
ar_PRN_ind = find(ar_PRN_logical);
for k = 1:sum(ar_PRN_logical)
    strLegendRINEX{k} = sprintf('PRN %d', ar_PRN_ind(k));
end

%% display results

c = 299792458;

figure;
ax = subplot(211);
plot(timestamps-timestamps(1), (obs_code_phase-min(min(obs_code_phase)))/c, '.');
ylabel('Pseudorange - smallest pseudorange [s]');
legend(strLegendRINEX);
ax(2) = subplot(212);
plot(timestamps-timestamps(1), obs_doppler, '.');
ylabel('Doppler shift [Hz]');

% first_code_phase_for_each_sat = [];
% for k = 1:size(obs_code_phase, 1)
%     
%     first_ind_for_each_sat = find(~isnan(obs_code_phase(k, :)), 1);
%     first_code_phase_for_each_sat(k) = obs_code_phase(k, first_ind_for_each_sat);
% end

% figure;
% % ax = subplot(211);
% plot(timestamps-timestamps(1), (obs_code_phase-repmat(first_code_phase_for_each_sat.', 1, size(obs_code_phase, 2)))/c, '.-');
% ylabel('Pseudorange offset [s]');
% legend(strLegendRINEX);
% % ax(2) = subplot(212);
% % plot(timestamps-timestamps(1), obs_doppler, '.-');
% % ylabel('Doppler shift [Hz]');


