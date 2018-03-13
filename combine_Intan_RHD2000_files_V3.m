function combine_Intan_RHD2000_files_V3

% combine_Intan_RHD2000_files
% Edited from read_Intan_RHD2000_file.m by Stewart Holloway
%
% Version 1.1 2/12/16
%
% A MATLAB function to combine RHD files.
%


num_files = input('Number of Files: ');


[file, path, filterindex] = ...
    uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'off');

% Read most recent file automatically.
%path = 'C:\Users\Reid\Documents\RHD2132\testing\';
%d = dir([path '*.rhd']);
%file = d(end).name;

tic;
filename = [path,file];
fid = fopen(filename, 'r+');

s = dir(filename);
filesize = s.bytes;


% Check 'magic number' at beginning of file to make sure this is an Intan
% Technologies RHD2000 data file.
magic_number = fread(fid, 1, 'uint32');
if magic_number ~= hex2dec('c6912702')
    error('Unrecognized file type.');
end

% Read version number.
data_file_main_version_number = fread(fid, 1, 'int16');
data_file_secondary_version_number = fread(fid, 1, 'int16');

fprintf(1, '\n');
fprintf(1, 'Reading Intan Technologies RHD2000 Data File, Version %d.%d\n', ...
    data_file_main_version_number, data_file_secondary_version_number);
fprintf(1, '\n');

% Read information of sampling rate and amplifier frequency settings.
sample_rate = fread(fid, 1, 'single');
dsp_enabled = fread(fid, 1, 'int16');
actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
actual_lower_bandwidth = fread(fid, 1, 'single');
actual_upper_bandwidth = fread(fid, 1, 'single');

desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
desired_lower_bandwidth = fread(fid, 1, 'single');
desired_upper_bandwidth = fread(fid, 1, 'single');

% This tells us if a software 50/60 Hz notch filter was enabled during
% the data acquisition.
notch_filter_mode = fread(fid, 1, 'int16');
notch_filter_frequency = 0;
if (notch_filter_mode == 1)
    notch_filter_frequency = 50;
elseif (notch_filter_mode == 2)
    notch_filter_frequency = 60;
end

desired_impedance_test_frequency = fread(fid, 1, 'single');
actual_impedance_test_frequency = fread(fid, 1, 'single');

% Place notes in data strucure
notes = struct( ...
    'note1', fread_QString(fid), ...
    'note2', fread_QString(fid), ...
    'note3', fread_QString(fid) );
    
% If data file is from GUI v1.1 or later, see if temperature sensor data
% was saved.
num_temp_sensor_channels = 0;
if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 1) ...
    || (data_file_main_version_number > 1))
    num_temp_sensor_channels = fread(fid, 1, 'int16');
end

% If data file is from GUI v1.3 or later, load eval board mode.
eval_board_mode = 0;
if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 3) ...
    || (data_file_main_version_number > 1))
    eval_board_mode = fread(fid, 1, 'int16');
end

% Place frequency-related information in data structure.
frequency_parameters = struct( ...
    'amplifier_sample_rate', sample_rate, ...
    'aux_input_sample_rate', sample_rate / 4, ...
    'supply_voltage_sample_rate', sample_rate / 60, ...
    'board_adc_sample_rate', sample_rate, ...
    'board_dig_in_sample_rate', sample_rate, ...
    'desired_dsp_cutoff_frequency', desired_dsp_cutoff_frequency, ...
    'actual_dsp_cutoff_frequency', actual_dsp_cutoff_frequency, ...
    'dsp_enabled', dsp_enabled, ...
    'desired_lower_bandwidth', desired_lower_bandwidth, ...
    'actual_lower_bandwidth', actual_lower_bandwidth, ...
    'desired_upper_bandwidth', desired_upper_bandwidth, ...
    'actual_upper_bandwidth', actual_upper_bandwidth, ...
    'notch_filter_frequency', notch_filter_frequency, ...
    'desired_impedance_test_frequency', desired_impedance_test_frequency, ...
    'actual_impedance_test_frequency', actual_impedance_test_frequency );

% Define data structure for spike trigger settings.
spike_trigger_struct = struct( ...
    'voltage_trigger_mode', {}, ...
    'voltage_threshold', {}, ...
    'digital_trigger_channel', {}, ...
    'digital_edge_polarity', {} );

new_trigger_channel = struct(spike_trigger_struct);
spike_triggers = struct(spike_trigger_struct);

% Define data structure for data channels.
channel_struct = struct( ...
    'native_channel_name', {}, ...
    'custom_channel_name', {}, ...
    'native_order', {}, ...
    'custom_order', {}, ...
    'board_stream', {}, ...
    'chip_channel', {}, ...
    'port_name', {}, ...
    'port_prefix', {}, ...
    'port_number', {}, ...
    'electrode_impedance_magnitude', {}, ...
    'electrode_impedance_phase', {} );

new_channel = struct(channel_struct);

% Create structure arrays for each type of data channel.
amplifier_channels = struct(channel_struct);
aux_input_channels = struct(channel_struct);
supply_voltage_channels = struct(channel_struct);
board_adc_channels = struct(channel_struct);
board_dig_in_channels = struct(channel_struct);
board_dig_out_channels = struct(channel_struct);

amplifier_index = 1;
aux_input_index = 1;
supply_voltage_index = 1;
board_adc_index = 1;
board_dig_in_index = 1;
board_dig_out_index = 1;

% Read signal summary from data file header.

number_of_signal_groups = fread(fid, 1, 'int16');

for signal_group = 1:number_of_signal_groups
    signal_group_name = fread_QString(fid);
    signal_group_prefix = fread_QString(fid);
    signal_group_enabled = fread(fid, 1, 'int16');
    signal_group_num_channels = fread(fid, 1, 'int16');
    signal_group_num_amp_channels = fread(fid, 1, 'int16');

    if (signal_group_num_channels > 0 && signal_group_enabled > 0)
        new_channel(1).port_name = signal_group_name;
        new_channel(1).port_prefix = signal_group_prefix;
        new_channel(1).port_number = signal_group;
        for signal_channel = 1:signal_group_num_channels
            new_channel(1).native_channel_name = fread_QString(fid);
            new_channel(1).custom_channel_name = fread_QString(fid);
            new_channel(1).native_order = fread(fid, 1, 'int16');
            new_channel(1).custom_order = fread(fid, 1, 'int16');
            signal_type = fread(fid, 1, 'int16');
            channel_enabled = fread(fid, 1, 'int16');
            new_channel(1).chip_channel = fread(fid, 1, 'int16');
            new_channel(1).board_stream = fread(fid, 1, 'int16');
            new_trigger_channel(1).voltage_trigger_mode = fread(fid, 1, 'int16');
            new_trigger_channel(1).voltage_threshold = fread(fid, 1, 'int16');
            new_trigger_channel(1).digital_trigger_channel = fread(fid, 1, 'int16');
            new_trigger_channel(1).digital_edge_polarity = fread(fid, 1, 'int16');
            new_channel(1).electrode_impedance_magnitude = fread(fid, 1, 'single');
            new_channel(1).electrode_impedance_phase = fread(fid, 1, 'single');
            
            if (channel_enabled)
                switch (signal_type)
                    case 0
                        amplifier_channels(amplifier_index) = new_channel;
                        spike_triggers(amplifier_index) = new_trigger_channel;
                        amplifier_index = amplifier_index + 1;
                    case 1
                        aux_input_channels(aux_input_index) = new_channel;
                        aux_input_index = aux_input_index + 1;
                    case 2
                        supply_voltage_channels(supply_voltage_index) = new_channel;
                        supply_voltage_index = supply_voltage_index + 1;
                    case 3
                        board_adc_channels(board_adc_index) = new_channel;
                        board_adc_index = board_adc_index + 1;
                    case 4
                        board_dig_in_channels(board_dig_in_index) = new_channel;
                        board_dig_in_index = board_dig_in_index + 1;
                    case 5
                        board_dig_out_channels(board_dig_out_index) = new_channel;
                        board_dig_out_index = board_dig_out_index + 1;
                    otherwise
                        error('Unknown channel type');
                end
            end
            
        end
    end
end

% Summarize contents of data file.
num_amplifier_channels = amplifier_index - 1;
num_aux_input_channels = aux_input_index - 1;
num_supply_voltage_channels = supply_voltage_index - 1;
num_board_adc_channels = board_adc_index - 1;
num_board_dig_in_channels = board_dig_in_index - 1;
num_board_dig_out_channels = board_dig_out_index - 1;

fprintf(1, 'Found %d amplifier channel%s.\n', ...
    num_amplifier_channels, plural(num_amplifier_channels));
fprintf(1, 'Found %d auxiliary input channel%s.\n', ...
    num_aux_input_channels, plural(num_aux_input_channels));
fprintf(1, 'Found %d supply voltage channel%s.\n', ...
    num_supply_voltage_channels, plural(num_supply_voltage_channels));
fprintf(1, 'Found %d board ADC channel%s.\n', ...
    num_board_adc_channels, plural(num_board_adc_channels));
fprintf(1, 'Found %d board digital input channel%s.\n', ...
    num_board_dig_in_channels, plural(num_board_dig_in_channels));
fprintf(1, 'Found %d board digital output channel%s.\n', ...
    num_board_dig_out_channels, plural(num_board_dig_out_channels));
fprintf(1, 'Found %d temperature sensors channel%s.\n', ...
    num_temp_sensor_channels, plural(num_temp_sensor_channels));
fprintf(1, '\n');

% Determine how many samples the data file contains.

% Each data block contains 60 amplifier samples.
bytes_per_block = 60 * 4;  % timestamp data
bytes_per_block = bytes_per_block + 60 * 2 * num_amplifier_channels;
% Auxiliary inputs are sampled 4x slower than amplifiers
bytes_per_block = bytes_per_block + 15 * 2 * num_aux_input_channels;
% Supply voltage is sampled 60x slower than amplifiers
bytes_per_block = bytes_per_block + 1 * 2 * num_supply_voltage_channels;
% Board analog inputs are sampled at same rate as amplifiers
bytes_per_block = bytes_per_block + 60 * 2 * num_board_adc_channels;
% Board digital inputs are sampled at same rate as amplifiers
if (num_board_dig_in_channels > 0)
    bytes_per_block = bytes_per_block + 60 * 2;
end
% Board digital outputs are sampled at same rate as amplifiers
if (num_board_dig_out_channels > 0)
    bytes_per_block = bytes_per_block + 60 * 2;
end
% Temp sensor is sampled 60x slower than amplifiers
if (num_temp_sensor_channels > 0)
   bytes_per_block = bytes_per_block + 1 * 2 * num_temp_sensor_channels; 
end

% How many data blocks remain in this file?
data_present = 0;
bytes_remaining = filesize - ftell(fid);
header_size = ftell(fid);
if (bytes_remaining > 0)
    data_present = 1;
end

num_data_blocks = bytes_remaining / bytes_per_block;

num_amplifier_samples = 60 * num_data_blocks;
num_aux_input_samples = 15 * num_data_blocks;
num_supply_voltage_samples = 1 * num_data_blocks;
num_board_adc_samples = 60 * num_data_blocks;
num_board_dig_in_samples = 60 * num_data_blocks;
num_board_dig_out_samples = 60 * num_data_blocks;

record_time = num_amplifier_samples / sample_rate;

if (data_present)
    fprintf(1, 'File contains %0.3f seconds of data.  Amplifiers were sampled at %0.2f kS/s.\n', ...
        record_time, sample_rate / 1000);
    fprintf(1, '\n');
else
    fprintf(1, 'Header file contains no data.  Amplifiers were sampled at %0.2f kS/s.\n', ...
        sample_rate / 1000);
    fprintf(1, '\n');
end

if (data_present)
    
    % Pre-allocate memory for data.
    fprintf(1, 'Allocating memory for data...\n');

    % Read sampled data from file.
    fprintf(1, 'Reading data from file...\n');

    for i=1:num_data_blocks
        % In version 1.2, we moved from saving timestamps as unsigned
        % integeters to signed integers to accomidate negative (adjusted)
        % timestamps for pretrigger data.
        if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 2) ...
        || (data_file_main_version_number > 1))
            fread(fid, 60, 'int32');
        else
            fread(fid, 60, 'uint32');
        end
        if (num_amplifier_channels > 0)
            fread(fid, [60, num_amplifier_channels], 'uint16');
        end
        if (num_aux_input_channels > 0)
            fread(fid, [15, num_aux_input_channels], 'uint16');
        end
        if (num_supply_voltage_channels > 0)
            fread(fid, [1, num_supply_voltage_channels], 'uint16');
        end
        if (num_temp_sensor_channels > 0)
            fread(fid, [1, num_temp_sensor_channels], 'int16');
        end
        if (num_board_adc_channels > 0)
            fread(fid, [60, num_board_adc_channels], 'uint16');
        end
        if (num_board_dig_in_channels > 0)
            fread(fid, 60, 'uint16');
        end
        if (num_board_dig_out_channels > 0)
            fread(fid, 60, 'uint16');
        end

    end

    % Make sure we have read exactly the right amount of data.
    bytes_remaining = filesize - ftell(fid);
    if (bytes_remaining ~= 0)
        %error('Error: End of file not reached.');
    end
    
% Close data file.
fclose(fid);
    %%
    %THIS IS THE PART WHERE WE ADD ADDITIONAL FILES
    for i=1:(num_files-1)
        [file, path, filterindex] = ...
            uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'off');

        tic;
        filename_add = [path,file];
        fid_add = fopen(filename_add, 'r');

        s = dir(filename);
        filesize = s.bytes;

        % Check 'magic number' at beginning of file to make sure this is an Intan
        % Technologies RHD2000 data file.
        magic_number = fread(fid_add, 1, 'uint32');
        if magic_number ~= hex2dec('c6912702')
            error('Unrecognized file type.');
        end
        
        if i==1
        copyfile(filename,'CombinedFiles.rhd');
        fid = fopen('CombinedFiles.rhd', 'a+');
        end
        
        fseek(fid_add,header_size,'bof');
        
        for j=1:num_data_blocks
            % In version 1.2, we moved from saving timestamps as unsigned
            % integeters to signed integers to accomidate negative (adjusted)
            % timestamps for pretrigger data.
            if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 2) ...
            || (data_file_main_version_number > 1))
                toAdd = fread(fid_add, 60, 'int32');
                fwrite(fid,toAdd,'int32');
            else
                toAdd = fread(fid_add, 60, 'uint32');
                fwrite(fid,toAdd,'uint32');
            end
            if (num_amplifier_channels > 0)
                toAdd = fread(fid_add, [60, num_amplifier_channels], 'uint16');
                fwrite(fid,toAdd,'uint16');
            end
            if (num_aux_input_channels > 0)
                toAdd = fread(fid_add, [15, num_aux_input_channels], 'uint16');
                fwrite(fid,toAdd,'uint16');
            end
            if (num_supply_voltage_channels > 0)
                toAdd = fread(fid_add, [1, num_supply_voltage_channels], 'uint16');
                fwrite(fid,toAdd,'uint16');
            end
            if (num_temp_sensor_channels > 0)
                toAdd = fread(fid_add, [1, num_temp_sensor_channels], 'int16');
                fwrite(fid,toAdd,'int16');
            end
            if (num_board_adc_channels > 0)
                toAdd = fread(fid_add, [60, num_board_adc_channels], 'uint16');
                fwrite(fid,toAdd,'uint16');
            end
            if (num_board_dig_in_channels > 0)
                toAdd = fread(fid_add, 60, 'uint16');
                fwrite(fid,toAdd,'uint16');
            end
            if (num_board_dig_out_channels > 0)
                toAdd = fread(fid_add, 60, 'uint16');
                fwrite(fid,toAdd,'uint16');
            end

        end
     
        fclose(fid_add);
    end
    fclose(fid);
    
end


function a = fread_QString(fid)

% a = read_QString(fid)
%
% Read Qt style QString.  The first 32-bit unsigned number indicates
% the length of the string (in bytes).  If this number equals 0xFFFFFFFF,
% the string is null.

a = '';
length = fread(fid, 1, 'uint32');
if length == hex2num('ffffffff')
    return;
end
% convert length from bytes to 16-bit Unicode words
length = length / 2;

for i=1:length
    a(i) = fread(fid, 1, 'uint16');
end

return


function s = plural(n)

% s = plural(n)
% 
% Utility function to optionally plurailze words based on the value
% of n.

if (n == 1)
    s = '';
else
    s = 's';
end

return



