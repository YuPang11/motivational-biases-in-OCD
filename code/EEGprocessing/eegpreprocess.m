%%
clear all; clc;
root_dir = 'H:\newdata';
out_dir = [root_dir, '\rename'];
cd(root_dir);
data_list = dir('*.set');
for di = 1:length(data_list)
    EEG = pop_loadset('filename', data_list(di).name, 'filepath', root_dir);

    % Extract event fields
    mffkey_tar1 = {EEG.event.mffkey_tar1};
    mffkey_tar2 = {EEG.event.mffkey_tar2};
    mffkey_tar3 = {EEG.event.mffkey_tar3};
    mffkey_tar4 = {EEG.event.mffkey_tar4};
    mffkey_tar5 = {EEG.event.mffkey_tar5};
    mffkey_tar6 = {EEG.event.mffkey_tar6};

    targetValues = {'111', '112', '220', '321', '322', '410'};

    for i = 1:length(EEG.event)
        % Check mffkey_tar fields of the current event
        for j = 1:length(targetValues)
            if strcmp(mffkey_tar1{i}, targetValues{j}) || ...
               strcmp(mffkey_tar2{i}, targetValues{j}) || ...
               strcmp(mffkey_tar3{i}, targetValues{j}) || ...
               strcmp(mffkey_tar4{i}, targetValues{j}) || ...
               strcmp(mffkey_tar5{i}, targetValues{j}) || ...
               strcmp(mffkey_tar6{i}, targetValues{j})
               
                % Update event type with the matched value
                EEG.event(i).type = targetValues{j};
                break; % Exit inner loop after a match
            end
        end
    end
    
    res_mark = {EEG.event.mffkey_rsp};
    respIndexList = find(strcmp(res_mark, '1'));
    
    for i = 1:length(respIndexList)
        resp_label = strcat(EEG.event(respIndexList(i)).mffkey_cel, 'resp1');
        EEG.event(respIndexList(i)).type = resp_label;
    end
    
    fed_mark = {EEG.event.mffkey_argu};
    fedIndexList1 = find(strcmp(fed_mark, '66'));
    
    for i = 1:length(fedIndexList1)
        fed_label1 = strcat(EEG.event(fedIndexList1(i)).mffkey_cel, 'fed66');
        EEG.event(fedIndexList1(i)).type = fed_label1;
    end
    
    fedIndexList2 = find(strcmp(fed_mark, '88'));
    
    for i = 1:length(fedIndexList2)
        fed_label2 = strcat(EEG.event(fedIndexList2(i)).mffkey_cel, 'fed88');
        EEG.event(fedIndexList2(i)).type = fed_label2;
    end
    
    EEG = eeg_checkset(EEG);
    
    close all
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    cd(out_dir)
    EEG = pop_saveset(EEG, 'filename', [data_list(di).name]);
    EEG = eeg_checkset(EEG);
    cd(root_dir)
end




%%
%% Filtering
clear all; clc;
% Set folder paths
input_folder = 'H:\newdata\rename';
output_folder = 'H:\newdata\rename\filiter';

% Set filter parameters
low_freq = 0.1; % High-pass cutoff frequency
high_freq = 50; % High-frequency cutoff

% Get all .set files in the folder
file_pattern = fullfile(input_folder, '*.set');
files = dir(file_pattern);

% Create output folder
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Process each .set file
for i = 1:length(files)
    file_name = files(i).name;
    full_file_path = fullfile(input_folder, file_name);
    
    % Open file with EEGLAB
    EEG = pop_loadset('filename', file_name, 'filepath', input_folder);
    
    % Band-pass filter
    EEG = pop_eegfiltnew(EEG, low_freq, high_freq, [], 0);

    % Remove 50 Hz line noise
    EEG = pop_eegfiltnew(EEG, 50 - 1, 50 + 1, [], 1);
    
    % Rename file
    new_file_name = [file_name];
    
    % Save filtered data
    EEG = pop_saveset(EEG, 'filename', new_file_name, 'filepath', output_folder);
end



%%
%% Re-reference
clear all; clc;
% Set folder paths
root_dir = 'H:\newdata\rename\filiter';
out_dir = [root_dir, '\reference'];

cd(root_dir)
data_list=dir('*.set');

for di=1:length(data_list)
    % Load data
    EEG.etc.eeglabvers = '2022.0'; % This tracks the EEGLAB version and can be ignored
     EEG = pop_loadset('filename',data_list(di).name);
    
    % Average reference to bilateral mastoids
    EEG = pop_reref(EEG, [57 100]);
    
    
    EEG = eeg_checkset( EEG );
    
    % Save processed data to the new folder
close all
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    cd(out_dir)
    EEG = pop_saveset( EEG, 'filename',[data_list(di).name]);
    EEG = eeg_checkset( EEG );
    cd(root_dir)
end



%%
%% Epoching + baseline correction
clear all; clc;
root_dir='H:\newdata\rename\filiter\reference';
out_dir=[root_dir,'\cue_epoch'];
cd(root_dir)
data_list=dir('*.set');
for di=1:length(data_list)
    EEG.etc.eeglabvers = '2022.0'; % This tracks the EEGLAB version and can be ignored
    EEG = pop_loadset('filename',data_list(di).name);
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {'111' '112' '220' '321' '322' '410'}, [-1.75    4.5], 'newname', 'EEProbe continuous data epochs', 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [-200 0] ,[]);% Baseline correction
    EEG = eeg_checkset( EEG );


close all
    if~exist(out_dir,'file')
        mkdir(out_dir);
    end
    cd(out_dir)
    EEG = pop_saveset( EEG, 'filename',[data_list(di).name]);
    EEG = eeg_checkset( EEG );
    cd(root_dir)
end

%%
%% ICA
clear all; clc;

root_dir='H:\newdata\rename\filiter\reference\cue_epoch';
out_dir=[root_dir,'\out_put_ICA'];
cd(root_dir)
data_list=dir('*.set');
for di=1:length(data_list)

    EEG.etc.eeglabvers = '2022.0'; % This tracks the EEGLAB version and can be ignored
    EEG = pop_loadset('filename',data_list(di).name);

     EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on','pca',64);
     EEG = eeg_checkset( EEG );


close all
    if~exist(out_dir,'file')
        mkdir(out_dir);
    end
    cd(out_dir)
    EEG = pop_saveset( EEG, 'filename',[data_list(di).name]);
    EEG = eeg_checkset( EEG );
    cd(root_dir)
end
