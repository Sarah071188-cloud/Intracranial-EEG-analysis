%% CNN - Step 1 - Extract FOOOF

% PURPOSE: implement fooof (to serve as input features to 1D-SCNN)

% Input: data_preproc_fooof.mat, trials_by_dimension.mat
% Output: fooof variable saved per subject per session, ready for CNN

clear all


%% Set Directories

% Add path to Fieldtrip
%addpath('');
%ft_defaults;

% Set main directory and add path
maindir = ('/Users/Data/');
addpath(genpath('/Users/Data/'))

% Set other directories
cnndir = [maindir 'stats/cnn/'];
featuresdir = [cnndir 'features/'];
if ~exist(cnndir); mkdir(cnndir); end
if ~exist(featuresdir); mkdir(featuresdir); end
dlm = '/';


%% Set up relevant variables

% List subject, session and dimension names
subjnames = {'s101','s102','s103','s104','s105','s106','s107'};
sessnames = {'d1','d2','d3','d4','d5','d6','d7'};
dimnames = {'off','free','stuck','goal','int','self','other','visual','auditory'};
dimnames_class1 = {'off_trials','free_trials','stuck_trials','goal_trials','int_trials','self_trials','other_trials'};
dimnames_class2 = {'on_trials','notfree_trials','notstuck_trials','notgoal_trials','ext_trials','notself_trials','notother_trials'};

% Enter subjects, sessions and dimensions to process
subj_list = 1:7;
sess_list = 1:7;
dim_list = 1:7;
chan_list = 1:32;


%% Extract fooof

for subj = subj_list
    
    subjectid = subjnames{subj};
    
    for sess = sess_list
        
        % set paths
        sessionid = sessnames{sess};
        preprocdir = [maindir dlm subjectid dlm sessionid dlm 'Task/Preprocessing'];
        cd(preprocdir);
        
        % load preprocessed data
        load([preprocdir dlm 'data_preproc_fooof.mat'])
        data = data_preproc_car; clear data_preproc*
        
        % select relevant time points of interest
        cfg             = [];
        cfg.latency     = [0 2];
        data            = ft_selectdata(cfg, data);
        
        % Implement FOOOF at trial level 
        % ~ this yields a fooof output per trial, which then gets
        % cateroized into classes of thought dimensions in the next step
        % ~ loop through each trial
        for tr = 1:size(data.trial,2)
            
            % ~~ extract info required by fieldtrip data structure
            % ~~ reformat single trial into fieldtrip format
            trial{1,1}      = data.trial{1,tr};
            time{1,1}       = data.time{1,1};
            data_ft.trial   = trial;
            data_ft.label   = data.label;
            data_ft.time    = time;
            data_ft.fsample = data.fsample;
            data_ft.sampleinfo = data.sampleinfo(tr,:);
            
            % ~ compute fractal and original
            cfg             = [];
            cfg.channel     = 'all';
            cfg.foi         = 4:1:30;
            cfg.taper       = 'hanning';
            cfg.pad         = 'nextpow2';
            cfg.method      = 'mtmfft';
            cfg.output      = 'fooof_aperiodic';
            fractal         = ft_freqanalysis(cfg,data_ft);
            cfg.output      = 'pow';
            original        = ft_freqanalysis(cfg,data_ft);
            
            % ~ compute 'true' oscillatory component (subtract fractal component from power spectrum)
            cfg                 = [];
            cfg.parameter       = 'powspctrm';
            cfg.operation       = 'x2-x1';
            oscillatory{1,tr}  	= ft_math(cfg, fractal, original);
            
            % extract power spectrum field
            oscillatory_ps(tr,:,:) = oscillatory{1,tr}.powspctrm;
            
            clear trial time data_ft fractal original
            
        end
        
        % save data
        filename = strcat(subjnames{subj},'_',sessnames{sess},'_fooof_pertrial.mat');
        save([featuresdir filename],'oscillatory','oscillatory_ps');
        clear data oscillatory* preprocdir
        
    end
    
end


%% Categorize trials into thought dimension classes
% store fooof data per subject (for all thought dimensions and across all sessions)

for dim = dim_list
    
    for subj = subj_list
        
        subjectid = subjnames{subj};
        class1_data=[]; class2_data=[];
        
        for sess = sess_list
            
            sessionid = sessnames{sess};
            preprocdir = [maindir dlm subjectid dlm sessionid dlm 'Task/Preprocessing'];
            cd(preprocdir);
            
            % load trials for both classes of given dimension
            load([preprocdir dlm 'trials_by_dimension.mat'],dimnames_class1{dim},dimnames_class2{dim});
            class1_trials = evalin('base',dimnames_class1{dim});
            class2_trials = evalin('base',dimnames_class2{dim});
            
            % load FOOOF power spectrum
            filename = strcat(subjnames{subj},'_',sessnames{sess},'_fooof_pertrial.mat');
            load([featuresdir filename],'oscillatory_ps');
            
            % Select relevant trials from oscillory_ps
            if ~isempty(class1_trials) && ~isempty(class2_trials)
                tmp_data1 = oscillatory_ps(class1_trials(:,1),:,:);
                tmp_data2 = oscillatory_ps(class2_trials(:,1),:,:);
            elseif isempty(class1_trials) && ~isempty(class2_trials)
                tmp_data2 = oscillatory_ps(class2_trials(:,1),:,:);
                tmp_data1 = [];
            elseif ~isempty(class1_trials) && isempty(class2_trials)
                tmp_data1 = oscillatory_ps(class1_trials(:,1),:,:);
                tmp_data2 = [];
            elseif isempty(class2_trials) && isempty(class1_trials)
                tmp_data1 = [];
                tmp_data2 = [];
            end
            
            % store data across sessions into one variable per subject
            class1_data = [class1_data; tmp_data1];
            class2_data = [class2_data; tmp_data2];
            clear tmp_data* *trials filename oscillatory_ps
        end
        
        class1_dimname = dimnames_class1{dim};
        class2_dimname = dimnames_class2{dim};
        
        % Save data
        filename = strcat(subjnames{subj},'_',dimnames{dim},'_fooof_forcnn.mat');
        save([featuresdir filename],'class1_dimname','class2_dimname','class1_data','class2_data');
        clear class*dimname class*data
        
    end
    
end
cd(featuresdir)
