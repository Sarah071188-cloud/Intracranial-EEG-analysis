% HFA SWITCHING ANALYSIS | CONSISTENT RESEARCH VARIABLE NAMING
% Working variables use descriptive camelCase names (as in the IRASA script).
% Original MAT-file output variables and HFA helper functions are preserved.
% IMPORTANT: Not validated numerically against the original in MATLAB.
% Dependencies: FieldTrip and ZeroOnesCount.m

clearvars
close all
%% Inputs
baselineName = {'endpoint200s'}; % Name of the BL approach
subjectIDs=arrayfun(@num2str, 101:120, 'UniformOutput', false); % Subject number
siteNames = [repmat({'Calgary'}, 1, 14), {'CPMC'}, repmat({'WashU'}, 1, 5)];
%% Begining of HFA approach
% Subject-level reporting arrays (same dimensions as the original output).
all_trial_removal = nan(numel(subjectIDs), 4);
all_chan_removal = nan(numel(subjectIDs), 4);
all_chan_removalfinal = nan(numel(subjectIDs), 4);
allArtifactResults = cell(1, numel(subjectIDs));

for subjectIdx = 1:numel(subjectIDs)
    artifactResults = struct();
    chan_resp_all = cell(1, 2);
    trial_avg_all = cell(1, 2);
    freq_avg_all = cell(1, 2);
    % Mention the baseline time period
    if str2double(subjectIDs{subjectIdx}) < 116
        baselineWindowOptions = {[1.3 1.5]}; % BL time period
    else
        baselineWindowOptions = {[0.3 0.5]}; % BL time period
    end
    % Load conditionData
    load(strcat("Data\",siteNames{1,subjectIdx},subjectIDs{1,subjectIdx},"\Switching\Pre-processed\S",subjectIDs{1,subjectIdx},"_Switching_preprocessed_data_ItoE_external_switching.mat"))
    externalswitchData=data_preproc;
    clear data_preproc
    load(strcat("Data\",siteNames{1,subjectIdx},subjectIDs{1,subjectIdx},"\Switching\Pre-processed\S",subjectIDs{1,subjectIdx},"_Switching_preprocessed_data_EtoI_internal_switching.mat"))
    internalswitchData=data_preproc;
    clear data_preproc

    % save location
    outputDirectory=strcat("Analysis\",siteNames{1,subjectIdx},subjectIDs{1,subjectIdx},"\Switching\Analysis_3sec_updated_blockbootstrap_200rand_10BS_NAMING_TEST");
    % Check the existence of the Save folder
    if ~exist(outputDirectory, 'dir')
        mkdir(outputDirectory)
    end

    % save name
    outputPrefix=strcat('S',subjectIDs{1,subjectIdx},'_Switching_',baselineName{1,1});
    baselineWindow = baselineWindowOptions{1,1}; % [Start End]
    switchWindow = [0 3]; % [Start End]

    %% Select the needed trials
    data1_ind_abc_option =  externalswitchData.external_switching_trials; % Selecting only trials in which participants successfully switched to external attention
    data1_ind_abc_option(isnan(data1_ind_abc_option))=[];
    data2_ind_abc_option =  internalswitchData.internal_switching_trials; % Selecting only trials in which participants successfully switched to internal attention
    data2_ind_abc_option(isnan(data2_ind_abc_option))=[];
    %% Select trials
    % Condition 1 conditionData
    cfg=[];
    cfg.trials=data1_ind_abc_option;
    externalswitchData = ft_selectdata(cfg, externalswitchData);
    % Condition 2 conditionData
    cfg=[];
    cfg.trials=data2_ind_abc_option;
    internalswitchData = ft_selectdata(cfg, internalswitchData);
    % Time indices for switching conditionData
    switchStartIndex = nearest(externalswitchData.time{1,1}, switchWindow(1,1));
    switchEndIndex = nearest(internalswitchData.time{1,1}, switchWindow(1,2))-1;

    %% Getting their equivalent network names
    yeoNetworkLabels = externalswitchData.yeonetwork;
    %%
    analysisMode = 2; % Original workflow: analyze conditions separately.
    assert(analysisMode == 2, ...
        "Combined-condition mode requires separate baseline handling and is not supported by this refactor.");

    for conditionIdx = 1:analysisMode
        % Combine trials across conditions to do HFA analysis
        if analysisMode == 1
            conditionData=externalswitchData;
            conditionData.trial=[conditionData.trial,internalswitchData.trial];
            conditionData.time=[conditionData.time,internalswitchData.time];
            conditionData.trialinfo=[conditionData.trialinfo;internalswitchData.trialinfo];
            conditionData.sampleinfo=[conditionData.sampleinfo;internalswitchData.sampleinfo];
        elseif analysisMode==2 && conditionIdx==1
            conditionData = externalswitchData;
            % Baseline
            if str2double(subjectIDs{subjectIdx}) ~= 109 && str2double(subjectIDs{subjectIdx}) < 116
                % Load baseline
                load(strcat("Data\",siteNames{1,subjectIdx},subjectIDs{1,subjectIdx},"\Switching\Pre-processed\S",subjectIDs{1,subjectIdx},"_Switching_preprocessed_data_baseline.mat"))
                baselineData=data_preproc;
                % Select baseline time period
                baselineStartIndex = nearest(baselineData.time{1,1}, baselineWindow(1,1))+1;
                baselineEndIndex = nearest(baselineData.time{1,1}, baselineWindow(1,2));
                baselineIndices=baselineStartIndex:baselineEndIndex;
            else
                % select number of artifact trials in each session for Washu dataset
                if str2double(subjectIDs{subjectIdx}) > 115
                    if conditionData.trialsfirstsession+conditionData.trialssecondsession ~= size(data1_ind_abc_option,1)
                        tmp=1:(conditionData.trialsfirstsession+conditionData.trialssecondsession);
                        [uniqueVals, ~, idx]=unique([tmp';data1_ind_abc_option]);

                        counts = histcounts(idx, 1:numel(uniqueVals)+1);
                        nonRepeating = uniqueVals(counts == 1);
                        for tmp_t=1:size(nonRepeating,1)
                            if nonRepeating(tmp_t)<=conditionData.trialsfirstsession
                                conditionData.trialsfirstsession=conditionData.trialsfirstsession-1;
                            else
                                conditionData.trialssecondsession=conditionData.trialssecondsession-1;

                            end
                        end
                    end
                end
                % Load baseline
                load(strcat("Data\",siteNames{1,subjectIdx},subjectIDs{1,subjectIdx},"\Switching\Pre-processed\S",subjectIDs{1,subjectIdx},"_Switching_preprocessed_data_baseline_external_switch.mat"))
                baselineData=data_preproc;
                % Select baseline time period
                baselineStartIndex = nearest(baselineData.time{1,1}, baselineWindow(1,1))+1;
                baselineEndIndex = nearest(baselineData.time{1,1}, baselineWindow(1,2));
                baselineIndices=baselineStartIndex:baselineEndIndex;
            end
            clear data_preproc
        elseif analysisMode==2 && conditionIdx==2
            conditionData = internalswitchData;
            % Baseline
            if str2double(subjectIDs{subjectIdx}) ~= 109 && str2double(subjectIDs{subjectIdx}) < 116
                % Load baseline
                load(strcat("Data\",siteNames{1,subjectIdx},subjectIDs{1,subjectIdx},"\Switching\pre-processed\S",subjectIDs{1,subjectIdx},"_Switching_preprocessed_data_baseline.mat"))
                baselineData=data_preproc;
                % Select baseline time period
                baselineStartIndex = nearest(baselineData.time{1,1}, baselineWindow(1,1))+1;
                baselineEndIndex = nearest(baselineData.time{1,1}, baselineWindow(1,2));
                baselineIndices=baselineStartIndex:baselineEndIndex;
            else
                % select number of artifact trials in each session for Washu dataset
                if str2double(subjectIDs{subjectIdx}) > 115
                    if conditionData.trialsfirstsession+conditionData.trialssecondsession ~= size(data2_ind_abc_option,1)
                        tmp=1:(conditionData.trialsfirstsession+conditionData.trialssecondsession);
                        [uniqueVals, ~, idx]=unique([tmp';data2_ind_abc_option]);

                        counts = histcounts(idx, 1:numel(uniqueVals)+1);
                        nonRepeating = uniqueVals(counts == 1);
                        for tmp_t=1:size(nonRepeating,1)
                            if nonRepeating(tmp_t)<=conditionData.trialsfirstsession
                                conditionData.trialsfirstsession=conditionData.trialsfirstsession-1;
                            else
                                conditionData.trialssecondsession=conditionData.trialssecondsession-1;

                            end
                        end
                    end
                end
                % Load baseline
                load(strcat("Data\",siteNames{1,subjectIdx},subjectIDs{1,subjectIdx},"\Switching\pre-processed\S",subjectIDs{1,subjectIdx},"_Switching_preprocessed_data_baseline_internal_switch.mat"))
                baselineData=data_preproc;
                % Select baseline time period
                baselineStartIndex = nearest(baselineData.time{1,1}, baselineWindow(1,1))+1;
                baselineEndIndex = nearest(baselineData.time{1,1}, baselineWindow(1,2));
                baselineIndices=baselineStartIndex:baselineEndIndex;
            end
        end
        % calling function file to do HFA analysis
        if str2double(subjectIDs{subjectIdx}) < 116
            [chan_resp_all{1,conditionIdx},trial_avg_all{1,conditionIdx},freq_avg_all{1,conditionIdx},data_st_bl] = HFA_analysis(conditionData,baselineData,baselineIndices,[switchStartIndex:switchEndIndex]);%#ok<AGROW>
        else
            % Additional parameters needed for conditionData collected from WASHU site, since the task was performed by splitting into four different sessions
            % Switch to external split into 2 sessions and switch to internal split into two sessions
            % conditionData.trialsfirstsession = No. of trials selected from first half of the switching to external/internal session
            % conditionData.trialssecondsession = No. of trials selected from second half of the switching to external/internal session

            clear tmp uniqueVals idx counts nonRepeating
            [chan_resp_all{1,conditionIdx},trial_avg_all{1,conditionIdx},freq_avg_all{1,conditionIdx},data_st_bl] = HFA_analysis_WashU(conditionData,baselineData,baselineIndices,[switchStartIndex:switchEndIndex],conditionData.trialsfirstsession,conditionData.trialssecondsession);
        end
    end
    clear chan_resp trial_avg HFA_avg_data1 HFA_avg_data2
    %%  Determine the active electrodes
    % For separate analysis of switching we pick the common active electrodes
    % across both conditions.
    % For combined analysis (One dataset of switching) the channel responsive
    % electrodes are just the results we got from the above function file.
    if analysisMode == 1
        chan_resp= chan_resp_all{1,1};
        trial_avg = trial_avg_all{1,1};
        % grouping of trials into different conditions
        HFA_avg_data1 = freq_avg_all{1,1}(1:size(externalswitchData.trial,2),:,:);
        HFA_avg_data2 = freq_avg_all{1,1}(size(externalswitchData.trial,2)+1:end,:,:);
    elseif analysisMode == 2
        % If the HFA analysis is for both condions are performed separately, then the
        % unique active electrodes of both the conditions are selected.
        chan_resp_temp= [chan_resp_all{1,1},chan_resp_all{1,2}];
        UniqueC = unique(chan_resp_temp);
        chan_resp = UniqueC;
        trial_avg = trial_avg_all{1,1};
        HFA_avg_data1 = freq_avg_all{1,1};
        HFA_avg_data2 = freq_avg_all{1,2};
    end


    %% Map responsive electrodes to the original Yeo network definitions
    % Indices refer to positions in chan_resp, not directly to conditionData.label.
    [ch_resp_ind, yeo_net_sel, networkIndex] = ...
        classify_responsive_channels(chan_resp, conditionData.label, yeoNetworkLabels);

    VANA_ind_Ori = networkIndex.VANA;
    VANB_ind_Ori = networkIndex.VANB;
    FPCNA_ind_Ori = networkIndex.FPCNA;
    FPCNB_ind_Ori = networkIndex.FPCNB;
    DMNA_ind_Ori = networkIndex.DMNA;
    DMNB_ind_Ori = networkIndex.DMNB;
    DMNC_ind_Ori = networkIndex.DMNC;
    DMND_ind_Ori = networkIndex.DMND;
    DANA_ind_Ori = networkIndex.DANA;
    DANB_ind_Ori = networkIndex.DANB;
    DMN_ind_Ori = networkIndex.DMN;
    DAN_ind_Ori = networkIndex.DAN;
    VAN_ind_Ori = networkIndex.VAN;
    FPCN_ind_Ori = networkIndex.FPCN;

    %% saving the output
    data_time=externalswitchData.time; % This information is needed to plot the results
    save(strcat(outputDirectory,'\',outputPrefix,'chan_responsive_Yeo_label_check'),"chan_resp","yeo_net_sel","ch_resp_ind","VANA_ind_Ori","VANB_ind_Ori","FPCNB_ind_Ori","FPCNA_ind_Ori","FPCN_ind_Ori","DMNA_ind_Ori", "DMNB_ind_Ori","DMNC_ind_Ori","DMND_ind_Ori","trial_avg", "HFA_avg_data1","HFA_avg_data2", "freq_avg_all","DANA_ind_Ori","DANB_ind_Ori","data_time","DMN_ind_Ori","VAN_ind_Ori","DAN_ind_Ori","data1_ind_abc_option","data2_ind_abc_option");

    %% Remove outlier electrodes
    switchTimeIndices = switchStartIndex:switchEndIndex;
    networkNames={'VAN','FPCN','FPCNA','FPCNB','VANA','VANB','DMN','DAN','DMNA','DMNB','DMNC','DMND','DANA','DANB'};

        networkLookup = struct('VAN', VAN_ind_Ori, 'FPCN', FPCN_ind_Ori, ...
            'FPCNA', FPCNA_ind_Ori, 'FPCNB', FPCNB_ind_Ori, ...
            'VANA', VANA_ind_Ori, 'VANB', VANB_ind_Ori, ...
            'DMN', DMN_ind_Ori, 'DMNA', DMNA_ind_Ori, ...
            'DMNB', DMNB_ind_Ori, 'DMNC', DMNC_ind_Ori, ...
            'DMND', DMND_ind_Ori, 'DAN', DAN_ind_Ori, ...
            'DANA', DANA_ind_Ori, 'DANB', DANB_ind_Ori);
    for l=1:size(networkNames,2)
        % Select the network's responsive-electrode indices.
        networkChannelIndices = networkLookup.(networkNames{l});

        maxExternalswitch=max(HFA_avg_data1(:,ch_resp_ind(networkChannelIndices),switchTimeIndices),[],3);
        maxInternalswitch=max(HFA_avg_data2(:,ch_resp_ind(networkChannelIndices),switchTimeIndices),[],3);
        % A channel is rejected if >5 trials exceed 25 in either condition.
        rejectedChannels = find(sum(maxExternalswitch > 25, 1) > 5 | ...
                                sum(maxInternalswitch > 25, 1) > 5);
        rejectedChannelMatrix = nan(1, max(25, numel(rejectedChannels)));
        rejectedChannelMatrix(1:numel(rejectedChannels)) = rejectedChannels;

        nNetworkChannels = size(maxExternalswitch, 2);
        all_elec(l, :) = nNetworkChannels; %#ok<AGROW>
        maxExternalswitch(:, rejectedChannels) = [];
        maxInternalswitch(:, rejectedChannels) = [];

        % A trial is rejected if any remaining channel exceeds 40.
        % Preserve original 25-column output layout when possible.
        rejectedExtTrials = find(any(maxExternalswitch > 40, 2));
        rejectedIntTrials = find(any(maxInternalswitch > 40, 2));
        rejectedTrialMatrix = nan(2, max([25, numel(rejectedExtTrials), ...
                              numel(rejectedIntTrials)]));
        rejectedTrialMatrix(1, 1:numel(rejectedExtTrials)) = rejectedExtTrials.';
        rejectedTrialMatrix(2, 1:numel(rejectedIntTrials)) = rejectedIntTrials.';

        artifactResults.(networkNames{1,l}).chan_rem=rejectedChannelMatrix;
        artifactResults.(networkNames{1,l}).chansize=nNetworkChannels;
        artifactResults.(networkNames{1,l}).finalchansize=size(maxExternalswitch,2);
        artifactResults.(networkNames{1,l}).trial_rem=rejectedTrialMatrix;



    end

    allArtifactResults{1,subjectIdx}=artifactResults;
    summaryNetworkIndices=[1,2,7,8];
    all_trial_removal_tmp1num=[]; all_trial_removal_tmp2num=[];
    for summaryIdx=1:4
        all_trial_removal_tmp1(1,summaryIdx)=sum(~isnan(artifactResults.(networkNames{1,summaryNetworkIndices(summaryIdx)}).trial_rem(1,:)));
        all_trial_removal_tmp2(1,summaryIdx)=sum(~isnan(artifactResults.(networkNames{1,summaryNetworkIndices(summaryIdx)}).trial_rem(2,:)));
        all_chan_removal_tmp(1,summaryIdx)=artifactResults.(networkNames{1,summaryNetworkIndices(summaryIdx)}).chansize;
        all_chan_removalfinal_tmp(1,summaryIdx)=artifactResults.(networkNames{1,summaryNetworkIndices(summaryIdx)}).finalchansize;
        all_trial_removal_tmp1num = [all_trial_removal_tmp1num;artifactResults.(networkNames{1,summaryNetworkIndices(summaryIdx)}).trial_rem(1,~isnan(artifactResults.(networkNames{1,summaryNetworkIndices(summaryIdx)}).trial_rem(1,:)))'];
        all_trial_removal_tmp2num = [all_trial_removal_tmp2num;artifactResults.(networkNames{1,summaryNetworkIndices(summaryIdx)}).trial_rem(2,~isnan(artifactResults.(networkNames{1,summaryNetworkIndices(summaryIdx)}).trial_rem(2,:)))'];

    end
    all_trial_removal_tmp1num=unique(all_trial_removal_tmp1num);
    all_trial_removal_tmp2num=unique(all_trial_removal_tmp2num);
    VAN_ind_Ori(allArtifactResults{1,subjectIdx}.VAN.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.VAN.chan_rem)),:)=[];
    VANA_ind_Ori(allArtifactResults{1,subjectIdx}.VANA.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.VANA.chan_rem)),:)=[];
    VANB_ind_Ori(allArtifactResults{1,subjectIdx}.VANB.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.VANB.chan_rem)),:)=[];
    FPCN_ind_Ori(allArtifactResults{1,subjectIdx}.FPCN.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.FPCN.chan_rem)),:)=[];
    FPCNA_ind_Ori(allArtifactResults{1,subjectIdx}.FPCNA.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.FPCNA.chan_rem)),:)=[];
    FPCNB_ind_Ori(allArtifactResults{1,subjectIdx}.FPCNB.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.FPCNB.chan_rem)),:)=[];
    DMN_ind_Ori(allArtifactResults{1,subjectIdx}.DMN.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.DMN.chan_rem)),:)=[];
    DMNA_ind_Ori(allArtifactResults{1,subjectIdx}.DMNA.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.DMNA.chan_rem)),:)=[];
    DMNB_ind_Ori(allArtifactResults{1,subjectIdx}.DMNB.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.DMNB.chan_rem)),:)=[];
    DMNC_ind_Ori(allArtifactResults{1,subjectIdx}.DMNC.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.DMNC.chan_rem)),:)=[];
    DMND_ind_Ori(allArtifactResults{1,subjectIdx}.DMND.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.DMND.chan_rem)),:)=[];
    DAN_ind_Ori(allArtifactResults{1,subjectIdx}.DAN.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.DAN.chan_rem)),:)=[];
    DANA_ind_Ori(allArtifactResults{1,subjectIdx}.DANA.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.DANA.chan_rem)),:)=[];
    DANB_ind_Ori(allArtifactResults{1,subjectIdx}.DANB.chan_rem(1,~isnan(allArtifactResults{1,subjectIdx}.DANB.chan_rem)),:)=[];




    all_trial_removal(subjectIdx,1)=size(all_trial_removal_tmp1num,1);
    all_trial_removal(subjectIdx,2)=size(HFA_avg_data1,1);
    all_trial_removal(subjectIdx,3)=size(all_trial_removal_tmp2num,1);
    all_trial_removal(subjectIdx,4)=size(HFA_avg_data2,1);
    all_chan_removal(subjectIdx,:) = all_chan_removal_tmp(1,:);
    all_chan_removalfinal(subjectIdx,:) = all_chan_removalfinal_tmp(1,:);

    data1_ind_abc_option(all_trial_removal_tmp1num)=[];
    data2_ind_abc_option(all_trial_removal_tmp2num)=[];
  
    save(strcat(outputDirectory,'\',outputPrefix,'chan_responsive_Yeo_label_meanoutlierremoved'),"chan_resp","yeo_net_sel","ch_resp_ind","VANA_ind_Ori","VANB_ind_Ori","FPCNB_ind_Ori","FPCNA_ind_Ori","FPCN_ind_Ori","DMNA_ind_Ori", "DMNB_ind_Ori","DMNC_ind_Ori","DMND_ind_Ori","trial_avg", "HFA_avg_data1","HFA_avg_data2", "freq_avg_all","DANA_ind_Ori","DANB_ind_Ori","data_time","DMN_ind_Ori","VAN_ind_Ori","DAN_ind_Ori","data1_ind_abc_option","data2_ind_abc_option","all_trial_removal_tmp1num","all_trial_removal_tmp2num");
end

%% Local helper: classify responsive channels without growing arrays
function [channelIndices, networkLabels, groups] = ...
    classify_responsive_channels(responsiveLabels, allLabels, yeoLabels)
% Use the original first substring match for compatibility. A stricter
% strcmp match would be preferable but may alter legacy channel assignments.
    n = numel(responsiveLabels);
    channelIndices = zeros(n, 1);
    networkLabels = cell(n, 2);
    names = {'VANA','VANB','FPCNA','FPCNB','DMNA','DMNB', ...
             'DMNC','DMND','DANA','DANB','DMN','DAN','VAN','FPCN'};
    groups = struct();
    for g = 1:numel(names)
        groups.(names{g}) = zeros(n, 1);
    end
    counts = zeros(1, numel(names));

    dmnPairs = {'DMN-A-DMN-D','DMN-A-DMN-C','DMN-A-DMN-B', ...
        'DMN-B-DMN-C','DMN-B-DMN-D','DMN-C-DMN-D', ...
        'DMN-D-DMN-A','DMN-C-DMN-A','DMN-B-DMN-A', ...
        'DMN-C-DMN-B','DMN-D-DMN-B','DMN-D-DMN-C'};
    pairLabels = {dmnPairs, {'DAN-A-DAN-B','DAN-B-DAN-A'}, ...
        {'VAN-A-VAN-B','VAN-B-VAN-A'}, ...
        {'FPCN-A-FPCN-B','FPCN-B-FPCN-A'}};
    pairNames = {'DMN','DAN','VAN','FPCN'};
    singleLabels = {'VAN-A','VAN-B','FPCN-A','FPCN-B', ...
        'DMN-A','DMN-B','DMN-C','DMN-D','DAN-A','DAN-B'};
    singleNames = {'VANA','VANB','FPCNA','FPCNB', ...
        'DMNA','DMNB','DMNC','DMND','DANA','DANB'};

    for ch = 1:n
        match = find(~cellfun('isempty', strfind(allLabels, responsiveLabels{ch})), 1, 'first');
        assert(~isempty(match), 'Responsive channel not found: %s', responsiveLabels{ch});
        channelIndices(ch) = match;
        networkLabels{ch,1} = yeoLabels{match};
        networkLabels{ch,2} = yeoLabels{match};
        label = char(string(yeoLabels{match}));
        groupName = '';
        index = find(strcmp(label, singleLabels), 1);
        if ~isempty(index)
            groupName = singleNames{index};
        elseif numel(label) > 6
            for p = 1:numel(pairNames)
                if any(strcmp(label, pairLabels{p}))
                    groupName = pairNames{p};
                    break;
                end
            end
        end
        if ~isempty(groupName)
            g = find(strcmp(names, groupName), 1);
            counts(g) = counts(g) + 1;
            values = groups.(groupName);
            values(counts(g)) = ch;
            groups.(groupName) = values;
        end
    end
    for g = 1:numel(names)
        values = groups.(names{g});
        groups.(names{g}) = values(1:counts(g));
    end
    groups.DMN = [groups.DMN;groups.DMNA;groups.DMNB;groups.DMNC;groups.DMND];
    groups.DAN = [groups.DAN;groups.DANA;groups.DANB];
    groups.VAN = [groups.VAN;groups.VANA;groups.VANB];
    groups.FPCN = [groups.FPCN;groups.FPCNA;groups.FPCNB];
end

%% Bandpass filter and hilbert transformation
function [chan_resp,trial_avg,freq_avg,data1] = HFA_analysis(data,data_bl,t_BL,t_TT)
%% HFA analysis
% Preallocate double-precision Hilbert-amplitude arrays. Retain the original
% nine frequency bins (70-80 through 150-160 Hz) and filter settings.
nBands = numel(70:10:150);
y = zeros(numel(data.trial), numel(data.label), nBands, size(data.trial{1,1},2));
ybl = zeros(1, numel(data_bl.label), nBands, size(data_bl.trial{1,1},2));
% Hilbert transform
for j=1:size(data.trial,2)
    c=1;
    for i = 70:10:150 % High frequency band
       
        x1=ft_preproc_bandpassfilter(data.trial{1,j}(:,:), 1000,...
            [i i+10],...
            4, 'but', 'twopass', 'no');

        y(j,:,c,:)=ft_preproc_hilbert(x1); % Hilbert transform data for each trial

        x2=ft_preproc_bandpassfilter(data_bl.trial{1,1}(:,:), 1000,...
            [i i+10],...
            4, 'but', 'twopass', 'no');
        ybl(1,:,c,:)=ft_preproc_hilbert(x2); % Hilbert transform data for baseline


        c=c+1; % increment
    end
end

% Baseline correction
baseline_mtx = ybl(:,:,:,t_BL);
%% Bootstrap
iter = 1000; % Permutation
nrand = 200; % Number of Randomized values
ntrial=size(y,1); % Number of trials
nchannel=size(baseline_mtx,2); % Number ofchannels
nfrq=size(baseline_mtx,3); % Number of frequency
y2=y; % Duplicating the values to maintain the padding value of each trial after normalization
t=size(y,4); % Size of the trial time stamps
block_size=10; % Determine block size for blockbootstrap
n_blocks = ceil(nrand / block_size); % Select number of blocks
bsd_baseline = zeros(1, iter); % Reused for every channel and frequency
bl_meanstd = zeros(nchannel, nfrq, 2);
% Block bootstrap
for i=1:nchannel
    for j=1:nfrq
        baseline=baseline_mtx(:,i,j,:);
        baseline=baseline(:);
        N = length(baseline);
        for k = 1:iter
            block_starts = randi([1, N-block_size+1], n_blocks, 1);
            rand_baseline = zeros(n_blocks * block_size, 1);
            for b = 1:n_blocks
                idx = block_starts(b):(block_starts(b)+block_size-1);
                rand_baseline((b-1)*block_size+1:b*block_size) = baseline(idx);
            end
            rand_baseline = rand_baseline(1:nrand);
            bsd_baseline(k) = mean(rand_baseline);
        end
        baseline_mean = repmat(mean(bsd_baseline), [ntrial 1 1 t]);
        baseline_std  = repmat(std(bsd_baseline), [ntrial 1 1 t]);
        bl_meanstd(i,j,1) = mean(bsd_baseline);
        bl_meanstd(i,j,2) = std(bsd_baseline);
        y2(:,i,j,:)=(y(:,i,j,:)-baseline_mean)./baseline_std;
    end
end
data1.bl_meanstd = bl_meanstd;
data1.blc_data = y2;% Baseline corrected data
data1.data = y; % Actual Hilbert transform data
%% Mean zscore value sustained over 50ms
chans_resp = cell(1, 0); % Handles subjects with no responsive electrodes
k=1;
freq_avg = squeeze(mean(y2(:,:,:,:),3));
trial_avg = squeeze(mean(freq_avg ,1));
% trial_avg(trial_avg>25)=nan;
zthr = 1.96;
dur = 0.05;
mask = trial_avg (:,t_TT)>= zthr;
st=data.time{1,1}(1,t_TT);
for c=1:size(mask,1)
    [start, len, k1] = ZeroOnesCount(mask(c,:)); % ZeroOnesCount.m publicly available function file was used for the analysis (by Derek O'Connor 2011)
    for l=1:k1
        tmp_dur =  (st(start(1,l)+len(1,l)-1)) - (st(start(1,l))) ;
        if tmp_dur >= dur

            chans_resp(k) = data.label(c);
            k = k+1;
        end
    end
end
chan_resp=unique(chans_resp);
end
%% Bandpass filter and hilbert transformation
function [chan_resp,trial_avg,freq_avg,data1] = HFA_analysis_WashU(data,data_bl,t_BL,t_TT,t1,t2)
%% HFA analysis
% Preallocate double arrays without altering the frequency bins or precision.
nBands = numel(70:10:150);
y = zeros(numel(data.trial), numel(data.label), nBands, size(data.trial{1,1},2));
ybl = zeros(numel(data_bl.trial), numel(data_bl.label), nBands, size(data_bl.trial{1,1},2));
% Hilbert transform for switching data
for j=1:size(data.trial,2)
    c=1;
    for i = 70:10:150 % High frequency band

        x1=ft_preproc_bandpassfilter(data.trial{1,j}(:,:), 1000,...
            [i i+10],...
            4, 'but', 'twopass', 'no');

        y(j,:,c,:)=ft_preproc_hilbert(x1); % Hilbert transform data for each trial

        c=c+1;
    end
end

% Hilbert transform for Baseline
for j=1:size(data_bl.trial,2)
    c=1;
    for i = 70:10:150 % High frequency band

        x2=ft_preproc_bandpassfilter(data_bl.trial{1,j}(:,:), 1000,...
            [i i+10],...
            4, 'but', 'twopass', 'no');
        ybl(j,:,c,:)=ft_preproc_hilbert(x2); % Hilbert transform data for baseline

        c=c+1;
    end
end
% Baseline correction
baseline_mtx = ybl(:,:,:,t_BL); % Select the baseline time point

%% Bootstrap
iter = 1000; % Permutation
nrand = 200; % Number of Randomized values
ntrial=size(y,1); % Number of trials
nchannel=size(baseline_mtx,2); % Number ofchannels
nfrq=size(baseline_mtx,3); % Number of frequency
y2=y; % Duplicating the values to maintain the padding value of each trial after normalization
t=size(y,4); % Size of the trial time stamps
ntrial=size([1:t1],2);
block_size=10; % Determine block size for blockbootstrap
n_blocks = ceil(nrand / block_size); % Select number of blocks
bsd_baseline = zeros(1, iter); % Reused for every channel and frequency
bl_meanstd = zeros(nchannel, nfrq, 2);
% Block bootstrap using baseline corresponding to first half of the switch to external/internal data
for i=1:nchannel
    for j=1:nfrq
        baseline=baseline_mtx(1,i,j,:);
        baseline=baseline(:);
        N = length(baseline);
        for k = 1:iter
            block_starts = randi([1, N-block_size+1], n_blocks, 1);
            rand_baseline = zeros(n_blocks * block_size, 1);
            for b = 1:n_blocks
                idx = block_starts(b):(block_starts(b)+block_size-1);
                rand_baseline((b-1)*block_size+1:b*block_size) = baseline(idx);
            end
            rand_baseline = rand_baseline(1:nrand);
            bsd_baseline(k) = mean(rand_baseline);

        end
        baseline_mean = repmat(mean(bsd_baseline), [ntrial 1 1 t]);
        baseline_std  = repmat(std(bsd_baseline), [ntrial 1 1 t]);
        bl_mean = repmat(mean(bsd_baseline), [1 1 1 size(t_BL,2)]);
        bl_std  = repmat(std(bsd_baseline), [1 1 1 size(t_BL,2)]);
        bl_meanstd(i,j,1) = mean(bsd_baseline);
        bl_meanstd(i,j,2) = std(bsd_baseline);
        y2(1:t1,i,j,:)=(y(1:t1,i,j,:)-baseline_mean)./baseline_std; % Baseline standardisation
        baseline_mtx_real (:,i,j,:) = (baseline_mtx(1,i,j,:)-bl_mean)./bl_std;
        baseline_mtx_std (:,i,j,:) = (baseline_mtx(:,i,j,:)-bl_mean)./bl_std;
    end
end

clear baseline_mtx_real baseline_mtx_std b1_mean b1_std rand_baseline baseline t
% Block bootstrap using baseline corresponding to first half of the switch to external/internal data
t=size(y,4); % Size of the trial time stamps
ntrial=size([t1+1:t1+t2],2); % Trial details of second half
for i=1:nchannel
    for j=1:nfrq
        baseline=baseline_mtx(2,i,j,:);
        baseline=baseline(:);
        N = length(baseline);
        for k = 1:iter
            block_starts = randi([1, N-block_size+1], n_blocks, 1);
            rand_baseline = zeros(n_blocks * block_size, 1);
            for b = 1:n_blocks
                idx = block_starts(b):(block_starts(b)+block_size-1);
                rand_baseline((b-1)*block_size+1:b*block_size) = baseline(idx);
            end
            rand_baseline = rand_baseline(1:nrand);
            bsd_baseline(k) = mean(rand_baseline);

        end
        baseline_mean = repmat(mean(bsd_baseline), [ntrial 1 1 t]);
        baseline_std  = repmat(std(bsd_baseline), [ntrial 1 1 t]);
        bl_mean = repmat(mean(bsd_baseline), [1 1 1 size(t_BL,2)]);
        bl_std  = repmat(std(bsd_baseline), [1 1 1 size(t_BL,2)]);
        bl_meanstd(i,j,1) = mean(bsd_baseline);
        bl_meanstd(i,j,2) = std(bsd_baseline);
        y2(t1+1:end,i,j,:)=(y(t1+1:end,i,j,:)-baseline_mean)./baseline_std; % Baseline standardisation
        baseline_mtx_real (:,i,j,:) = (baseline_mtx(1,i,j,:)-bl_mean)./bl_std;
        baseline_mtx_std (:,i,j,:) = (baseline_mtx(:,i,j,:)-bl_mean)./bl_std;
    end
end
data1.bl_meanstd = bl_meanstd;
data1.bl_data = y2;
data1.data = y;
data1.bl_bStrap_arti=baseline_mtx_std;
%% Mean zscore value sustained over 50ms
chans_resp = cell(1, 0); % Handles subjects with no responsive electrodes
k=1;

freq_avg = squeeze(mean(y2(:,:,:,:),3));
trial_avg = squeeze(mean(freq_avg ,1));
trial_avg(trial_avg>25)=nan;
zthr = 1.96;
dur = 0.05;
mask = trial_avg (:,t_TT)>= zthr;
st=data.time{1,1}(1,t_TT);
for c=1:size(mask,1)
    [start, len, k1] = ZeroOnesCount(mask(c,:)); % ZeroOnesCount.m publicly available function file was used for the analysis (by Derek O'Connor 2011)
    for l=1:k1
        tmp_dur =  (st(start(1,l)+len(1,l)-1)) - (st(start(1,l))) ;
        if tmp_dur >= dur

            chans_resp(k) = data.label(c);
            k = k+1;
        end
    end
end
chan_resp=unique(chans_resp);
end








