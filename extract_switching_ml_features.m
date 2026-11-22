%% Extract network-averaged mean-power features for switching ML
% Output compatibility: Final.<network>.all_feat.c1/c2 and label_name.
clearvars;

subjectIDs = arrayfun(@num2str, 101:120, 'UniformOutput', false);
siteNames = [repmat({'Calgary'}, 1, 14), {'CPMC'}, repmat({'WashU'}, 1, 5)];
frequencyNames = {'Theta','Alpha','Beta','LGamma','HGamma'};
frequencyIndices = {14:20, 21:24, 25:32, 33:41, []}; % To precisely get the frequency ranges
networkNames = {'FPCN','VAN','DAN','DMN'};
analysisFolder = 'Analysis_3sec_updated_blockbootstrap_200rand_10BS';

for subjectIdx = 1:numel(subjectIDs)
    subjectID = subjectIDs{subjectIdx};
    outputDirectory = fullfile('Analysis', [siteNames{subjectIdx}, subjectID], ...
        'Switching', analysisFolder);
    hfaFile = fullfile(outputDirectory, ...
        ['S', subjectID, '_Switching_endpoint200schan_responsive_Yeo_label_meanoutlierremoved.mat']);
    hfa = load(hfaFile);
    tfFile = fullfile(outputDirectory, 'LogR_last300ms_48log.mat');
    tf = load(tfFile, 'tfa_BL_corrected_class1', 'tfa_BL_corrected_class2');

    % HFA: remove rejected trials; retain original sample indices 2001:5000.
    externalHFA = hfa.HFA_avg_data1;
    internalHFA = hfa.HFA_avg_data2;
    externalHFA(hfa.all_trial_removal_tmp1num,:,:) = [];
    internalHFA(hfa.all_trial_removal_tmp2num,:,:) = [];

    for frequencyIdx = 1:numel(frequencyNames)
        frequencyName = frequencyNames{frequencyIdx};
        if frequencyIdx == 5
            externalSwitch = externalHFA;
            internalSwitch = internalHFA;
            timeIndices = 2001:5000;
        else
            bandIndices = frequencyIndices{frequencyIdx};
            % Mean over frequency dimension (trials x channels x time).
            externalSwitch = mean(tf.tfa_BL_corrected_class1.powspctrm(:,:,bandIndices,:), 3, 'omitnan');
            internalSwitch = mean(tf.tfa_BL_corrected_class2.powspctrm(:,:,bandIndices,:), 3, 'omitnan');
            externalSwitch = reshape(externalSwitch, size(externalSwitch,1), size(externalSwitch,2), size(externalSwitch,4));
            internalSwitch = reshape(internalSwitch, size(internalSwitch,1), size(internalSwitch,2), size(internalSwitch,4));
            timeIndices = 41:100;
        end

        Final = struct();
        for networkIdx = 1:numel(networkNames)
            networkName = networkNames{networkIdx};
            Final.(networkName) = [];
            originalIndices = hfa.([networkName, '_ind_Ori']);
            if isempty(originalIndices)
                continue;
            end
            channelIndices = hfa.ch_resp_ind(originalIndices,1);
            % Trial x channel x time -> trial x channel -> trial x 1.
            externalChannelMeans = mean(externalSwitch(:,channelIndices,timeIndices), 3, 'omitnan');
            internalChannelMeans = mean(internalSwitch(:,channelIndices,timeIndices), 3, 'omitnan');
            externalFeature = mean(externalChannelMeans, 2, 'omitnan');
            internalFeature = mean(internalChannelMeans, 2, 'omitnan');

            Final.(networkName).all_feat.c1 = externalFeature;
            Final.(networkName).all_feat.c2 = internalFeature;
            Final.(networkName).label_name = hfa.chan_resp(:,originalIndices);
        end
        save(fullfile(outputDirectory, ['S',subjectID,'_ML_feat_Switching_',frequencyName,'.mat']), 'Final');
    end
    fprintf('Finished subject %s\n', subjectID);
end
