%% PSI analysis for attention switching (IRASA peak-centred frequencies)
% Requires FieldTrip and Mike Cohen's data2psiX on the MATLAB path.
% Original analysis: 3 s window, +/-3 Hz per channel, 2 s zero padding.
clearvars;
close all;

subjectIDs = arrayfun(@num2str, 101:120, 'UniformOutput', false);
siteNames = [repmat({'Calgary'}, 1, 14), {'CPMC'}, repmat({'WashU'}, 1, 5)];
frequencyBands = {'Theta','Alpha','Beta'};
timeNames = {'0to3sec'};
timeWindows = {2001:5000};
samplingRate = 1000;
zeroPaddingSeconds = 2;
frequencyHalfWidth = 3;
analysisFolder = 'Analysis_3sec_updated_blockbootstrap_200rand_10BS';

for timeIdx = 1:numel(timeNames)
    sampleIndices = timeWindows{timeIdx};
    for bandIdx = 1:numel(frequencyBands)
        bandName = frequencyBands{bandIdx};
        bandField = lower(bandName);
        peakField = [bandField 'peak'];
        for subjectIdx = 1:numel(subjectIDs)
            subjectID = subjectIDs{subjectIdx};
            siteName = siteNames{subjectIdx};
            subjectName = ['S' subjectID];
            analysisDirectory = fullfile('Analysis', [siteName subjectID], ...
                'Switching', analysisFolder);
            saveDirectory = analysisDirectory;
            if ~exist(saveDirectory, 'dir')
                mkdir(saveDirectory);
            end
            % Keep the legacy filename (including its spelling and spacing).
            outputFile = fullfile(saveDirectory, sprintf( ...
                'S%s_ %s_PSI_output_for_swiching_IRASA_channels_1000iterations%s.mat', ...
                subjectID, bandName, timeNames{timeIdx}));
            psi_out = struct('data1', [], 'data2', [], ...
                'data1_statpsi', [], 'data2_statpsi', [], ...
                'data1_zpsi', [], 'data2_zpsi', []);

            hfaFile = fullfile(analysisDirectory, ...
                [subjectName '_Switching_endpoint200schan_responsive_Yeo_label_meanoutlierremoved.mat']);
            networkFile = fullfile(analysisDirectory, ...
                [subjectName '_IRASA_active_electrodes_network_meanoutlierremoved_Switching.mat']);
            irasaFile = fullfile(analysisDirectory, ...
                [subjectName '_Switching_IRASA_output.mat']);

            hfa = load(hfaFile, 'ch_resp_ind', 'yeo_net_sel');
            network = load(networkFile, 'out_osci_elec');
            irasa = load(irasaFile, 'oscillatory_class_all');

            networkNames = {'VAN','FPCN','DMN','DAN'};
            responsiveIndices = zeros(0,1);
            for networkIdx = 1:numel(networkNames)
                fieldName = [networkNames{networkIdx} '_' bandField];
                if ~isfield(network.out_osci_elec, fieldName) || ...
                        ~isfield(network.out_osci_elec.(fieldName), 'out')
                    error('Missing network field: %s in %s', fieldName, networkFile);
                end
                responsiveIndices = [responsiveIndices; ...
                    network.out_osci_elec.(fieldName).out(:)]; %#ok<AGROW>
            end
            responsiveIndices = unique(responsiveIndices);
            if numel(responsiveIndices) < 2
                warning('%s %s: fewer than two oscillatory channels. Saving empty output.', ...
                    subjectName, bandName);
                save(outputFile, 'psi_out', '-v7.3');
                continue;
            end

            channelIndices = hfa.ch_resp_ind(responsiveIndices, 1);
            channelLabels = irasa.oscillatory_class_all.label;
            peakFrequencies = network.out_osci_elec.(peakField);
            selectedPeaks = nan(numel(channelIndices), 1);
            externalPath = fullfile('Data', [siteName subjectID], 'Switching', ...
                'pre-processed', [subjectName '_Switching_preprocessed_data_ItoE_external_Switching.mat']);
            internalPath = fullfile('Data', [siteName subjectID], 'Switching', ...
                'pre-processed', [subjectName '_Switching_preprocessed_data_EtoI_internal_Switching.mat']);
            externalLoaded = load(externalPath, 'data_preproc');
            internalLoaded = load(internalPath, 'data_preproc');
            cfg = [];
            cfg.trials = externalLoaded.data_preproc.external_switching_trials;
            externalSwitch = ft_selectdata(cfg, externalLoaded.data_preproc);
            cfg = [];
            cfg.trials = internalLoaded.data_preproc.internal_switching_trials;
            internalSwitch = ft_selectdata(cfg, internalLoaded.data_preproc);

            for channelIdx = 1:numel(channelIndices)
                matches = find(strcmp(externalSwitch.label{channelIndices(channelIdx)}, channelLabels));
                if numel(matches) ~= 1
                    error('%s: expected one IRASA match for channel %s, found %d.', ...
                        subjectName, externalSwitch.label{channelIndices(channelIdx)}, numel(matches));
                end
                selectedPeaks(channelIdx) = peakFrequencies(matches);
            end
            if any(~isfinite(selectedPeaks))
                warning('%s %s: nonfinite IRASA peak; saving empty PSI output.', ...
                    subjectName, bandName);
                save(outputFile, 'psi_out', '-v7.3');
                continue;
            end
            frequencyBins = [selectedPeaks-frequencyHalfWidth, ...
                             selectedPeaks+frequencyHalfWidth];
            nChannels = numel(channelIndices);
            nSamples = numel(sampleIndices);
            nExternal = numel(externalSwitch.trial);
            nInternal = numel(internalSwitch.trial);
            if nExternal == 0 || nInternal == 0
                warning('%s: missing trials for one condition.', subjectName);
                save(outputFile, 'psi_out', '-v7.3');
                continue;
            end
            padSamples = samplingRate * zeroPaddingSeconds;
            externalArray = zeros(nChannels, nSamples+2*padSamples, nExternal);
            internalArray = zeros(nChannels, nSamples+2*padSamples, nInternal);
            for trialIdx = 1:nExternal
                externalArray(:, padSamples+(1:nSamples), trialIdx) = ...
                    externalSwitch.trial{trialIdx}(channelIndices, sampleIndices);
            end
            for trialIdx = 1:nInternal
                internalArray(:, padSamples+(1:nSamples), trialIdx) = ...
                    internalSwitch.trial{trialIdx}(channelIndices, sampleIndices);
            end

            externalPSI = data2psiX(externalArray, samplingRate, frequencyBins, 0);
            internalPSI = data2psiX(internalArray, samplingRate, frequencyBins, 0);
            externalZPSI = data2psiX(externalArray, samplingRate, frequencyBins, 1);
            internalZPSI = data2psiX(internalArray, samplingRate, frequencyBins, 1);

            labels = externalSwitch.label(channelIndices);
            networkLabels = hfa.yeo_net_sel(responsiveIndices);
            % Preserve the original saved metadata and output field names.
            externalResult = struct('label', {labels}, 'psispctrm', externalPSI, ...
                'freq', selectedPeaks, 'time', {externalSwitch.time{1}}, ...
                'nw', {networkLabels}, 'trialinfo', externalSwitch.trialinfo);
            internalResult = struct('label', {labels}, 'psispctrm', internalPSI, ...
                'freq', selectedPeaks, 'time', {internalSwitch.time{1}}, ...
                'nw', {networkLabels}, 'trialinfo', internalSwitch.trialinfo);
            psi_out.data1 = externalResult;
            psi_out.data2 = internalResult;
            psi_out.data1_zpsi = externalResult;
            psi_out.data1_zpsi.psispctrm = externalZPSI;
            psi_out.data1_zpsi.freqbins = frequencyBins;
            psi_out.data2_zpsi = internalResult;
            psi_out.data2_zpsi.psispctrm = internalZPSI;
            psi_out.data2_zpsi.freqbins = frequencyBins;
            save(outputFile, 'psi_out', '-v7.3');
            fprintf('Saved %s (%s, %d channels)\n', subjectName, bandName, nChannels);
        end
    end
end
