%% IRASA analysis of attention-state switching
% Computes oscillatory, fractal, and original spectra from combined
% external- and internal-switching trials. Existing analytical settings
% and output variable names are retained for downstream compatibility.
%
% Dependencies: FieldTrip, pow_irasa local function, Statistics Toolbox
% (randsample). The bootstrap operates across frequencies, not trials.

clearvars;

subjectIDs = arrayfun(@num2str, 101:120, 'UniformOutput', false);
siteNames = [repmat({'Calgary'}, 1, 14), {'CPMC'}, repmat({'WashU'}, 1, 5)];

% Analysis settings (preserved from original script)
analysisConfig.pad = 'nextpow2';
analysisConfig.method = 'irasa';
analysisConfig.foi = [3, 30];
analysisWindow = [0, 3];
bootstrapIterations = 1000;
bootstrapSampleCount = 5;

for subjectIdx = 1:numel(subjectIDs)
    subjectID = subjectIDs{subjectIdx};
    siteName = siteNames{subjectIdx};
    subjectFolder = [siteName, subjectID];
    fprintf('Processing subject %s (%s) [%d/%d]\n', ...
        subjectID, siteName, subjectIdx, numel(subjectIDs));

    saveDir = fullfile('Analysis', subjectFolder, 'Switching', ...
        'Analysis_3sec_updated_blockbootstrap_200rand_10BS');
    preprocessedDir = fullfile('Data', subjectFolder, 'Switching', 'pre-processed');
    externalswitchFile = fullfile(preprocessedDir, ...
        ['S', subjectID, '_Switching_preprocessed_data_ItoE_external_Switching.mat']);
    internalswitchFile = fullfile(preprocessedDir, ...
        ['S', subjectID, '_Switching_preprocessed_data_EtoI_internal_Switching.mat']);
    selectionFile = fullfile(saveDir, ...
        ['S', subjectID, '_Switching_endpoint200schan_responsive_Yeo_label_meanoutlierremoved.mat']);
    outputFile = fullfile(saveDir, ['S', subjectID, '_Switching_IRASA_output.mat']);

    % Load channel selection and retained trial indices from HFA analysis.
    selection = load(selectionFile, 'VAN_ind_Ori', 'FPCN_ind_Ori', ...
        'DMN_ind_Ori', 'DAN_ind_Ori', 'ch_resp_ind', ...
        'data1_ind_abc_option', 'data2_ind_abc_option');
    networkChannelPositions = unique([selection.VAN_ind_Ori; ...
        selection.FPCN_ind_Ori; selection.DMN_ind_Ori; selection.DAN_ind_Ori]);
    analysisConfig.chan = selection.ch_resp_ind(networkChannelPositions);

    % Select the same trials as in the original HFA-derived selection.
    externalswitchLoaded = load(externalswitchFile, 'data_preproc');
    cfg = [];
    cfg.trials = selection.data1_ind_abc_option;
    externalswitchData = ft_selectdata(cfg, externalswitchLoaded.data_preproc);
    clear externalswitchLoaded;

    internalswitchLoaded = load(internalswitchFile, 'data_preproc');
    cfg = [];
    cfg.trials = selection.data2_ind_abc_option;
    internalswitchData = ft_selectdata(cfg, internalswitchLoaded.data_preproc);
    clear internalswitchLoaded;

    % Crop each trial to [0, 3) seconds using the original nearest-index rule.
    externalswitchData = crop_switch_trials(externalswitchData, analysisWindow);
    internalswitchData = crop_switch_trials(internalswitchData, analysisWindow);

    % Combine conditions before spectral estimation (original behavior).
    combinedData = externalswitchData;
    combinedData.trial = [externalswitchData.trial, internalswitchData.trial];
    combinedData.time = [externalswitchData.time, internalswitchData.time];
    combinedData.trialinfo = [externalswitchData.trialinfo; internalswitchData.trialinfo];
    combinedData.sampleinfo = [externalswitchData.sampleinfo; internalswitchData.sampleinfo];

    % Preserve original sub-segmentation settings; verify overlap separately.
    [~, shortestTrialIdx] = min(cellfun('size', combinedData.trial, 2));
    trialDuration = combinedData.time{shortestTrialIdx}(end) - ...
        combinedData.time{shortestTrialIdx}(1);
    cfg = [];
    cfg.length = trialDuration * 0.9;
    cfg.overlap = 1 - ((trialDuration - cfg.length) / (10 - 1));
    combinedData = ft_redefinetrial(cfg, combinedData);

    % IRASA decomposition and original-spectrum calculation.
    [oscillatory_class_all, fractal_class_all, original_class_all] = ...
        pow_irasa(analysisConfig, combinedData);

    % Average across repetitions, as in the original nanmean call.
    oscillatory_class_all.powspctrm = ...
        squeeze(mean(oscillatory_class_all.powspctrm, 1, 'omitnan'));

    % Frequency-wise bootstrap z-scoring for each channel.
    oscillatory_class_all_bs = bs_psd(oscillatory_class_all, ...
        bootstrapIterations, bootstrapSampleCount);

    if ~isfolder(saveDir)
        mkdir(saveDir);
    end
    save(outputFile, 'oscillatory_class_all', 'oscillatory_class_all_bs', ...
        'original_class_all', 'fractal_class_all');
    fprintf('Saved: %s\n', outputFile);
end

%% Crop individual switching trials
function data = crop_switch_trials(data, timeWindow)
    for trialIdx = 1:numel(data.trial)
        startIdx = nearest(data.time{trialIdx}, timeWindow(1));
        endIdx = nearest(data.time{trialIdx}, timeWindow(2)) - 1;
        data.trial{trialIdx} = data.trial{trialIdx}(:, startIdx:endIdx);
        data.time{trialIdx} = data.time{trialIdx}(1, startIdx:endIdx);
    end
    data.label = data.label(:, 1); % Preserve original channel-label selection.
end

%% Bootstrap normalization of oscillatory spectrum
function zSpectrum = bs_psd(data, nIterations, nSamples)
    spectrum = data.powspctrm;
    [nChannels, nFrequencies] = size(spectrum);
    zSpectrum = zeros(nChannels, nFrequencies);
    bootstrapMeans = zeros(nIterations, 1);

    for channelIdx = 1:nChannels
        channelSpectrum = spectrum(channelIdx, :).';
        for iterationIdx = 1:nIterations
            sampledSpectrum = randsample(channelSpectrum, nSamples);
            bootstrapMeans(iterationIdx) = mean(sampledSpectrum);
        end
        baselineMean = mean(bootstrapMeans);
        baselineStd = std(bootstrapMeans);
        zSpectrum(channelIdx, :) = ...
            (spectrum(channelIdx, :) - baselineMean) ./ baselineStd;
    end
end

%% IRASA spectral decomposition
function [oscillatory, fractal, original] = pow_irasa(analysisConfig, data)
    cfg = [];
    cfg.foilim = analysisConfig.foi;
    cfg.taper = 'hanning';
    cfg.keeptrials = 'yes';
    cfg.channel = analysisConfig.chan;
    cfg.pad = analysisConfig.pad;
    cfg.method = analysisConfig.method;
    fractalSegments = ft_freqanalysis(cfg, data);

    cfg.method = 'mtmfft';
    originalSegments = ft_freqanalysis(cfg, data);

    % Average sub-segments by the first trialinfo column, as originally.
    trialGroups = unique(fractalSegments.trialinfo(:, 1));
    nGroups = numel(trialGroups);
    fractalByTrial = cell(1, nGroups);
    originalByTrial = cell(1, nGroups);
    for groupIdx = 1:nGroups
        cfg = [];
        cfg.trials = find(fractalSegments.trialinfo(:, 1) == trialGroups(groupIdx));
        cfg.avgoverrpt = 'yes';
        fractalByTrial{groupIdx} = ft_selectdata(cfg, fractalSegments);
        originalByTrial{groupIdx} = ft_selectdata(cfg, originalSegments);
    end

    fractal = ft_appendfreq([], fractalByTrial{:});
    original = ft_appendfreq([], originalByTrial{:});

    % Oscillatory residual = original spectrum - fractal spectrum.
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'x2-x1';
    oscillatory = ft_math(cfg, fractal, original);
end
