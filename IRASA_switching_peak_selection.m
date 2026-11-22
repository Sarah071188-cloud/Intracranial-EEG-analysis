%% IRASA peak selection | Attention switching
% INPUT:  Subject-specific IRASA output from the switching analysis.
% OUTPUT: osci_elec, retaining the original fields and selection criteria.
%
% Required functions: findpeaks (Signal Processing Toolbox).
% The IRASA output must contain oscillatory_class_all and
% oscillatory_class_all_bs.
%
% Peak inclusion: positive bootstrap peak, prominence >= 1.64, and
% electrode-wide maximum oscillatory power < 55.
% For each band, retain the frequency of the most prominent eligible peak.
% No eligible peak is represented by NaN.

close all;

subjectIDs = arrayfun(@num2str, 101:120, 'UniformOutput', false);
siteNames = [repmat({'Calgary'}, 1, 14), {'CPMC'}, repmat({'WashU'}, 1, 5)];

prominenceThreshold = 1.64;
maximumOscillatoryPower = 55;
frequencyBands = struct( ...
    'theta', [4, 8], ...   % 4 <= f < 8 Hz
    'alpha', [8, 12], ...  % 8 <= f <= 12 Hz
    'beta',  [13, 30]);   % 13 <= f <= 30 Hz

for subjectIdx = 1:numel(subjectIDs)
    subjectID = subjectIDs{subjectIdx};
    siteName = siteNames{subjectIdx};
    subjectPrefix = ['S', subjectID];

    outputDirectory = fullfile('Analysis', [siteName, subjectID], ...
        'Switching', 'Analysis_3sec_updated_blockbootstrap_200rand_10BS');
    irasaFile = fullfile(outputDirectory, ...
        [subjectPrefix, '_Switching_IRASA_output.mat']);

    % Load only the required IRASA variables.
    irasaResults = load(irasaFile, ...
        'oscillatory_class_all', 'oscillatory_class_all_bs');
    oscillatorySpectrum = irasaResults.oscillatory_class_all;
    bootstrapSpectrum = irasaResults.oscillatory_class_all_bs;

    nChannels = numel(oscillatorySpectrum.label);
    frequencies = oscillatorySpectrum.freq;

    assert(size(bootstrapSpectrum, 1) == nChannels, ...
        'Bootstrap spectrum channel count does not match labels for subject %s.', subjectID);
    assert(size(bootstrapSpectrum, 2) == numel(frequencies), ...
        'Bootstrap spectrum frequency count does not match frequency vector for subject %s.', subjectID);
    assert(size(oscillatorySpectrum.powspctrm, 1) == nChannels, ...
        'Oscillatory spectrum channel count does not match labels for subject %s.', subjectID);

    % Preserve the legacy output structure and its field names.
    osci_elec = struct();
    osci_elec.pk.pr = cell(nChannels, 1);
    osci_elec.pk.fr = cell(nChannels, 1);
    osci_elec.pk.wd = cell(nChannels, 1);
    osci_elec.pk.pk = cell(nChannels, 1);
    osci_elec.theta = nan(nChannels, 1);
    osci_elec.alpha = nan(nChannels, 1);
    osci_elec.beta = nan(nChannels, 1);

    for channelIdx = 1:nChannels
        [peakValues, peakFrequencies, peakWidths, peakProminences] = ...
            findpeaks(bootstrapSpectrum(channelIdx, :), frequencies);

        % Original rule: reject all peaks for an electrode if its maximum
        % oscillatory-spectrum power is 55 or greater.
        electrodeMaxPower = max(oscillatorySpectrum.powspctrm(channelIdx, :));
        eligiblePeaks = (peakValues > 0) & ...
            (electrodeMaxPower < maximumOscillatoryPower) & ...
            (peakProminences >= prominenceThreshold);

        peakValues = peakValues(eligiblePeaks);
        peakFrequencies = peakFrequencies(eligiblePeaks);
        peakWidths = peakWidths(eligiblePeaks);
        peakProminences = peakProminences(eligiblePeaks);

        % Match the original sorting: prominence descending, then frequency
        % ascending. Retain the corresponding values in every peak field.
        [~, prominenceOrder] = sort(peakProminences, 'descend');
        peakValues = peakValues(prominenceOrder);
        peakFrequencies = peakFrequencies(prominenceOrder);
        peakWidths = peakWidths(prominenceOrder);
        peakProminences = peakProminences(prominenceOrder);

        [~, frequencyOrder] = sort(peakFrequencies, 'ascend');
        peakValues = peakValues(frequencyOrder);
        peakFrequencies = peakFrequencies(frequencyOrder);
        peakWidths = peakWidths(frequencyOrder);
        peakProminences = peakProminences(frequencyOrder);

        osci_elec.pk.pr{channelIdx} = peakProminences;
        osci_elec.pk.fr{channelIdx} = peakFrequencies;
        osci_elec.pk.wd{channelIdx} = peakWidths;
        osci_elec.pk.pk{channelIdx} = peakValues;

        osci_elec.theta(channelIdx) = select_band_peak( ...
            peakFrequencies, peakProminences, frequencyBands.theta, false);
        osci_elec.alpha(channelIdx) = select_band_peak( ...
            peakFrequencies, peakProminences, frequencyBands.alpha, true);
        osci_elec.beta(channelIdx) = select_band_peak( ...
            peakFrequencies, peakProminences, frequencyBands.beta, true);
    end

    saveFile = fullfile(outputDirectory, ...
        [subjectPrefix, '_IRASA_peak_electrode_selection_switching.mat']);
    save(saveFile, 'osci_elec');
    fprintf('IRASA peak selection completed: %s (%d channels)\n', ...
        subjectPrefix, nChannels);
end

%% Local function: choose the most prominent peak in one frequency band
function selectedFrequency = select_band_peak(peakFrequencies, peakProminences, bandLimits, includeUpperLimit)
    selectedFrequency = NaN;

    if includeUpperLimit
        inBand = peakFrequencies >= bandLimits(1) & ...
            peakFrequencies <= bandLimits(2);
    else
        inBand = peakFrequencies >= bandLimits(1) & ...
            peakFrequencies < bandLimits(2);
    end

    bandFrequencies = peakFrequencies(inBand);
    bandProminences = peakProminences(inBand);

    if isempty(bandFrequencies)
        return;
    end

    % The original script chooses the first maximum prominence after
    % frequency sorting; max preserves that behavior for equal prominences.
    [~, highestProminenceIdx] = max(bandProminences);
    selectedFrequency = bandFrequencies(highestProminenceIdx);
end
