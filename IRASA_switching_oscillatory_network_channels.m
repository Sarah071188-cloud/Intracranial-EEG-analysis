%% IRASA oscillatory electrodes by network | Switching task
% Identify HFA-responsive electrodes with an IRASA peak in theta/alpha/beta.
% Retains the legacy MAT-file structures and summary-column order.
% Dependencies: output from HFA switching and IRASA peak-selection scripts.

clearvars;
close all;

subjectIDs = arrayfun(@num2str, 101:120, 'UniformOutput', false);
siteNames = [repmat({'Calgary'}, 1, 14), {'CPMC'}, repmat({'WashU'}, 1, 5)];
saveResults = true;

% Preserve the network order and the 31-column legacy summary layout.
networkNames = {'FPCNA', 'FPCNB', 'VAN', 'DAN', 'DMN', 'FPCN'};
frequencyBands = {'theta', 'alpha', 'beta'};
summaryHeaders = cell(1, 1 + numel(networkNames) * 5);
summaryHeaders{1} = 'Subject';
for networkIdx = 1:numel(networkNames)
    networkName = networkNames{networkIdx};
    firstCol = 2 + (networkIdx - 1) * 5;
    summaryHeaders(firstCol:firstCol+4) = { ...
        [networkName ' elec'], [networkName ' active elec'], ...
        [networkName ' theta elec'], [networkName ' alpha elec'], ...
        [networkName ' beta elec']};
end
final_out = cell(numel(subjectIDs) + 1, numel(summaryHeaders));
final_out(1,:) = summaryHeaders;

% Retain all 12 network groupings in the saved output, including subnets.
allNetworkNames = {'VAN', 'FPCN', 'DMN', 'DMNA', 'DMNB', 'DMNC', ...
                   'DMND', 'DAN', 'FPCNA', 'FPCNB'};

for subjectIdx = 1:numel(subjectIDs)
    subjectID = subjectIDs{subjectIdx};
    subjectName = ['S' subjectID];
    outputDirectory = fullfile('Analysis', [siteNames{subjectIdx} subjectID], ...
        'Switching', 'Analysis_3sec_updated_blockbootstrap_200rand_10BS');

    % Load the legacy HFA results: post-rejection trial selections and
    % pre-rejection channel/network indices (as in the original script).
    trialFile = fullfile(outputDirectory, ...
        [subjectName '_Switching_endpoint200schan_responsive_Yeo_label_meanoutlierremoved.mat']);
    networkFile = fullfile(outputDirectory, ...
        [subjectName '_Switching_endpoint200schan_responsive_Yeo_label_check.mat']);
    irasaFile = fullfile(outputDirectory, [subjectName '_Switching_IRASA_output.mat']);
    peakFile = fullfile(outputDirectory, ...
        [subjectName '_IRASA_peak_electrode_selection_switching.mat']);

    % The original script loads the trial indices, but never uses them.
    % Avoid loading those large/unneeded variables into the workspace.
    hfa = load(networkFile, 'ch_resp_ind', 'chan_resp', ...
        'HFA_avg_data1', 'HFA_avg_data2', ...
        'VAN_ind_Ori', 'FPCN_ind_Ori', 'DMN_ind_Ori', ...
        'DMNA_ind_Ori', 'DMNB_ind_Ori', 'DMNC_ind_Ori', 'DMND_ind_Ori', ...
        'DAN_ind_Ori', 'FPCNA_ind_Ori', 'FPCNB_ind_Ori');
    irasa = load(irasaFile, 'oscillatory_class_all');
    peaks = load(peakFile, 'osci_elec');

    % The post-rejection file was loaded by the original but its selected
    % trial indices were not referenced. Check that it exists for parity.
    assert(isfile(trialFile), 'Missing HFA result: %s', trialFile);

    out_all_actv_elec = struct();
    out_osci_elec = struct();
    for networkIdx = 1:numel(allNetworkNames)
        networkName = allNetworkNames{networkIdx};
        originalIndexField = [networkName '_ind_Ori'];
        out_all_actv_elec.(originalIndexField) = hfa.(originalIndexField);

        for bandIdx = 1:numel(frequencyBands)
            bandName = frequencyBands{bandIdx};
            fieldName = [networkName '_' bandName];
            out_osci_elec.(fieldName) = struct('out', [], ...
                'name', {cell(0,1)}, 'osci_elec', []);

            selectedNetworkIndices = hfa.(originalIndexField);
            bandPeakFrequencies = peaks.osci_elec.(bandName);
            channelLabels = irasa.oscillatory_class_all.label;

            % Preserve the original matching order: traverse IRASA channels,
            % keep only finite peak frequencies and exact label matches,
            % then restrict to the specified responsive-network indices.
            selectedIRASAChannels = zeros(numel(bandPeakFrequencies), 1);
            selectedResponsiveIndices = zeros(numel(bandPeakFrequencies), 1);
            selectedCount = 0;
            for channelIdx = 1:numel(bandPeakFrequencies)
                if isnan(bandPeakFrequencies(channelIdx))
                    continue;
                end
                responsiveIndex = find(strcmp(channelLabels{channelIdx}, ...
                    hfa.chan_resp));
                if numel(responsiveIndex) ~= 1
                    continue;
                end
                if ~ismember(responsiveIndex, selectedNetworkIndices)
                    continue;
                end
                selectedCount = selectedCount + 1;
                selectedIRASAChannels(selectedCount) = channelIdx;
                selectedResponsiveIndices(selectedCount) = responsiveIndex;
            end

            selectedIRASAChannels = selectedIRASAChannels(1:selectedCount);
            selectedResponsiveIndices = selectedResponsiveIndices(1:selectedCount);
            out_osci_elec.(fieldName).out = selectedResponsiveIndices;
            out_osci_elec.(fieldName).name = ...
                reshape(hfa.chan_resp(selectedResponsiveIndices), [], 1);
            out_osci_elec.(fieldName).osci_elec = selectedIRASAChannels;
        end
    end

    % Preserve legacy peak vectors and MAT-file field names.
    out_osci_elec.thetapeak = peaks.osci_elec.theta;
    out_osci_elec.alphapeak = peaks.osci_elec.alpha;
    out_osci_elec.betapeak = peaks.osci_elec.beta;

    if saveResults
        save(fullfile(outputDirectory, ...
            [subjectName '_IRASA_active_electrodes_network_meanoutlierremoved_Switching.mat']), ...
            'out_osci_elec', 'out_all_actv_elec');
    end

    % The original first column of each block counts all HFA channels,
    % not the total number of channels in the named network. Preserve it.
    final_out{subjectIdx+1,1} = subjectName;
    for networkIdx = 1:numel(networkNames)
        networkName = networkNames{networkIdx};
        firstCol = 2 + (networkIdx-1)*5;
        if ismember(networkName, {'DMN', 'FPCN'})
            final_out{subjectIdx+1,firstCol} = size(hfa.HFA_avg_data2,2);
        else
            final_out{subjectIdx+1,firstCol} = size(hfa.HFA_avg_data1,2);
        end
        final_out{subjectIdx+1,firstCol+1} = ...
            size(out_all_actv_elec.([networkName '_ind_Ori']), 1);
        for bandIdx = 1:numel(frequencyBands)
            fieldName = [networkName '_' frequencyBands{bandIdx}];
            final_out{subjectIdx+1,firstCol+1+bandIdx} = ...
                size(out_osci_elec.(fieldName).out, 1);
        end
    end
    fprintf('Processed %s (%s)\n', subjectName, siteNames{subjectIdx});
end
