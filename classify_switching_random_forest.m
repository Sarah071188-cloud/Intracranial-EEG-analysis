%% Within-subject random forest classification of switching and attention
% Requires Statistics and Machine Learning Toolbox and per_metric.m.
% IMPORTANT: Feature ranks below are externally fixed. If learned from
% these same observations, use nested training-only feature selection.
% This refactor replaces column-wise SMOTE with TRAINING-ONLY random
% oversampling of whole rows. It is NOT numerically identical to SMOTE.
clearvars;

subjectIDs = arrayfun(@num2str, 101:120, 'UniformOutput', false);
siteNames = [repmat({'Calgary'},1,14), {'CPMC'}, repmat({'WashU'},1,5)];
networkNames = {'VAN','DMN','DAN','FPCN'};
frequencyNames = {'Theta','Alpha','Beta','HGamma','LGamma'};
class1Names = {'Switch to external','Switch to external','Switch to external', ...
               'Switch to internal','Switch to internal'};
class2Names = {'Switch to internal','Post switch external attention', ...
               'Pre switch internal attention','Post switch internal attention', ...
               'Pre switch external attention'};
featureRanks = { ...
    [20 15 5 13 9 2 7 11 16 6 18 17 4 10 12 19 3 8 14 1], ... % Switch to ext vs switch to int]
    [7 20 11 16 14 17 9 19 2 6 18 10 13 3 12 15 1 5 4 8], ... % Extswitch vs post-extswitch
    [15 12 8 13 20 19 18 6 7 3 10 9 16 1 17 5 11 14 2 4], ... % Extswitch vs pr-intswitch
    [3 10 11 19 16 1 7 15 9 14 20 18 12 4 6 2 5 13 8 17], ... % Intswitch vs post-int attn
    [15 7 11 6 16 13 14 18 10 2 1 3 19 9 4 20 5 12 17 8]}; % Intswitch vs pre-ext attn

minimumTrials = 5;
numberOfFeatures = 10;
numberOfFolds = 5;
numberOfIterations = 10;
numberOfTrees = 50;
resultsDirectory = fullfile('Analysis','ML_results','VAN_FPCN_DMN_DAN');
if ~exist(resultsDirectory,'dir'), mkdir(resultsDirectory); end
rng(1,'twister'); % Reproducible random partitions and ensembles.

% Assemble each subject once rather than loading the same files 50 times.
subjectFeatures = cell(numel(subjectIDs),1);
for subjectIdx = 1:numel(subjectIDs)
    subjectID = subjectIDs{subjectIdx};
    subjectFolder = fullfile('Analysis',[siteNames{subjectIdx},subjectID],'Switching');
    stageNames = {'Switching','Postvideo','Prevideo'};
    stageFolders = {'Analysis_3sec_updated_blockbootstrap_200rand_10BS', ...
                    'Analysis_5sec_updated_blockbootstrap_200rand_10BS', ...
                    'Analysis_5sec_updated_blockbootstrap_200rand_10BS'};
    subjectData = struct();
    for stageIdx = 1:3
        stageName = stageNames{stageIdx};
        stageFeatures = cell(1,numel(frequencyNames));
        validStage = true;
        for frequencyIdx = 1:numel(frequencyNames)
            frequencyName = frequencyNames{frequencyIdx};
            % Legacy source files may use theta/Theta or Prevido/Prevideo.
            stem = ['S',subjectID,'_ML_feat_',stageName,'_',frequencyName];
            candidateFiles = {fullfile(subjectFolder,stageFolders{stageIdx},[stem,'.mat']), ...
                fullfile(subjectFolder,stageFolders{stageIdx}, ...
                ['S',subjectID,'_ML_feat_',stageName,'_',lower(frequencyName),'.mat'])};
            if strcmp(stageName,'Postvideo') || strcmp(stageName,'Prevideo')
                misspelledStage = strrep(stageName,'video','vido');
                candidateFiles{end+1} = fullfile(subjectFolder,stageFolders{stageIdx}, ...
                    ['S',subjectID,'_ML_feat_',misspelledStage,'_',frequencyName,'.mat']); %#ok<SAGROW>
                candidateFiles{end+1} = fullfile(subjectFolder,stageFolders{stageIdx}, ...
                    ['S',subjectID,'_ML_feat_',misspelledStage,'_',lower(frequencyName),'.mat']); %#ok<SAGROW>
            end
            fileIdx = find(cellfun(@(p) exist(p,'file')==2,candidateFiles),1);
            if isempty(fileIdx)
                warning('Missing feature file: %s (%s, %s)',subjectID,stageName,frequencyName);
                validStage = false;
                break;
            end
            loaded = load(candidateFiles{fileIdx},'Final');
            if ~isfield(loaded,'Final')
                validStage = false;
                break;
            end
            stageFeatures{frequencyIdx} = loaded.Final;
        end
        if ~validStage
            subjectData.(stageName) = [];
            continue;
        end
        stageData = struct('c1',[],'c2',[]);
        for frequencyIdx = 1:numel(frequencyNames)
            features = stageFeatures{frequencyIdx};
            for networkIdx = 1:numel(networkNames)
                networkName = networkNames{networkIdx};
                if ~isfield(features,networkName) || isempty(features.(networkName)) || ...
                   ~isfield(features.(networkName),'all_feat')
                    validStage = false;
                    break;
                end
                stageData.c1 = [stageData.c1, features.(networkName).all_feat.c1(:)]; %#ok<AGROW>
                stageData.c2 = [stageData.c2, features.(networkName).all_feat.c2(:)]; %#ok<AGROW>
            end
            if ~validStage, break; end
        end
        if validStage && size(stageData.c1,2)==20 && size(stageData.c2,2)==20
            subjectData.(stageName) = stageData;
        else
            subjectData.(stageName) = [];
        end
    end
    subjectFeatures{subjectIdx} = subjectData;
end

% Results: 16 legacy columns (subject, 8 metrics, AUC, kappa,
% training class 1/0 counts, test class 1/0 counts, balanced accuracy).
for comparisonIdx = 1:5
    allResultIterMean = nan(numberOfIterations,16);
    allResultIterStd = nan(numberOfIterations,16);
    comparisonSubjectResults = cell(numberOfIterations,numel(subjectIDs));
    for iterationIdx = 1:numberOfIterations
        subjectMeans = nan(numel(subjectIDs),16);
        subjectStds = nan(numel(subjectIDs),16);
        for subjectIdx = 1:numel(subjectIDs)
            subjectData = subjectFeatures{subjectIdx};
            if isempty(subjectData.Switching), continue; end
            switch comparisonIdx
                case 1
                    firstClass = subjectData.Switching.c1;
                    secondClass = subjectData.Switching.c2;
                case 2
                    if isempty(subjectData.Postvideo), continue; end
                    firstClass = subjectData.Switching.c1;
                    secondClass = subjectData.Postvideo.c1;
                case 3
                    if isempty(subjectData.Prevideo), continue; end
                    firstClass = subjectData.Switching.c1;
                    secondClass = subjectData.Prevideo.c2;
                case 4
                    if isempty(subjectData.Postvideo), continue; end
                    firstClass = subjectData.Switching.c2;
                    secondClass = subjectData.Postvideo.c2;
                case 5
                    if isempty(subjectData.Prevideo), continue; end
                    firstClass = subjectData.Switching.c2;
                    secondClass = subjectData.Prevideo.c1;
            end
            if size(firstClass,1)<max(minimumTrials,numberOfFolds) || ...
               size(secondClass,1)<max(minimumTrials,numberOfFolds)
                continue;
            end
            selectedFeatures = featureRanks{comparisonIdx}(1:numberOfFeatures);
            predictorData = [firstClass(:,selectedFeatures);secondClass(:,selectedFeatures)];
            labels = [ones(size(firstClass,1),1);zeros(size(secondClass,1),1)];
            if any(~isfinite(predictorData(:)))
                warning('Skipping %s: NaN/Inf in features.',subjectIDs{subjectIdx});
                continue;
            end
            partition = cvpartition(labels,'KFold',numberOfFolds,'Stratify',true);
            foldResults = nan(numberOfFolds,16);
            for foldIdx = 1:numberOfFolds
                trainingRows = training(partition,foldIdx);
                testingRows = test(partition,foldIdx);
                trainData = predictorData(trainingRows,:);
                trainLabels = labels(trainingRows);
                testData = predictorData(testingRows,:);
                testLabels = labels(testingRows);

                % Whole-row random oversampling of minority training class.
                countPositive = sum(trainLabels==1);
                countNegative = sum(trainLabels==0);
                if countPositive~=countNegative
                    minorityLabel = double(countPositive < countNegative);
                    minorityRows = find(trainLabels==minorityLabel);
                    difference = abs(countPositive-countNegative);
                    replicatedRows = minorityRows(randi(numel(minorityRows),difference,1));
                    trainData = [trainData;trainData(replicatedRows,:)];
                    trainLabels = [trainLabels;trainLabels(replicatedRows)];
                end
                model = fitcensemble(trainData,trainLabels,'Method','Bag', ...
                    'Learners','tree','NumLearningCycles',numberOfTrees);
                [predictedLabels,scores] = predict(model,testData);
                confusion = confusionmat(testLabels,predictedLabels,'Order',[0 1]);
                tn = confusion(1,1); fp = confusion(1,2);
                fn = confusion(2,1); tp = confusion(2,2);
                metrics = per_metric(tp,tn,fp,fn);
                positiveColumn = find(model.ClassNames==1,1);
                [~,~,~,auc] = perfcurve(testLabels,scores(:,positiveColumn),1);
                observedAgreement = (tp+tn)/sum(confusion(:));
                expectedAgreement = ((tp+fp)*(tp+fn)+(tn+fn)*(tn+fp))/sum(confusion(:))^2;
                if expectedAgreement==1
                    kappa = 0;
                else
                    kappa = (observedAgreement-expectedAgreement)/(1-expectedAgreement);
                end
                foldResults(foldIdx,:) = [subjectIdx,metrics,auc,kappa, ...
                    sum(trainLabels==1),sum(trainLabels==0), ...
                    sum(testLabels==1),sum(testLabels==0),mean(metrics(2:3))];
            end
            subjectMeans(subjectIdx,:) = mean(foldResults,1,'omitnan');
            subjectStds(subjectIdx,:) = std(foldResults,0,1,'omitnan');
            comparisonSubjectResults{iterationIdx,subjectIdx} = foldResults;
        end
        allResultIterMean(iterationIdx,:) = mean(subjectMeans,1,'omitnan');
        allResultIterStd(iterationIdx,:) = std(subjectStds,0,1,'omitnan');
    end
    class1Name = class1Names{comparisonIdx};
    class2Name = class2Names{comparisonIdx};
    resultFile = fullfile(resultsDirectory,sprintf( ...
        'RF_results_across_iterations_k_%d_lim_dp_%d_comparison_%d.mat', ...
        numberOfFolds,minimumTrials,comparisonIdx));
    save(resultFile,'allResultIterMean','allResultIterStd', ...
        'comparisonSubjectResults','class1Name','class2Name', ...
        'featureRanks','subjectIDs');
    fprintf('Completed comparison %d: %s vs %s\n',comparisonIdx,class1Name,class2Name);
end
