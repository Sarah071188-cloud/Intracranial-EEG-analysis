%% CNN - Step 2 - Implement deep learning

% PURPOSE: implement 1D-SCNN with k-fold 

% Input: xxx(subject)_xxx(session)_fooof_forcnn.mat (cnn_1_extractfeatures.m)
% Output: 1) result = classification performance for each fold
%         2) data_size = information on number of data points in each class for training, test, and, validation sets.
%         3) subj = subject ID
%         4) network = trained deep learning model for each fold

clear all

%% Set Directories

% Set main directory and add path
maindir = ('/Users/Data/');
addpath(genpath('/Users/Data/'))

% Set other directories
cnndir = [maindir 'stats/cnn/'];
featuresdir = [cnndir 'features/'];
dlm = '/';


%% Set up relevant variables

% List subject, session and dimension names
subjnames = {'s101','s102','s103','s104','s105','s106','s107'};
sessnames = {'d1','d2','d3','d4','d5','d6','d7'};
dimnames = {'off','free','stuck','goal','int','self','other'};

% Enter subjects, sessions and dimensions to process
subj_list = 1:7;
sess_list = 1:7;
dim_list = 1:7;
chan_list = 1:32;


%% Implement 1D-SCNN

for dim = dim_list
    
    for subj = subj_list
        
        % Load data
        filename = strcat(subjnames{subj},'_',dimnames{dim},'_fooof_forcnn.mat');
        load([featuresdir,filename]);
        if ~isempty(class1_data) && ~isempty(class2_data)
            class1_trials = class1_data; 
            class2_trials = class2_data; 
        elseif ~isempty(class1_data) && isempty(class2_data)
            class1_trials = class1_data;
            class2_trials = [];
        elseif isempty(class1_data) && ~isempty(class2_data)
            class1_trials = [];
            class2_trials = class2_data;
        end
        clear class*data
        
        
        %% Classification
        
        % Group data into a single variable and create labels
        data = [class1_trials;class2_trials];
        label = [zeros(size(class1_trials,1),1);ones(size(class2_trials,1),1)];
        class1 = 0;
        class2 = 1;
        
        % Implement classification only if number of data points is => threshold (50 data points).
        if size(class1_trials,1)>=50 && size(class2_trials,1)>=50
            % specify k value
            kfold = 5;
            label_predict = [];
            result = zeros(kfold,9);
            % k-fold cross-validation
            cv = cvpartition(label,'KFold',kfold,"Stratify",true);
            numFolds = cv.NumTestSets; % check the number of data points in each class for each fold
            nTestData = zeros(numFolds,size(categories(categorical(label)),1));

            for i = 1:numFolds
                
                % Obtain training and test class
                testclasses = (cv.test(i));
                trainclasses = ~testclasses;
                nCounts = countcats(categorical(label(testclasses))); % number of test set observations in each class
                nTestData(i,:) = nCounts';
                % Train data
                train_data = data(trainclasses,:,:);
                train_label= label(trainclasses,:);

                % Split into validation data
                cv1 = cvpartition(train_label,"Holdout",0.1,"Stratify",true);
                validation_class = (cv1.test);
                train_class =  ~validation_class;
                train_data_tmp = train_data;
                train_label_tmp = train_label;
                XValidation = train_data_tmp (validation_class,:,:);
                YValidation = train_label_tmp(validation_class,1);
                train_data = train_data_tmp (train_class,:,:);
                train_label = train_label_tmp(train_class,1);

                % Test data
                test_data = data(testclasses,:,:);
                test_label = label(testclasses,:);

                % Differentiate between original and balanced (see below) dataset
                train_data_ori = [train_data];
                train_label_ori = [train_label];

                
                %% Balancing Data
                
                for chan = chan_list
                    
                    dataset = array2table(squeeze(train_data_ori(:,chan,:)));
                    dataset = addvars(dataset, [string(train_label_ori)],'NewVariableNames','label');
                    
                    % determine which class has less data
                    if size(find(train_label_ori==class1),1) < size(find(train_label_ori==class2),1)
                        min_class = class1;
                    elseif size(find(train_label_ori==class1),1) == size(find(train_label_ori==class2),1)
                        min_class = [];
                    elseif size(find(train_label_ori==class1),1) > size(find(train_label_ori==class2),1)
                        min_class = class2;
                    end
                    % implement SMOTE (upsample data points in class with less data)
                    if ~isempty(min_class)
                        [newdata,visdata] = cnn_2a_smote((dataset), min_class,abs([size(find(train_label_ori==class1),1)-size(find(train_label_ori==class2),1)]), "NumNeighbors",8);
                        balanced_traindata(:,chan,:) = [table2array(newdata(1:abs([size(find(train_label_ori==class1),1)-size(find(train_label_ori==class2),1)]),1:size(newdata,2)-1));squeeze(train_data_ori(:,chan,:))];
                        balanced_trainclass = [str2double(table2array(newdata(1:abs([size(find(train_label_ori==class1),1)-size(find(train_label_ori==class2),1)]),size(newdata,2))));train_label_ori];
                    else % if classes were balanced, original = balanced dataset
                        balanced_traindata(:,chan,:) = [train_data_ori];
                        balanced_trainclass = [train_label_ori];
                    end
                end
                
                % Convert data into a deep learning data format
                for k=1:size(balanced_traindata,1)
                    train_data_all{k,:}=squeeze(balanced_traindata(k,:,:));
                end
                for k=1:size(XValidation,1)
                    valid_data_all{k,:}=squeeze(XValidation(k,:,:));
                end
                for k=1:size(test_data,1)
                    test_data_all{k,:}=squeeze(test_data(k,:,:));
                end
                val_lb = categorical(YValidation);
                
                % CNN Architecture
                filterSize = 7;
                numFilters = 25;
                layers_cnn= [ ...
                    sequenceInputLayer(size(train_data_all{1,1},1),Normalization="zscore")
                    convolution1dLayer(filterSize,numFilters,Padding="same",WeightsInitializer="glorot")
                    reluLayer
                    convolution1dLayer(5,20,Padding="same",Stride=1,WeightsInitializer="he")
                    reluLayer
                    globalAveragePooling1dLayer
                    dropoutLayer(0.1)
                    fullyConnectedLayer(2)
                    softmaxLayer
                    classificationLayer];

                options = trainingOptions("adam", ...
                    ValidationData={valid_data_all,val_lb}, ...
                    ValidationFrequency=50, ...
                    ExecutionEnvironment="cpu", ...
                    GradientThreshold=1, ...
                    MaxEpochs=100, ...
                    MiniBatchSize=100, ...
                    Shuffle='every-epoch', ...
                    OutputNetwork='best-validation-loss',... 
                    Verbose=true);
                
                % Train CNN
                dt = train_data_all;
                lb = categorical(balanced_trainclass);
                net = trainNetwork(dt,lb,layers_cnn,options);
                % Testing the network
                [testPred2,scores2] = classify(net,test_data_all);
                % Confusion matrix
                confmat2 = confusionmat(double(categorical(test_label)),double(testPred2));
                % AUC computation
                [rx,ry,~,result(i,9)] = perfcurve(double(categorical(test_label)),scores2(:,2),2);
                % Classifier Performance
                [result(i,1:8)] = cnn_2b_compute_metrics(confmat2(2,2),confmat2(1,1),confmat2(1,2),confmat2(2,1));

                % Final verification of validation, test and train set size
                % and distribution of data points between two classes
                data_size.train(i,1)                    = size(train_data_ori,1);
                data_size.test(i,1)                     = size(test_data,1);
                data_size.train_class1_size(i,1)        = size(find(balanced_trainclass==0),1);
                data_size.train_class2_size(i,1)        = size(find(balanced_trainclass==1),1);
                data_size.validation(i,1)               = size(XValidation,1);
                data_size.validation_class1_size(i,1)   = size(find(YValidation==0),1);
                data_size.validation_class2_size(i,1)   = size(find(YValidation==1),1);
                data_size.test_class1_size(i,1)         = size(find(test_label==0),1);
                data_size.test_class2_size(i,1)         = size(find(test_label==1),1);

                data_size.confmat2{i,1} = confmat2;
                
                % Saving the trained network for all the k folds
                network{i,1} = net;
                clear train_label* train_class*i train_data* dt lb val* balanced* confmat2 testPred2 scores2
                
            end
        end
        
        % Save results
        filename = strcat('best_classification_cnn_model_',subjnames{subj},'_',dimnames{dim});
        save([cnndir filename],"data_size","subj","result","network")
        
    end
end


