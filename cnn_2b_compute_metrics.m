function [metrics] = cnn_2b_compute_metrics(tps,tns,fps,fns)

% PURPOSE: calculate classification performance
% part of cnn_2_kfold

% Function Inputs: true/false positive/negative
% Output: classification performance metrics

% Accuracy and Confusion Matrix values
sas = ((tps+tns)/(tps+tns+fps+fns)); % Accuracy
sspes = (tns/(tns+fps)); % Sensitivity
ssens = (tps/(tps+fns)); % Specificity
sppv = (tps/(tps+fps)); % PPV
snpv = (tns/(tns+fns)); % NPV

% MCC 
mcc_numer = ((tps*tns)-(fps*fns));
mcc_denom = sqrt((tps+fns)*(tps+fps)*(tns+fns)*(tns+fps));
mcc = (mcc_numer/mcc_denom);

 % F1-Score 
F1_Score_numer = 2*(ssens)*(sppv);
F1_Score_denom = (ssens)+(sppv);
F1_Score = (F1_Score_numer/F1_Score_denom); % F1_score
Gmean = (sqrt((ssens)*(sspes))); % G-mean

metrics = [sas,ssens,sspes,sppv,snpv,mcc,F1_Score,Gmean];

end