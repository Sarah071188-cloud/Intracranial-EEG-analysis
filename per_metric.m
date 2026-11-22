function resultCase = per_metric(tp, tn, fp, fn)
%PER_METRIC Eight binary classification metrics in legacy output order.
% [accuracy, sensitivity, specificity, PPV, NPV, MCC, F1, GMean]
% Positive class is the class represented by tp/fn.
    safeDivide = @(numerator, denominator) numerator ./ max(denominator, eps);
    accuracy = safeDivide(tp + tn, tp + tn + fp + fn);
    sensitivity = safeDivide(tp, tp + fn);
    specificity = safeDivide(tn, tn + fp);
    precision = safeDivide(tp, tp + fp);
    negativePredictiveValue = safeDivide(tn, tn + fn);
    mccNumerator = tp*tn - fp*fn;
    mccDenominator = sqrt((tp+fn)*(tp+fp)*(tn+fn)*(tn+fp));
    if mccDenominator == 0
        mcc = 0;
    else
        mcc = mccNumerator / mccDenominator;
    end
    f1 = safeDivide(2*precision*sensitivity, precision+sensitivity);
    geometricMean = sqrt(sensitivity*specificity);
    resultCase = [accuracy,sensitivity,specificity,precision,negativePredictiveValue,mcc,f1,geometricMean];
end
