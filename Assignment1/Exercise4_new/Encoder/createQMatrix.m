
function adaptiveQP = getAdaptiveQP(block, baseQP)
    blockVariance = var(double(block(:)));
    if blockVariance > 1000
        adaptiveQP = max(0, baseQP - 2);  % Reduce QP for high-variance blocks
    elseif blockVariance < 100
        adaptiveQP = min(51, baseQP + 2);  % Increase QP for low-variance blocks
    else
        adaptiveQP = baseQP;
    end
    end