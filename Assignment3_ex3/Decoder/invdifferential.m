function [originalBlock] = invdifferential(diffedBlock,currBlock)
    originalBlock = int32(diffedBlock) + int32(currBlock);
end