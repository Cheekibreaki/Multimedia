function [originalBlock] = invdifferential(diffedBlock,currBlock)
    originalBlock = diffedBlock + currBlock;
end