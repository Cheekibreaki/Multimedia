function diffEncoded = differentialEncode(modes)
    diffEncoded = zeros(size(modes));
    diffEncoded(1) = modes(1);
    for i = 2:numel(modes)
        diffEncoded(i) = modes(i) - modes(i-1);
    end
end