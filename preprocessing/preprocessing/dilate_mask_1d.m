function maskOut = dilate_mask_1d(maskIn, halfWidth)
%DILATE_MASK_1D Expand a logical vector by a half-width in frames.

maskIn = logical(maskIn(:));
if halfWidth <= 0
    maskOut = maskIn;
    return
end
kernel = ones(2*halfWidth + 1, 1);
maskOut = conv(double(maskIn), kernel, 'same') > 0;
maskOut = logical(maskOut);
end
