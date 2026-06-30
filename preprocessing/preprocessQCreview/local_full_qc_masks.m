function masks = local_full_qc_masks(qa, nFrames, nNodes)
masks.lowConf = local_mask_field(qa, 'lowConfMask', nFrames, nNodes);
masks.jump = local_mask_field(qa, 'jumpMask', nFrames, nNodes);
masks.interp = local_mask_field(qa, 'interpMask', nFrames, nNodes);
masks.geom = local_mask_field(qa, 'geomMask', nFrames, nNodes);
masks.arena = local_mask_field(qa, 'arenaMask', nFrames, nNodes);
masks.finalNan = local_mask_field(qa, 'finalNanMask', nFrames, nNodes);
end