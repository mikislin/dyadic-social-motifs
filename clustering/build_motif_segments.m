function Segments = build_motif_segments(anchorTable, labels, varargin)
%BUILD_MOTIF_SEGMENTS Form anchor-level motif segments within each session.

p = inputParser;
p.addParameter('MergeGapFrames', 0, @(x)isscalar(x) && x >= 0);
p.parse(varargin{:});
P = p.Results;

labels = labels(:);
assert(height(anchorTable) == numel(labels), 'anchorTable / labels size mismatch.');

[~, ord] = sortrows([anchorTable.session_index, anchorTable.anchor_frame], [1 2]);
A = anchorTable(ord,:);
z = labels(ord);

rows = struct([]);
r = 0;
sessions = unique(A.session_index, 'stable');
for s = 1:numel(sessions)
    idx = find(A.session_index == sessions(s));
    if isempty(idx)
        continue
    end
    fr = A.anchor_frame(idx);
    lab = z(idx);
    startLocal = 1;
    for i = 2:numel(idx)
        breakHere = (lab(i) ~= lab(i-1)) || ((fr(i) - fr(i-1)) > (P.MergeGapFrames + 1));
        if breakHere
            stopLocal = i - 1;
            r = r + 1;
            rows(r).session_index = sessions(s); %#ok<AGROW>
            rows(r).cluster = lab(startLocal);
            rows(r).start_frame = fr(startLocal);
            rows(r).stop_frame = fr(stopLocal);
            rows(r).length_frames = fr(stopLocal) - fr(startLocal) + 1;
            rows(r).n_anchors = stopLocal - startLocal + 1;
            startLocal = i;
        end
    end
    stopLocal = numel(idx);
    r = r + 1;
    rows(r).session_index = sessions(s); %#ok<AGROW>
    rows(r).cluster = lab(startLocal);
    rows(r).start_frame = fr(startLocal);
    rows(r).stop_frame = fr(stopLocal);
    rows(r).length_frames = fr(stopLocal) - fr(startLocal) + 1;
    rows(r).n_anchors = stopLocal - startLocal + 1;
end

Segments = struct();
Segments.table = struct2table(rows);
Segments.sortedAnchorTable = A;
Segments.sortedLabels = z;
end
