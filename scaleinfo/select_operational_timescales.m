function Selected = select_operational_timescales(ScaleScore, varargin)
%SELECT_OPERATIONAL_TIMESCALES Choose a hierarchical subset of useful scales.
%
% Input
%   ScaleScore.scaleTable : output table from score_multiscale_chunk_bank
%
% Name-value pairs
%   'nMicro'   : number of micro scales to retain (default 3)
%   'nMotif'   : number of motif scales to retain (default 4)
%   'nContext' : number of context scales to retain (default 2)
%   'minLogGap': minimum separation between retained scales within band (default 0.12)
%
% Output
%   Selected   : table of selected scales and hierarchical assignments

p = inputParser;
p.addParameter('nMicro', 3, @(x)isscalar(x) && x >= 1);
p.addParameter('nMotif', 4, @(x)isscalar(x) && x >= 1);
p.addParameter('nContext', 2, @(x)isscalar(x) && x >= 1);
p.addParameter('minLogGap', 0.12, @(x)isscalar(x) && x >= 0);
p.parse(varargin{:});
P = p.Results;

T = ScaleScore.scaleTable;
assert(istable(T) && height(T) >= 3, 'ScaleScore.scaleTable must be a non-empty table.');

bands = ["micro","motif","context"];
want = [P.nMicro, P.nMotif, P.nContext];
out = table();
for b = 1:numel(bands)
    Tb = T(T.initial_band == bands(b), :);
    if isempty(Tb)
        continue
    end
    switch bands(b)
        case "micro"
            scoreName = 'composite_micro';
        case "motif"
            scoreName = 'composite_motif';
        otherwise
            scoreName = 'composite_context';
    end
    Tb = sortrows(Tb, scoreName, 'descend');
    keep = false(height(Tb),1);
    chosenLog = [];
    for i = 1:height(Tb)
        ls = log10(Tb.chunk_sec(i));
        if isempty(chosenLog) || all(abs(ls - chosenLog) >= P.minLogGap)
            keep(i) = true;
            chosenLog(end+1,1) = ls; %#ok<AGROW>
        end
        if sum(keep) >= want(b)
            break
        end
    end
    Tk = Tb(keep,:);
    Tk.hierarchical_role = repmat(bands(b), height(Tk), 1);
    Tk.rank_within_role = (1:height(Tk))';
    out = [out; Tk]; %#ok<AGROW>
end

% Add a recommended processing order.
out.processing_priority = nan(height(out),1);
for i = 1:height(out)
    switch string(out.hierarchical_role(i))
        case "micro"
            out.processing_priority(i) = 1;
        case "motif"
            out.processing_priority(i) = 2;
        case "context"
            out.processing_priority(i) = 3;
    end
end
out = sortrows(out, {'processing_priority','rank_within_role','chunk_sec'});
Selected = out;
end
