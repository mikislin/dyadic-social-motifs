function scaleBankSec = make_logscale_chunk_bank(varargin)
%MAKE_LOGSCALE_CHUNK_BANK Generate a dense log-spaced chunk-scale bank.
%
% Example
%   scaleBankSec = make_logscale_chunk_bank('minSec',0.2,'maxSec',8.0,'nScales',25)

p = inputParser;
p.addParameter('minSec', 0.2, @(x)isscalar(x) && x > 0);
p.addParameter('maxSec', 8.0, @(x)isscalar(x) && x > 0);
p.addParameter('nScales', 25, @(x)isscalar(x) && x >= 3);
p.addParameter('roundDigits', 4, @(x)isscalar(x) && x >= 0);
p.parse(varargin{:});
P = p.Results;

assert(P.maxSec > P.minSec, 'maxSec must exceed minSec.');
scaleBankSec = logspace(log10(P.minSec), log10(P.maxSec), P.nScales);
scaleBankSec = unique(round(scaleBankSec, P.roundDigits));
end
