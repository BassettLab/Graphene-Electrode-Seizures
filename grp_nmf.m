function [W,H,FMT,fid,fset,trange,fastrange,wsort] = grp_nmf(EPH, IMG, fcE, Fs, whichk, goodk, whichset)

% This function puts together a feature matrix for running the non-negative
% matrix factorisation, and then runs it either across a range of k-values
% and/or for a specified 'optimum' k
%
%   whichk - toggles scanning through different values for the parameter
%   goodk  - defines the finally chosen k 

warning('off', 'all')
% whichset    = 'all';

mxid        = 5000; 
[val, id]   = max(mean(vertcat(EPH.plhg)));
srt         = id + 200; 
trange      = srt:mxid; 
fastrange   = trange * Fs * 0.1; 
[FMT, fid, fset]  = grp_featmatmaker(trange, whichset, EPH, IMG, fcE); 

% Run non-negative matrix factorisation
%--------------------------------------------------------------------------
% Run multiple times to find optimal number of components
%--------------------------------------------------------------------------
if whichk
R   = [];
for k = 1:1:25
    disp(['Calculating k = ' num2str(k) ' components'])
    [W H i t r] = nmfnnls(FMT, k);
    R(k)        = r;
end
end

% Run again with 'optimum' number of components
%--------------------------------------------------------------------------
rng(45)
[W H i t r] = nmfnnls(FMT, goodk);


% Find the time point with peak loading on factor and sort by that
%--------------------------------------------------------------------------
clear smH h

for h = 1:size(H,1)
    smH(h,:)    = smooth(H(h,:), 300);
    [pks ids]   = findpeaks(double(smH(h,:)));
    [val mid]   = max(pks);
    id(h)       = ids(mid); 
%     plot(smH(h,:), 'color', cols(h,:));     hold on
%     plot([id(h), id(h)], [0, pks(mid)], 'color', cols(h,:)); 
end
[std wsort] = sort(id);