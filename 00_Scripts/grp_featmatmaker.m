function [FMT, fid, fset] = grp_featmatmaker(trange, fset, EPH, IMG, fcE)

if ischar(fset) || isstring(fset)
switch fset
    case 'all',     fset = {'plhg', 'lopow', 'hipow', 'lofc', 'hifc', 'eprox', 'ampl'};
    case 'reduced', fset = {'plhg', 'lopow', 'hipow', 'eprox', 'ampl'};
    case 'ephys',   fset = {'plhg', 'lopow', 'hipow', 'lofc', 'hifc'}; 
    case 'calcium', fset = {'eprox', 'ampl'}; 
end
end

FMT = []; 
fid = []; 
for f = 1:length(fset)
if isfield(EPH(1), fset{f}) && isfield(IMG(1), fset{f}), disp(['Ambiguous feature names: ' fset{f}]); 
elseif isfield(EPH(1), fset{f}),        fmt = vertcat(EPH.(fset{f})); 
elseif isfield(IMG(1), fset{f}),        fmt = vertcat(IMG.(fset{f})); 
elseif strcmp(fset{f}, 'lofc'),         fmt = fcE(1).dat;     
elseif strcmp(fset{f}, 'hifc'),         fmt = fcE(2).dat; 
else,                                   disp(['Could not find the right feature: ' fset{f}]); 
end 

% rescale features to 0 - 1, then normalise number of measures (rows)
%--------------------------------------------------------------------------
fmt = (fmt - min(fmt(:))) ./ max(max(fmt - min(fmt(:)))); %./  size(fmt,1);     

fid = [fid; ones(size(fmt,1),1) * f]; 

if ~isempty(FMT), FMT = vertcat(FMT, fmt(:,trange)); 
else,             FMT = fmt(:,trange);
end


end


% 
% featmat = [];
% li      = []; 
% 
% for d = 1:length(D)
%     td = D{d}(:, trange); 
%     td = (td - min(td)) ./ (max(td)-min(td));
%     td(isnan(td)) = 0; 
%     li = [li; ones(size(td,1),1) * d]; 
%     featmat = vertcat(featmat, td); 
% end