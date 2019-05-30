fs  = filesep; 
F   = grp_housekeeping; 
load([F.data fs 'Electrophysiology' fs 'Channel_Data.mat'])

%% Get EEG Datafeatures
%==========================================================================
pullieeg = 0; 
physfeat = 1; 

if pullieeg,    C   = grp_electrophys(F);   end
if ~exist('C'), load([F.data fs 'Electrophysiology' fs 'Channel_Data.mat']); end
if physfeat,    [EPH, fcE] = grp_ephysfeat(C);     end

% Get Imaging Datafeatures
%==========================================================================
estwave = 0;
getfeat = 1;
 
if estwave,     F   = grp_imaging(F);       end
if getfeat,     IMG = grp_imgfeat(C, F);    end

%% Matrix decomposition
%==========================================================================
whichk = 0; 

mxid        = 8000; 
[val, id]   = max(mean(vertcat(EPH.plhg)));
srt         = id + 200; 
id          = srt:mxid; 

[FMT, fid]   = grp_featmatmaker(...
                {vertcat(EPH.plhg), ... 
                 vertcat(EPH.spamp),...
                 vertcat(EPH.spcnt),...
                 vertcat(EPH.lopow),...
                 vertcat(EPH.hipow),...
                 fcE(1).dat, fcE(2).dat, ...
                 vertcat(IMG.eprox),...
                 vertcat(IMG.ampl)}, id); 

featnames = {'PLHG', 'SPAmp', 'SPFrq', 'LoPow', 'HiPow', 'FcLo', 'FcHi', 'WvFrnt', 'CaAmp'}; 
feattypes = [1 1 1 1 1 2 2 3 3];    % 1 = ephys univar; 2 = ephys bivar; 3 = Calcium

% Run non-negative matrix factorisation
%--------------------------------------------------------------------------
% Run multiple times to find optimal number of components
%--------------------------------------------------------------------------
if whichk
R   = [];
for k = 1:50
    disp(['Calculating k = ' num2str(k) ' components'])
    [W H i t r] = nmfnnls(FMT, k);
    R(k)        = r;
end
end

% Run again with 'optimum' number of components
%--------------------------------------------------------------------------
rng(45)
[W H i t r] = nmfnnls(FMT, 7);

%% Sort the loadings by a temporal sequence 
%==========================================================================
% Examine weights of the onset factors
%--------------------------------------------------------------------------
whichone = 'corrweights'; 
whichone = 'univars'; 

Featid = [4,5];
if any(Featid == 6 | Featid == 7)
    whichone = 'corrweights';
else,   whichone = 'univars'; 
end


clear tweight
for h = 1:size(H,1)
    tweight(h) = sum(H(h,:) .* [1:size(H,2)]) / sum(H(h,:)); 
end
[std wsort] = sort(tweight);

% Identify time points of maximum component loading
%--------------------------------------------------------------------------
lincols = cbrewer('qual', 'Dark2', 12); 
clear hndl
for s = 1:length(wsort)
    subplot(length(Featid)+1,3,[1:3]); 
    smH = smooth(H(wsort(s),:),300);
	hndl{s} = ['Factor ' num2str(wsort(s))];
    plot(smH, 'color', lincols(wsort(s),:), 'linewidth', 1.5); hold on           % plot components
    [val id] = max(smH);         % find maximum loading for component
    mt(s) = id;                  % time indices of component peak
end

legend(hndl);
xlabel('Time');
ylabel('Factor loading'); 

% Correlation weights 
%--------------------------------------------------------------------------
if strcmp(whichone, 'corrweights')
dum     = ones(length(C)); dum = zeros(length(C)) + triu(dum,1); dum = find(dum); 
ddum    = zeros(16);

Featid  = [6 7];    % 6: low freq corr; 7: high freq corr
zcut    = 1;        % zscore cut off when looking across all weight loadings
fnames  = {'LFP correlation', 'HFO correlation'}; 

for f = 1:length(Featid)
    featid = Featid(f); 
    allw 	= W(fid == featid,:);
    allw(zscore(allw(:)) < zcut) = 0; 
    for s = 1:3
        ddum(dum)   = allw(:, wsort(s)); 
        A           = ddum + ddum';

        subplot(3,3,s+(3*f))       
            G           = graph(A);  
            jx = randn(1,16)*0;
            jy = randn(1,16)*0;
            plot(G,'XData',[C.x]+jx,'YData',[C.y]+jy, 'LineWidth', G.Edges.Weight*30, ...
                 'NodeLabel', [C.id], 'EdgeColor', lincols(wsort(s),:))
            set(gca, 'YDir', 'Reverse');
            axis square
            title(['Factor ' num2str(wsort(s))])
            ylabel(fnames{f}); 
    end
end
end

% Imaging measures
%--------------------------------------------------------------------------
if strcmp(whichone, 'univars')
zcut   = 0; 

for f = 1:length(Featid)
    featid = Featid(f); 
    allw    = W(fid == featid,:); allz = zeros(size(allw)); 
    allz(:) = zscore(allw(:)); 
    maxz    = max(allz(:)); 
    mat     = zeros(4);
    for s = 1:3
        for c = 1:length(C)
            mat(C(c).row, C(c).col) = allz(c,wsort(s)); 
        end
        subplot(length(Featid)+1,3,s+(3*f)),
        imagesc(mat, [min(allz(:)), max(allz(:))]);
        title(['Factor ' num2str(wsort(s))])
      	ylabel(featnames{featid}); 
        colorbar
        axis square
    end
    
end

end

%% Plot trajectories
%==========================================================================
f   = 5;
cols = log(mean(FMT(fid == f,:)));
cols = 1:size(H,2); 
scatter3(H(wsort(1),:), H(wsort(2),:), H(wsort(3),:), [], cols, 'filled', ...
         'markerfacealpha', 0.8); 
title(['Trajectory color coded by ' featnames{f}]);
xlabel('Baseline component');
ylabel('Seizure onset component'); 
zlabel('Early Ictal component'); 


%% GRAVEYARD

alli = unique(id); 
for i = 1:length(alli)
    tid = find(id == alli(i)); 
    Wg{i} = W(tid,:);
    if (alli(i) == 8) || (alli(i) == 9),  cstr = 'r',     
    else,                       cstr = 'k';     end
    
    plot(mean(Wg{i}), cstr), hold on
    mW(i,:) = mean(Wg{i}); 
end
%%
rng(45)
Z = linkage(W, 'complete');
T = cluster(Z,'maxclust',10);
% dendrogram(Z)
[std tsort] = sort(T);
imagesc(corr(W(tsort,:)'));

% Percentage of channel features in each 
caid = find(feattypes == 3);
for c = 1:length(caid)
    tcaid = find(fid == caid(c));
    allts = unique(T);
    for t = 1:length(allts)
        rat(c,t) = length(find(fid(T == t) == caid(c))) / length(tcaid);
    end
    
end

imagesc(W(tsort,:))
set(gca, 'YTick', find([1; diff(sort(T))]))
