function grp_plot_snapshot(F, C, trange, fastrange, Fs, fset, W, H, fid, wsort)
plotgrid = 0; 
fs       = filesep; 
srt      = trange(1); 

ts = [500, 950, 1250, 2200, 3500, 4500];            % These are derived from the NMF analysis
if plotgrid,    ts = fix(linspace(1,5000,36));   end

for ti = 1:length(ts)

t = ts(ti); 


% Load detected edges
%--------------------------------------------------------------------------
if ~exist('edgm'), load([F.analysis fs 'Wave map' fs 'Edge.mat']); end
% if ~plotgrid
% emap = edgm(t,:,:);
% imagesc(squeeze(edgm(t,:,:)))
% end

% Plot map of normalised fluorescence and seizure wavefront
%--------------------------------------------------------------------------
SZ      = [F.data fs 'Graphene Imaging' fs 'bath1']; 
tiflist = dir([SZ fs '*.tif']); 
tiflist = {tiflist.name}; 
if ~exist('semmap'),    load([F.data fs 'SEM_map.mat']);    end
if ~exist('avgmap'),    load([F.data fs 'AVG_map.mat']);    end

loadid  = srt + t; 
tif  	= imread([SZ fs tiflist{loadid}]); 
tif     = (single(tif) - avgmap)./semmap; 

% Plot map of calcium amplitude
%--------------------------------------------------------------------------
if plotgrid,    subplot(ceil(length(ts)/6), 6, ti);
else,           figure(ti), subplot(3,1,[1,2]);        end

disp(['I''m plotting ' num2str(t) ' at ' num2str(ti)])
colormap(flip(cbrewer('div', 'RdGy', 100))),
imagesc(tif, [0 1200]);     hold on

% Format 
%--------------------------------------------------------------------------
axis square
title(['Plot timeframe no ' num2str(t)]); 
set(gca, 'ydir', 'reverse')
xlim([0 512]);
ylim([0 512]); 
colorbar

seq = flip(cbrewer('div', 'PuOr', 6)); 
seq = seq([1 2 5 6],:); 


if ~plotgrid
    % Find and plot
    edgi = find(squeeze(edgm(loadid,:,:)));  
    [r, c] = ind2sub([size(edgm,2), size(edgm,3)], edgi);

    % Plot settings 
    cl = cbrewer('qual', 'Set1', 10); hold on

    % Plot
    scatter(c,r, 20, cl(2,:), 'filled');  
    hold on

    % Show segments of electrophysiology
    %--------------------------------------------------------------------------
    k = 0;
    seq = flip(cbrewer('div', 'PuOr', 6)); 
    seq = seq([1 2 5 6],:); 
    for c = [5 6 7 8]               % WHELP introducing subplot(3,1,[1,2]) here deletes the existing plots. What Matlab Gods have I angered today?
        k = k + 1;
        scatter(C(c).x, C(c).y, 30, seq(k,:), 'filled'); 
        hold on
        scatter(C(c).x, C(c).y, 30, 'k'); 
        xlim([-Inf Inf])
    end
    
    k = 0;
    for c = [5 6 7 8]
        k = k + 1;
        subplot(3,1,3)      % Plot example traces
            plot(C(c).filt{1}(fastrange(t):fastrange(t)+Fs * 3) + k*300, 'color', seq(k,:));
            hold on
            ylim([-1500,1500]); 
            xlim([-Inf Inf]); 
    end

    [xc,lg]     = xcorr(C(5).filt{1}(fastrange(t):fastrange(t)+Fs * 3), C(8).filt{1}(fastrange(t):fastrange(t)+Fs * 3));
    [val, id]   = max(xc);
    subplot(3,1,3), title(num2str(['Lag between left and right edge: ' num2str(lg(id))])); 

    set(gcf, 'Position', [200, 200, 300, 350])
    end 

end


% Correlation weights 
%--------------------------------------------------------------------------
featid  = find(strcmp(fset, 'hifc'));
dum     = ones(length(C)); dum = zeros(length(C)) + triu(dum,1); dum = find(dum); 
ddum    = zeros(16);
zcut    = -.5;        % zscore cut off when looking across all weight loadings
wr      = size(W,2); 

allw 	= W(find(fid == featid),:);
allz    = reshape(zscore(allw(:)), size(allw)) > zcut;
for s = 1:wr
    ddum(dum)   = allw(:, wsort(s)) .* allz(:,wsort(s)); 
    A           = ddum + ddum';

    subplot(2,3,s)
        G               = graph(A); 
        p               = plot(G, 'layout', 'force');
        p.Marker        = 'o';
        p.MarkerSize    = 10;
        gry             = [.5 .5 .5];
        p.NodeColor     = [gry; gry; gry; gry; ...
                           seq(1,:); seq(2,:); seq(3,:); seq(4,:); 
                           gry; gry; gry; gry; ...
                           gry; gry; gry; gry];
        p.NodeLabel     = [];            
    axis square
%     ylim([-4 4]); 
%     xlim([-4 4]); 
end

