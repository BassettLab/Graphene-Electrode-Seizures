function grp_plot_nmf(W,H,C,trange,fastrange,fset,fid)

% Sort the loadings by a temporal sequence 
%==========================================================================
% Examine weights of the onset factors
%--------------------------------------------------------------------------
wr = 6;     % range of factor weights to be examined
Featid = [2 3 6 7];

clear tweight
for h = 1:size(H,1)
    tweight(h) = sum(H(h,:) .* [1:size(H,2)]) / sum(H(h,:)); 
end

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

subplot(length(Featid)+2,wr,[1:wr] + wr*(length(Featid)+1))
    plot(C(13).dat(fastrange(1):fastrange(end)))
    xlim([-Inf Inf]);


% Plot factor loading
%--------------------------------------------------------------------------
lincols = cbrewer('qual', 'Dark2', 12); 
clear hndl
for s = 1:length(wsort)
    subplot(length(Featid)+2,wr,[1:wr]); 
    smH = smooth(H(wsort(s),:),10);
	hndl{s} = ['Factor ' num2str(s)];
    plot(smH, 'color', lincols(wsort(s),:), 'linewidth', 1.5); hold on           % plot components
    xlim([-Inf Inf]); 
    [val id] = max(smH);         % find maximum loading for component
    mt(s) = id;                  % time indices of component peak
end

legend(hndl);
xlabel('Time');
ylabel('Factor loading'); 


% Imaging measures
%--------------------------------------------------------------------------
zcut   = 0; 
plotz  = 0; 

for f = 1:length(Featid)
    
    featid = Featid(f); 
    allw    = W(fid == featid,:); allz = zeros(size(allw)); 
    allz(:) = zscore(allw(:)); 
    maxz    = max(allz(:)); 
    mat     = zeros(4);
    for s = 1:wr
        for c = 1:length(C)
            if plotz, mat(C(c).row, C(c).col) = allz(c,wsort(s));
            else,     mat(C(c).row, C(c).col) = allw(c,wsort(s));   end
        end
        subplot(length(Featid)+2,wr,s+(wr*f)),
        if plotz,   collim = [min(allz(:)), max(allz(:))];
        else,       collim = [0 0.1];   
        end
        imagesc(mat, collim);
        title(['Factor ' num2str(s)])
      	ylabel(fset{featid}); 
        colorbar
        axis square
    end
    
end
