%% Plot W matrix features
%--------------------------------------------------------------------------
colormap(flip(cbrewer('div', 'RdBu', '100'))); 
figure(1)
for f = 1:7
    subplot(7,1,f)
    imagesc(W(find(fid == f),wsort), [0 0.1])
end

figure(2)
colormap(flip(cbrewer('div', 'RdBu', '100'))); 
imagesc(H(wsort,:))


%% Plot correlation between modalities
%--------------------------------------------------------------------------
lo = FMT(fid == find(strcmp(fset, 'lopow')),:);
am = FMT(fid == find(strcmp(fset, 'ampl')),:); 
cols = flip(cbrewer('div', 'RdYlBu', size(lo,2))); 

for c = 1:size(lo,1)
    scatter(lo(c,:), am(c,:), 200, cols, 'filled', 'markerfacealpha', .02), hold on;
    
end

