function grp_plot_modcorr(FMT,fset,fid)
figure
doscatter = 1; 

% Plot correlation between modalities
%--------------------------------------------------------------------------
lo = FMT(fid == find(strcmp(fset, 'lopow')),:);
am = FMT(fid == find(strcmp(fset, 'ampl')),:); 
cols = flip(cbrewer('div', 'RdYlBu', size(lo,2))); 


% Make 2D Histogram
%--------------------------------------------------------------------------
h = hist3([lo(:),am(:)],[100,100]);
sc = 10; 
m = zeros(size(h)*sc); 
for r = 1:size(h,1)
for c = 1:size(h,2)
    m([1:sc]+(r-1)*sc, [1:sc]+(c-1)*sc) = h(r,c);
    
end
end

yl = [min(am(:)), .7];
xl = [min(lo(:)), max(lo(:))]; 

if doscatter   
subplot(2,1,2)
    for c = 1:size(lo,1)
        scatter(lo(c,:), am(c,:), 10, cols, 'filled', 'markerfacealpha', 1), hold on;

    end
    axis square
    ylim(yl); 
    xlim(xl); 
    
    subplot(2,1,1);     
end
lrange = linspace(min(lo(:)),max(lo(:)),size(m,1)); 
arange = linspace(min(am(:)),max(am(:)),size(m,2));
imagesc(lrange,arange,imgaussfilt(m', 30), [0 60]);

colormap(flip(cbrewer('div', 'RdBu', 100))); 
set(gca, 'ydir', 'normal')
axis square
ylim(yl); 
xlim(xl); 




