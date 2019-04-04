% Calculate z-score from baseline


f   = 3;
mua = zeros(size(fdat{f})); 

for c = 1:size(fdat{f},1)
  
dat  = fdat{f}(c,:);  
[zbas mu sig] = zscore(fbdat{f}(c,:)); 
mua(c,:)    = (dat - mu) > 1.5*sig;

end

imagesc(mua)
