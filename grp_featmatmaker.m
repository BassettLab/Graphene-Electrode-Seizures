function [featmat, li] = grp_featmatmaker(D,id)

featmat = [];
li      = []; 

for d = 1:length(D)
    td = D{d}(:, id); 
    td = (td - min(td)) ./ (max(td)-min(td));
    td(isnan(td)) = 0; 
    li = [li; ones(size(td,1),1) * d]; 
    featmat = vertcat(featmat, td); 
end