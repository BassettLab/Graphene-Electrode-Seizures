function IMG = grp_imgfeat(C, F)
fs = filesep; 

% Extract imaging features underlying individual channel locations
%==========================================================================
if ~exist('edgm'), load([F.analysis fs 'Wave map' fs 'Edge.mat']); end
disp('Calculating proximity to edge'); 
for e = 1:size(edgm,1)
    sme = imgaussfilt(squeeze(double(edgm(e,:,:))), 20);
    for c = 1:length(C)
        IMG(c).eprox(e) = sme(C(c).y, C(c).x);  %% ALERT row/col maps to *y* / *x* I got this wrong many times
	end
end

% Temporary solution to find location within the seizure
%==========================================================================
if ~exist('smtm'), load([F.analysis fs 'Wave map' fs 'Smooth.mat']); end
disp('Calculating locally averaged calcium signal'); 

for b = 1:size(smtm,1)
    for c = 1:length(C)
        IMG(c).ampl(b) = smtm(b, C(c).y, C(c).x);  %% ALERT row/col maps to * y * / * x * I got this wrong many times
	end
end

