function [fd,fh,BdBox,pfix] = Ellipse_Domain


a = 2; b = 1; % ≥§∂Ã∞Î÷·

fd = @(p) p(:,1).^2/a^2+p(:,2).^2/b^2-1;
fh = @huniform; 
BdBox = [-a a -b b];
pfix = [];