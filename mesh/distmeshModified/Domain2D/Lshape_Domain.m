function [fd,fh,BdBox,pfix] = Lshape_Domain

a = -1; b = 1; c = -1; d = 1;
xm = (a+b)/2;  ym = (c+d)/2;

fd = @fd1;
fh = @huniform;
BdBox = [a b c d];

% % cut off the lower right
% pfix = [a c; xm c; xm ym; b ym; b d; a d];
%     function dist = fd1(p)
%         % The following code is not correct
%         % ddiff can not remove the distance==0
%         % d1 = drectangle(p,a,b,c,d);
%         % d2 = drectangle(p,xm,b,c,ym);
%         % dist = ddiff(d1,d2); 
% 
%         d1 = drectangle(p,a,b,ym,d);
%         d2 = drectangle(p,a,xm,c,ym);
%         dist = dunion(d1,d2);
%         
% 
%     end

% cut off the upper right
pfix = [a c; b c; b ym; xm ym; xm d; a d];
    function dist = fd1(p)
        d1 = drectangle(p,a,xm,ym,d);
        d2 = drectangle(p,a,b,c,ym);
        dist = dunion(d1,d2);
    end
end