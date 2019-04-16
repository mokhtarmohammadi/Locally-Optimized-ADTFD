function [ml, imbm] = link_edges(imb, rx, ry, sch, M1, M2)
% Link edges 
%
% [ml, imbm] = link_edges(imb, rx, ry, sch, M1, M2)
%
%This function tracks a edge through the image via a particular
% neighbourhood search pattern
%
% INPUTS: imb - the binary image to be linked
%                  [rx, ry] - the starting point of the component search
%                  sch - the user-defined search region
%                  [M1, M2] - the size of the binary image
%
% OUTPUTS: ml - the [x,y] co-ordinates of the linked component
%                     imbm - the updated binary image with the linked
%                                   component removed.
%
% UQCCR, level 04
% Nathan Stevenson
% August 2008

c1 = 1; 
mlx = rx;
mly = ry;
flag = 1;
imb(ry,rx)=0;
while flag==1
    for kk = 1:length(sch)
           x = mlx(c1)+sch(2,kk);
           y = mly(c1)+sch(1,kk);
           ref(kk) = imb(y, x);
    end
    qq = find(ref==1, 1);
    if isempty(qq)~=1        
        rx = mlx(c1)+sch(2,qq); 
        ry = mly(c1)+sch(1,qq);
        imb(ry,rx)=0;
        c1 = c1+1;
        mlx(c1) = rx; mly(c1) = ry;
        flag = 1;
    else
        flag = 0;
    end
    clear qq ref
end
ml = [mly-M2 ; mlx-M1];
imbm = imb;

