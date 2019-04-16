function [el,ei] = edge_link(imb, len, sch)
% Edge linking for binary images.
%
%  [el,ei] = edge_link(imb, len, sch)
%
% This functions links together sequences of 1's in a binary image. It is a 
% modified version of the edge linking algorithm presented in 
% Farag A, Delp E, Edge linking by sequential search. Pattern
% Recognition. 1995; 28: 611--33. The modifications include an user defined
% search region.
%
% INPUT: imb - the binary image under analysis
%               len - the minimum length of an edge
%               sch - a matrix defining the search region
%
% OUTPUT: el - a Lx1 cell array containing vectors that define the L detected 
%                          edges in the image. Each cell contains the x and y co-ordinates 
%                          of each linked edge. 
%                   ei - is an image with all linked edge projected onto a
%                          blank image.
%
% Notes: uses the link_edges.m function
%
% Nathan Stevenson
% NBRG, September 2011

% Initialise binary image
N = size(imb);
M1 = max(abs(sch(2,:)));
M2 = max(abs(sch(1,:)));
imb1 = imb;
clear imb;
imb = zeros(N(1)+2*M2, N(2)+2*M1);
imb(M2+1:N(1)+M2, M1+1:N(2)+M1) = imb1; 
M = size(imb);

% Perform edge-inking procedure
el = cell(1); count = 1;
for ii = 1:M(1);
    for jj = 1:M(2);
        if imb(ii,jj)==1           
            [ifest, imb]= link_edges(imb, jj, ii, sch, M1, M2);
            if length(ifest)>len
                el{count} = ifest';
                count = count+1;
            end
        end
    end
end

clear ii jj

% Generate binary image of linked edges
ei = zeros(N);
for ii = 1:length(el)
    el1 = el{ii};
    for jj = 1:length(el1)
       ei(el1(jj,1), el1(jj,2)) =1; 
    end
    clear el1
end
