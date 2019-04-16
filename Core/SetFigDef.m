function SetFig(w,h,fnt,sz)
% Figure size, font and axes placement
% Usage:
%       SetFigureDefaults(width,height)
%       SetFigureDefaults(width,height,FontName)
%       SetFigureDefaults(width,height,FontName,FontSize)
%
%		Default FontName is Times, default FontSize is 8pt
%       width and height are in cm dimension of axis box (for single subplot)
%
if nargin<3
	fnt='Times';
end
if nargin<4
	sz=8;
end	
sfX=0.75;sfY=0.75;
set(0,'DefaultAxesPosition',[0.15,0.15,sfX,sfY])
set(0,'DefaultAxesFontName',fnt)
set(0,'DefaultAxesFontSize',sz)
set(0,'DefaultTextFontName',fnt)
set(0,'DefaultTextFontSize',sz)

set(gcf,'PaperUnits','centimeters','Units','Centimeters')
p1=get(gcf,'Position');
p2=get(gcf,'PaperPosition');
p1(2)=p1(2)+p1(4)-h/sfY;
p1([3,4])=[w/sfX,h/sfY];
p2([3,4])=[w/sfX,h/sfY];
set(gcf,'Position',p1,'PaperPosition',p2)