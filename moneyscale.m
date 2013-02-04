function g = moneyscale(m)
%WINTER Shades of blue and green color map
%   WINTER(M) returns an M-by-3 matrix containing a "winter" colormap.
%   WINTER, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%       colormap(winter)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.8.4.2 $  $Date: 2005/06/21 19:31:52 $

if nargin < 1, m = size(get(gcf,'colormap'),1); end
g = (0:m-1)'/max(m-1,1);
g = [zeros(m,1) g zeros(m,1)];
