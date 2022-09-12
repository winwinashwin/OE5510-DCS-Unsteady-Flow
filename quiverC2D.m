% Copyright (c) 2021, Jonas Krimmer
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% 
% * Neither the name of  nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function quiverC2D(x,y,u,v,varargin)
%quiverC2D creates a 2D quiver plot and adds a color coding. The color coding is
%given by the absolut values of the component vectors. Large values result in colors 
%from the upper end of the used colormap. Plotting parameters have to be changed within 
%the function in this version. Values equal to NaN or inf are set to zero.
%In case of complex arrays the absolute value is used.
% 
%   INPUT:
%       x - array, x components of initial points
%       y - array, y components of initial points
%       u - array, x components of arrows
%       v - array, y components of arrows
%       scale - number, if > 0, automatic scaling is used to avoid overlap
%           of adjacent arrows. If scale = 0, the scaling is disabled.
%       maxNumArrows - a positive integer (non-integer should work as well)
%           number limiting the maximum number of plotted arrows. Since vectors
%           of length zero aren't plotted and the matrices are sampled
%           uniformly, it's possible that fewer arrows will be visible in the
%           end. If maxNumArrows is not given its value is set to inf.
% 
%   OUTPUT:
%       none
% 
%   WARNING!: This function might not create arrows with the same length as
%   would be obtained using quiver(x, y, u, v). I do not know whether the
%   scaling method is the same is in quiver.
% 
% --------------------------------------------------------------------------------------
% 
%   EXAMPLE: 
%       load wind
%       quiverC2D(x(:,:,1),y(:,:,1),u(:,:,1),v(:,:,1),1);
%       set(gca, 'Color', 'black');
%       axis tight
%       axis equal
%   
% --------------------------------------------------------------------------------------
% 
%   See also: QUIVER, LINESPEC, COLORMAP.
% 
% 

%% parse arguments

p = inputParser;
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addRequired(p, 'u', @isnumeric);
addRequired(p, 'v', @isnumeric);
addOptional(p, 'scale', 1, @isscalar);
addOptional(p, 'MaxNumArrows', inf, @isscalar);
addParameter(p, 'Normalize', false);
addParameter(p, 'LineWidth', 0.5, @isscalar);
addParameter(p, 'MaxHeadSize', 1, @isscalar);
addParameter(p, 'Colorbar', false, @islogical);
addParameter(p, 'Parent', gca(), @isgraphics);

parse(p, x, y, u, v, varargin{:});

scale = p.Results.scale;
maxnumarrows = p.Results.MaxNumArrows;
normalize = p.Results.Normalize;
lw = p.Results.LineWidth;
hs = p.Results.MaxHeadSize;
ax = p.Results.Parent;

assert(scale >= 0);

if ~isequal(size(x),size(y),size(u),size(v))
    error('X,Y,Z,U,V,W have to be arrays of the same size.');
end

%% initialization
if numel(u) > maxnumarrows
    if isvector(u)
        N = ceil(numel(u)/maxnumarrows);
    else
        N = ceil(nthroot(numel(u)/maxnumarrows, ndims(u)));
    end
    
    if isvector(x)
        idx = {1:N:length(x)};
    else
        idx = arrayfun(@(x) 1:N:x, size(x), 'UniformOutput', false);
    end
    
    x = x(idx{:});
    y = y(idx{:});
    u = u(idx{:});
    v = v(idx{:});
end


%% colormap definition
C = colormap;
ncolors = size(C, 1);
I = sqrt(u.^2 + v.^2);

if verLessThan('matlab','9.6')
    Imax = max(I(:), [], 'omitnan');
    Imin = min(I(:), [], 'omitnan');
else
    Imax = max(I, [], 'all', 'omitnan');
    Imin = min(I, [], 'all', 'omitnan');
end
% assume that the colormap matrix contains the colors in its rows
Ic = round(I/Imax*ncolors);
Ic(Ic == 0) = 1;

% scale not to overlap
if scale > 0
    [~,dxi] = uniquetol(x);
    dx = diff(x(sort(dxi)));
    [~,dyi] = uniquetol(y);
    dy = diff(y(sort(dyi)));
    
    % avg difference
    dr_mean = mean(abs([dx; dy]),'all');
    
    lenscale = scale*dr_mean/Imax;
else
    lenscale = 1;
end

if normalize
    lenscale = lenscale./I*Imax;
end

u = u.*lenscale;
v = v.*lenscale;


%% plotting
ishold_flag = ishold(ax);
hold(ax, 'on');

if numel(u) > ncolors
    % let's take an intelligent approach: group all values which belong to
    % the same color value
    for k = 1:ncolors
        mask = (Ic == k);    
        plt = quiver(ax, x(mask), y(mask), u(mask), v(mask), 0, ...
            'Color', C(k, :), 'LineWidth', lw, 'MaxHeadSize', hs);
    end
else
    % if there are so few values, it will be more efficient to plot each
    % arrow individually
    % linear indexing!
    for k = 1:numel(u)
        plt = quiver(ax, x(k), y(k), u(k), v(k), 0, 'Color', C(Ic(k),:), ...
            'LineWidth', lw, 'MaxHeadSize', hs);
    end
end

if ~ishold_flag
    hold(ax, 'off');
end

if p.Results.Colorbar
    cb = colorbar(ax);
    tlabels = cb.TickLabels;
    tlabels = sprintfc('%1.3g ',linspace(Imin, Imax, length(tlabels)));
    cb.TickLabels = tlabels;
end

end