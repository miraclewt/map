function [Z, refvec] = onem(latlim, lonlim, scale)
%ONEM  Construct regular data grid of 1s
%
%   ONEM is not recommended and may be removed in a future release.
%   Instead, use georefcells to construct a geographic raster reference
%   object. Then use ONES to initialize a data grid of the appropriate
%   size:
%
%       R = georefcells(latlim,lonlim,1/scale,1/scale);
%       Z = ones(R.RasterSize);
%
%   [Z, REFVEC] = ONEM(LATLIM, LONLIM, SCALE) constructs a regular
%   data grid consisting entirely of 1s.  The two-element vectors LATLIM
%   and LONLIM define the latitude and longitude limits of the grid, in
%   degrees.  They should be of the form [south north] and [west east],
%   respectively.  The number of rows and columns per degree is set by
%   the scalar value SCALE.  REFVEC is the three-element referencing
%   vector for the data grid.
%
%   See also NANM, SPZEROM, ONES, ZEROM.

% Copyright 1996-2015 The MathWorks, Inc.

[nrows, ncols, refvec] = sizem(latlim, lonlim, scale);
Z = ones(nrows, ncols);
