function [N, refvec] = egm96geoid(varargin)
% EGM96GEOID  Geoid height from Earth Gravitational Model 1996 (EGM96)
%
%   N = egm96geoid(LAT,LON) returns geoid height in meters from the EGM96
%   geoid model at the locations specified by latitude LAT and longitude
%   LON. LAT and LON must match in size.  Specify LAT and LON in degrees.
%
%   [N,REFVEC] = egm96geoid(SAMPLEFACTOR) returns the grid N of geoid
%   heights from EGM96, sampled every SAMPLEFACTOR elements. To return the
%   entire grid of geoid heights, specify SAMPLEFACTOR as 1. At full
%   resolution, the grid specifies geoid heights at 15 minute intervals.
%   The output argument REFVEC is a referencing vector used to associate
%   each geoid height with a latitude and longitude. The grid of geoid
%   heights has a raster representation of 'postings'.
%   
%   [N,REFVEC] = egm96geoid(SAMPLEFACTOR,LATLIM,LONLIM) returns geoid
%   heights within the specified latitude and longitude limits. Specify
%   longitude limits in the range [-180 180] or [0 360]. For example,
%   LONLIM = [170 190] returns data centered on the 180-degree meridian,
%   while LONLIM = [-10 10] returns geoid heights centered on the Prime
%   Meridian.
%
%   Example
%   -------
%   % Use geoid height at the location of Mount Everest to convert
%   % orthometric height at its summit (relative to sea level, roughly
%   % speaking) to ellipsoidal height (relative to the WGS84 ellipsoid).
%   % All heights are in meters.
%   evlat = 27.988056;
%   evlon = 86.925278;
%   evHeight = 8848;              % Orthometric height of Everest summit
%   N = egm96geoid(evlat,evlon)   % Geoid height at Everest lat-lon
%   h = evHeight + N              % Ellipsoidal height of Everest summit

% Copyright 1996-2019 The MathWorks, Inc.

    narginchk(1,3)
    
    persistent G
    if isempty(G)
        G = map.geodesy.internal.egm96GeoidHeightInterpolant();
    end
    
    if nargin == 2
        % N = egm96geoid(lat,lon)
        lat = varargin{1};
        lon = varargin{2};
        
        validateLatitudeLongitudeAttributes(lat, lon)
        lon = mod(lon,360);
        
        if~isequal(size(lat),size(lon))
            error(message('map:validate:inconsistentSizes','LAT','LON'))
        end
        
        N = double(G(lat, lon));
        
        if nargout > 1
            warning(message('map:validate:unnecessaryOutput','REFVEC'))
            refvec = [];
        end
    else
        % [N, refvec] = egm96geoid(samplefactor)
        % [N, refvec] = egm96geoid(samplefactor, latlim, lonlim)
        [samplefactor, latlim, lonlim] = parseInputs(varargin);
        
        % Import geoid height grid N and referencing object R.
        [N, R] = importEGM96(G, latlim, lonlim, samplefactor);
        
        refvec = [];
        if ~isempty(N)
            % Convert referencing object to referencing vector.
            refvec = geopostings2refvec(R);
        end
    end
end


function validateLatitudeLongitudeAttributes(lat, lon)
% Validate lat, lon type and limits (see validation in geopeaks)
    try
        if isfloat(lat)
            validateattributes(lat(~isnan(lat)), ...
                {'double','single'},{'real','>=',-90,'<=',90,'nonsparse'},'','LAT')
        else
            validateattributes(lat,{'double','single'},{},'','LAT')
        end
        validateattributes(lon,{'double','single'},{'real','nonsparse'},'','LON')
        if any(isinf(lon))
            % Allow lon to include NaN, but not +/- Inf. (Can't use the
            % 'finite' attribute because it will reject NaN.)
            validateattributes(lon(~isnan(lon)), ...
                {'double','single'},{'real','finite','nonsparse'},'','LON')
        end
    catch e
        throwAsCaller(e)
    end
end


function [samplefactor, latlim, lonlim] = parseInputs(inputs)

numInputs = numel(inputs);
try
    % Validate SAMPLEFACTOR.
    samplefactor = inputs{1};
    maxSampleFactor = 4*180;
    validateattributes(samplefactor, {'numeric'}, ...
        {'positive','scalar','integer','<=', maxSampleFactor}, ...
        mfilename, 'SAMPLEFACTOR', 1)
    
    % Set or validate LATLIM.
    if numInputs < 2
        latlim = [-90 90];
    else
        latlim = inputs{2};
        
        if isempty(latlim)
            latlim = [-90 90];
        end
        
        validateattributes(latlim, {'double'}, ...
            {'real','vector','finite','>=',-90,'<=',90}, mfilename, 'LATLIM')
        
        map.internal.assert(numel(latlim) == 2, ...
            'map:validate:expectedTwoElementVector', 'LATLIM');
        
        latlim = sort(latlim);
        
        map.internal.assert(latlim(1) < latlim(2), ...
            'map:maplimits:expectedAscendingLatlim')
        
        latlim = latlim(:)';
    end
    
    % Set or validate LONLIM.
    if numInputs < 3
        lonlim = [0 360];
    else
        lonlim = inputs{3};
        
        if isempty(lonlim)
            lonlim = [0 360];
        end
        
        validateattributes(lonlim, {'double'}, ...
            {'real','vector','finite'}, mfilename, 'LONLIM')
        
        map.internal.assert(numel(lonlim) == 2, ...
            'map:validate:expectedTwoElementVector', 'LONLIM');
        
        map.internal.assert(all(lonlim <= 360), ...
            'map:validate:expectedRange', 'LONLIM',  '0', 'lonlim', '360')
        
        if any(lonlim < 0) && (any(lonlim < -180) || any(lonlim > 180))
            error(message('map:validate:expectedRange', ...
                'LONLIM', '-180', 'lonlim', '180'))
        end
        
        if lonlim(2) < lonlim(1)
            lonlim(2) = lonlim(2) + 360;
        end
        
        lonlim = lonlim(:)';
    end
catch me
    throwAsCaller(me)
end
end


function [N, R] = importEGM96(G, latlim, lonlim, samplefactor)
% Import, crop, and subsample the full geoid height grid, returning the
% output grid and a geographic postings raster reference object.

    % Extract the grid from the padded version in the gridded interpolant
    % for EGM96 height, G.
    glatlim = [-90 90];
    glonlim = [0 360];
    [nrows, ncols] = size(G.Values);
    p = 8; % Extent of padding on all 4 sides of G.Values
    globalN = G.Values(1+p : nrows-p, 1+p : ncols-p);
    globalR = georefpostings(glatlim, glonlim, size(globalN));
    [R, rows, cols] = setupCropAndSubsampleForGlobalPostings( ...
        globalR, latlim, lonlim, samplefactor);
    N = double(globalN(rows, cols));
end


function [R, rows, cols] = setupCropAndSubsampleForGlobalPostings( ...
    globalR, latlim, lonlim, samplefactor)

% This function can be very simple for two reasons: (1) the input globalR
% is assumed to correspond to a raster of global extent, and (2) both the
% input and output rasters are snapped to multiples of the input sample
% spacings.  These conditions guarantee that every output sample exactly
% coincides (in latitude and longitude) with an input sample.  The
% function also assumes that globalR.ColumnsStartFrom is 'south' and
% globalR.RowsStartFrom is 'west.'
%
% See also map.rasterref.GeographicRasterReference/setupCropAndSubsample

    % Pre-snap the limits because, in the special cases in which
    % diff(latlim) or diff(lonlim) is an exact multiple of latspace or
    % lonspace, georefpostings won't snap the limits for us.
    [latlim, lonlim] = snapLimitsToRaster(globalR, latlim, lonlim);

    latspace = samplefactor * globalR.SampleSpacingInLatitude;
    lonspace = samplefactor * globalR.SampleSpacingInLongitude;
    
    R = georefpostings(latlim, lonlim, latspace, lonspace);

    rows = 1:R.RasterSize(1);
    cols = 1:R.RasterSize(2);

    if ~isequal(R, globalR)
        lat = intrinsicYToLatitude(R, rows);
        lon = intrinsicXToLongitude(R, cols);

        rows = round(latitudeToIntrinsicY(globalR, lat));
        cols = round(longitudeToIntrinsicX(globalR, lon));
    end
end


function [latlim, lonlim] = snapLimitsToRaster(R, latlim, lonlim)
% Assume that R.ColumnsStartFrom is 'south', R.RowsStartFrom is 'west',
% R.LongitudeLimits(1) = 0, R.SampleSpacingInLongitude = 1/4:
    
    ylimi = latitudeToIntrinsicY(R, latlim);
    ylimi = [floor(ylimi(1)) ceil(ylimi(2))];
    latlim = intrinsicYToLatitude(R, ylimi);
    
    % Need to scale and shift without wrapping:
    % avoid the longitudeToIntrinsicX method.
    xlimi = 1 + 4 * lonlim;
    xlimi = [floor(xlimi(1)) ceil(xlimi(2))];
    lonlim = intrinsicXToLongitude(R, xlimi);
end


function refvec = geopostings2refvec(R)
%   This function assumes:
%
%       'postings'
%       columns start from west
%       rows start from south
%       sample density is the same in both dimensions
%
%   This is not a general purpose conversion. Instead, the limits are
%   extrapolated 1/2 sample north and west of their true locations, for
%   consistency with the previous behavior of egm96geoid.

    cellsPerDegree = sampleDensity(R);

    % Step north 1/2 cell/sample-spacing relative to northern limit.
    northernLimit = R.LatitudeLimits(2) + 0.5 / cellsPerDegree;

    % Step west 1/2 cell/sample-spacing relative to the western limit.
    westernLimit  = R.LongitudeLimits(1) - 0.5 / cellsPerDegree;

    refvec = [cellsPerDegree northernLimit westernLimit];
end
