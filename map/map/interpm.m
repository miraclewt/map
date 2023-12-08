function [lat,lon] = interpm(lat, lon, maxsep, method, angleunit)
%INTERPM  Densify latitude-longitude sampling in lines or polygons
%
%   [DENSELAT,DENSELON] = INTERPM(LAT,LON,MAXSEP) inserts additional points
%   into vectors latitude (LAT) and longitude (LON) if adjacent points are
%   separated by more than the tolerance MAXSEP. New points are spaced
%   linearly, in both latitude and longitude, between existing points.
%   LAT, LON, and MAXSEP are specified in degrees.
%
%   [DENSELAT,DENSELON] = INTERPM(__,METHOD), where METHOD is 'gc', inserts
%   new points along great circles. Where METHOD is 'rh', INTERPM inserts
%   new points along rhumb lines. The default method, linear spacing in
%   latitude and longitude, can be specified as 'lin'.
%
%   [DENSELAT,DENSELON] = INTERPM(__,METHOD,ANGLEUNIT), where ANGLEUNIT is
%   'radians', specifies LAT, LON, and MAXSEP in radians.
%
%  See also GEOINTERP, INTRPLAT, INTRPLON, LINSPACE

% Copyright 1996-2019 The MathWorks, Inc.

    narginchk(3,5)
    
    if ~isequal(size(lat),size(lon))
        error(message('map:validate:inconsistentSizes2','INTERPM','LAT','LON'))
    end
    
    validateattributes(maxsep, {'double'}, {'scalar'}, 'INTERPM', 'MAXSEP', 3)
    
    lat = ignoreComplex(lat, 'interpm', 'lat');
    lon = ignoreComplex(lon, 'interpm', 'lon');
    maxsep = ignoreComplex(maxsep, 'interpm', 'maxdiff');
    
    if nargin < 4
        method = 'lin';
    else
        method = validatestring(method, {'gc','rh','lin'}, 'INTERPM', 'METHOD', 4);
    end
    
    if nargin < 5
        angleunit = 'degrees';
    else
        angleunit = checkangleunits(angleunit);
    end
    
    switch method
        case 'gc'
            [lat, lon] = densifyGreatCircle(lat, lon, maxsep, angleunit);
        case 'rh'
            [lat, lon] = densifyRhumbline(lat, lon, maxsep, angleunit);
        otherwise
            [lat, lon] = densifyLinear(lat, lon, maxsep);
    end
end


function [lat, lon] = densifyGreatCircle(lat, lon, maxsep, angleUnit)
% Insert additional vertices where needed, along great circle arcs
% connecting adjacent pairs of vertices.

    % Ensure column vectors.
    lat = lat(:);
    lon = lon(:);
    
    % Compute max angular separation between each pair of adjacent vertices.
    separation = max([abs(diff(lat))'; abs(diff(lon))'])';
    
    % Find separations that exceed the specified limit.
    indx = find(separation > maxsep);
    if ~isempty(indx)
        steps = ceil(separation(indx)/maxsep);
        [lat, lon] = toDegrees(angleUnit, lat, lon);
        [lat, lon] = insertDenseVertices(lat, lon, indx, steps, @linspaceGreatCircle);
        lon = unwrapMultipart(lon,'degree');
        [lat, lon] = fromDegrees(angleUnit, lat, lon);
    end
end


function [latinsert, loninsert] = linspaceGreatCircle(lat, lon, e, step)
% Work in degrees using the great functions from map.geodesy.internal
    lat1 = lat(e);
    lon1 = lon(e);
    lat2 = lat(e+1);
    lon2 = lon(e+1);
    [S,az] = map.geodesy.internal.greatCircleDistance(lat1,lon1,lat2,lon2);
    S = (1:step-1) * S / step;
    [latinsert, loninsert] = map.geodesy.internal.greatCircleTrace(lat1,lon1,S,az);
end


function [lat, lon] = densifyRhumbline(lat, lon, maxsep, angleUnit)
% Insert additional vertices where needed, along rumbline arcs
% connecting adjacent pairs of vertices.

    % Ensure column vectors.
    lat = lat(:);
    lon = lon(:);
    
    % Compute max angular separation between each pair of adjacent vertices.
    separation = max([abs(diff(lat))'; abs(diff(lon))'])';
    
    % Find separations that exceed the specified limit.
    indx = find(separation > maxsep);
    
    if ~isempty(indx)
        % Interpolate additional vertices where needed.
        steps = ceil(separation(indx)/maxsep);
        [lat, lon] = toRadians(angleUnit, lat, lon);
        [lat, lon] = insertDenseVertices(lat, lon, indx, steps, @linspaceRhumbline);
        lon = unwrapMultipart(lon);
        [lat, lon] = fromRadians(angleUnit, lat, lon);
    end
end


function [latinsert, loninsert] = linspaceRhumbline(lat, lon, e, step)
    phi1 = lat(e);
    phi2 = lat(e+1);
    lambda1 = lon(e);
    lambda2 = lon(e+1);
    
    %  Compute distance and azimuth, then sample at intermediate distances.
    [S, az] = rhumblineinv(phi1, lambda1, phi2, lambda2, 1);
    S = (1:step-1)' * S / step;
    
    %  Pass column vectors to rhumblinefwd.
    c = ones(size(S));
    [latinsert, loninsert] = rhumblinefwd(phi1(c,1), lambda1(c,1), az(c,1), S, 1);
    loninsert = wrapToPi(loninsert);
end


function [lat, lon] = densifyLinear(lat, lon, maxsep)
% Insert additional vertices where needed, between adjacent pairs of
% vertices, with independent linear spacing in latitude and longitude.

    % Ensure column vectors.
    lat = lat(:);
    lon = lon(:);
    
    % Compute max angular separation between each pair of adjacent vertices.
    separation = max([abs(diff(lat))'; abs(diff(lon))'])';
    
    % Find separations that exceed the specified limit.
    indx = find(separation > maxsep);
    
    if ~isempty(indx)
        % Interpolate additional vertices where needed.
        steps = ceil(separation(indx)/maxsep);
        [lat, lon] = insertDenseVertices(lat, lon, indx, steps, @linspaceLatitudeLongitude);
    end
end


function [latinsert, loninsert] = linspaceLatitudeLongitude(lat, lon, e, step)
% Independent linear spacing in laitude and longitude
    factors = (1:step-1)';
    latinsert = ((lat(e+1)-lat(e))/step)*factors + lat(e);
    loninsert = ((lon(e+1)-lon(e))/step)*factors + lon(e);
end


function [ilat, ilon] = insertDenseVertices(lat, lon, indx, steps, linspacefcn)
% Insert new, intermediate vertices when densifying a curve -- supports
% densifyGreatCircle, densifyRhumbline, and densifyLinear, each of which
% provides its own "linspacefcn" function handle.

    % Length of output vectors
    N = numel(lat) + sum(steps) - numel(steps);
    
    % Pre-allocate outputs.
    ilat = zeros(N,1);
    ilon = zeros(N,1);
    
    s = 1;
    sNew = 1;
    for k = 1:numel(indx)
        % Insert vertices between lat(indx(k)), lon(indx(k))
        % and lat(1 + indx(k)), lon(1 + indx(k)).
        
        % Copy the elements that precede the current insertion point.
        e = indx(k);
        eNew = sNew + e - s;
        ilat(sNew:eNew) = lat(s:e);
        ilon(sNew:eNew) = lon(s:e);
        
        % Advance s to be ready for the next iteration (or final copy),
        % and also for use in the insertion step below.
        s = 1 + e;
        
        % Compute the new vertices that are to be inserted
        [latinsert, loninsert] = linspacefcn(lat, lon, e, steps(k));
        
        % Copy the new vertices into the output arrays.
        sNew = eNew + 1;
        eNew = eNew + numel(latinsert);
        ilat(sNew:eNew) = latinsert;
        ilon(sNew:eNew) = loninsert;
        
        % Advance sNew to be ready for the next iteration (or final copy).
        sNew = eNew + 1;
    end
    
    % Append any vertices that come after the last insertion point. Note
    % that if indx is empty, we'll have skipped the loop, s will equal 1, and
    % we'll be making an exact copy of the inputs.
    ilat(sNew:end) = lat(s:end);
    ilon(sNew:end) = lon(s:end);
end
