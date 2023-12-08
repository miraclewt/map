function [xEast, yNorth, zUp] = aer2enu(az, elev, slantRange, angleUnit)
%AER2ENU Local spherical AER to Cartesian ENU
%
%   [xEast, yNorth, zUp] = AER2ENU(az, elev, slantRange) transforms point
%   locations in 3-D from local spherical coordinates (azimuth angle,
%   elevation angle, slantRange) to local Cartesian coordinates (xEast,
%   yNorth, zUp). The input angles are assumed to be in degrees.
%
%   [...] = AER2ENU(..., angleUnit) uses angleUnit, which matches either
%   'degrees' or 'radians', to specify the units of the azimuth and
%   elevation angles.
%
%   See ENU2AER for a description of the east-north-up (ENU) local 3-D
%   Cartesian coordinate system, and for definitions of the azimuth angle,
%   elevation angle, and slant range.
%
%   The transformation is similar to SPH2CART, except that the input angles
%   are in degrees by default and the first input angle is measured
%   clockwise from the positive Y-axis rather than counterclockwise from
%   the positive X-axis.
%
%   Class support for inputs az, elev, slantRange:
%      float: double, single
%
%   See also AER2NED, ENU2AER

% Copyright 2012-2019 The MathWorks, Inc.

if nargin < 4 || map.geodesy.isDegree(angleUnit)
    sinfun = @sind;
    cosfun = @cosd;
else
    sinfun = @sin;
    cosfun = @cos;
end

[xEast, yNorth, zUp] = aer2enuFormula(az, elev, slantRange, sinfun, cosfun);


function [xEast, yNorth, zUp] = aer2enuFormula( ...
    azimuthAngle, elevationAngle, slantRange, sinfun, cosfun)
% Transform local spherical (AER convention) to Cartesian (ENU convention)

zUp = slantRange .* sinfun(elevationAngle);
r   = slantRange .* cosfun(elevationAngle);
xEast  = r .* sinfun(azimuthAngle);
yNorth = r .* cosfun(azimuthAngle);
