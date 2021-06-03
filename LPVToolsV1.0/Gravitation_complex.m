function G = Gravitation_complex(ecef)
%h An array of m altitudes, in meters, with respect to the WGS84 ellipsoid.
%lat An array of m geodetic latitudes, in degrees, where north latitude is positive, and south latitude is negative.
%lon An array of m geodetic longitudes, in degrees, where east longitude is positive, and west longitude is negative. This input is available only with method specified as 'CloseApprox'or'Exact'.
%jd A scalar value specifying Julian date used to calculate Julian Centuries from Epoch J2000.0. This input is available only with method specified as 'CloseApprox'or'Exact'.
g = 0;  %filtering out garbage values
gt = 0; %filtering out garbage values
gn = 0; %filtering out garbage values
jd = 0;
method = 'Exact';
noatm = false;
nocent = false;
prec = false;
action = 'Warning';
Warning = 0;
%longtitude doesn't need a parameter shift

%latitude might need a parameter shift!
lla = ecef2lla([ecef]');
%(latitude, longitude and altitude)
h = lla(end);
lat = lla(1);
lon = lla(2);

if h < 0 %check if altitude is lower than 0
    h = 0;
    Warning = 1;
end

g = gravitywgs84(h, lat, lon, method, [noatm, nocent, prec, jd], action); %The function

G = [g; gt; gn; Warning]; %The return values
        return;
%g An array of m gravity values in the direction normal to the Earth's surface at a specific lat lon location. 
%- A positive value indicates a downward direction.

%gt An array of m gravity values in the direction tangential to the Earth's surface at a specific lat lon location. 
%- A positive value indicates a northward direction. This option is available only with method specified as'Exact'.

%gn An array of m total gravity values in the direction normal to the Earth's surface at a specific lat lon location. 
%- A positive value indicates a downward direction. This option is available only with method specified as'Exact'.

% Warning, if 0 then no warning. if 1 then the altitude has been
% aproximated to 0 because it was below zero!
end