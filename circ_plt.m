function [xp,yp] = circ_plt(r,ctr,n)
%CIRC_PLT  Given the radius or radii and center or centers of circles,
%          CIRC_PLT calculates the coordinates for one or more circles
%          for plotting.
%
%          [XP,YP] = CIRC_PLT(R,CTR) Given the radius or radii (R)
%          and center or centers (CTR) of a circle or multiple circles,
%          CIRC_PLT calculates the coordinates XP and YP for plotting
%          the circles.  The center coordinate matrix CTR should have
%          two columns for the X and Y coordinates.
%
%          [XP,YP] = CIRC_PLT(R,CTR,N) Uses N points to represent the
%          cicles.  N must be >=10.  The default is 144 points (every
%          2.5 degrees).
%
%          NOTES:  1.  For two-dimensional (2-D) circles.
%
%          08-Mar-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<3||isempty(n)
  n = 144;
end
%
if n<10
  n = 144;
end
%
if nargin<2
  error(' *** ERROR in CIRC_PLT:  Two inputs are required!');
end
%
[ncirc,ncoord] = size(ctr);
if ncoord~=2
  error([' *** ERROR in CIRC_PLT:  Input circle center array CTR must', ...
        ' have two columns!']);
end
%
r = r(:);
nr = size(r,1);
%
if nr~=ncirc
  error([' *** ERROR in CIRC_PLT:  Number of radii and centers must', ...
        ' be the same!']);
end
%
% Calculate Circle Coordinates
%
t = (0:n)/n;            % Fraction of circumference
t = t*2*pi;             % Angles in radians
st = sin(t);
st(n+1) = 0;
xp = r.*repmat(cos(t),ncirc,1);
xp = xp+repmat(ctr(:,1),1,n+1);
yp = r.*repmat(st,ncirc,1);
yp = yp+repmat(ctr(:,2),1,n+1);
if ncirc>1
  xp = [xp NaN(ncirc,1)];
  yp = [yp NaN(ncirc,1)];
end
xp = xp';
xp = xp(:);
yp = yp';
yp = yp(:);
%
return