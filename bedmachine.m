function h = bedmachine(variable,varargin)
% bedmachine loads and plots bedmachine data on a polar stereographic map. 
% (Antarctica only)
% 
%% Syntax
% 
%  bedmachine
%  bedmachine(variable) 
%  bedmachine(variable,'contour',PropertyName,PropertyValue,...) 
%  bedmachine(outline,LineProperty,LineValue,...)
%  h = bedmachine(...) 
% 
%% Description 
% 
% bedmachine plots a gray grounding line and coast line from BedMachine. 
% 
% bedmachine(variable) plots any BedMachine variable as an imagesc plot.
% The variable can be: 
%    * 'gl'        grounding line
%    * 'coast'     coast line
%    * 'hl'        approximate hydrostatic line given by flex=0.99. 
%    * 'mask'      0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice, 4 = non-Greenland land or Vostok
%    * 'surface'   meters relative to EIGEN-EC4 geoid.
%    * 'thickness' meters
%    * 'bed'       meters relative to EIGEN-EC4 geoid.
%    * 'errbed'    meters 
%    * 'source'    Greenland: 0 = none, 1 = gimpdem, 2 = Mass conservation, 3 = synthetic, 4 = interpolation, 5 = hydrostatic equilibrium, 6 = kriging, 7 = RTOPO-2, 8 = gravity inversion, 10+ = bathymetry data)
%                  Antarctic: 1 = REMA/IBCSO, 2 = Mass conservation, 3 = interpolation, 4 = hydrostatic, 5 = Kriging, 6 = gravity inversion
%    * 'geoid'     meters above WGS84 ellipsoid
%    * 'base'      meters base of the ice sheet (bottom of ice shelves, but same as bed over grounded ice.) 
%    * 'wct'       meters water column thickness (derived, not an official BedMachine product.) 
%    * 'taf'       meters thickness above flotation (derived, not an official BedMachine product.) 
%    * 'flex'      dimensionless coefficient of tidal flexure (can slightly exceed 1; see Vaughan 1995 or Holdsworth 1969; derived, not an official BedMachine product; Requires Image Processing Toolbox. ) 
%  
% bedmachine(variable,'contour',PropertyName,PropertyValue,...) plots any 
% BedMachine variable as a contour plot. 
%
% bedmachine(outline,LineProperty,LineValue,...) plots any of the following, 
% with optional line formatting: 
%    * 'gl'      grounding line
%    * 'coast'   coast line
%    * 'hl'      approximate hydrostatic line given by flex=0.99. 
% 
% h = bedmachine(...) returns a handle h of the plotted object(s). 
% 
%% Examples 
% For examples, type 
% 
%   showdemo bedmachine_documentation 
%
%% Citations
% If you use BedMachine data, please cite the Morlighem paper listed below. 
% And if this function is useful for you, please do me a kindness and cite 
% my Antarctic Mapping Tools paper. 
% 
% Morlighem, M., E. Rignot, T. Binder, D. D. Blankenship, R. Drews, G. Eagles, 
% O. Eisen, F. Ferraccioli, R. Forsberg, P. Fretwell, V. Goel, J. S. Greenbaum,
% H. Gudmundsson, J. Guo, V. Helm, C. Hofstede, I. Howat, A. Humbert, W. Jokat,
% N. B. Karlsson, W. Lee, K. Matsuoka, R. Millan, J. Mouginot, J. Paden, F. Pattyn,
% J. L. Roberts, S. Rosier, A. Ruppel, H. Seroussi, E. C. Smith, D. Steinhage, 
% B. Sun, M. R. van den Broeke, T. van Ommen, M. van Wessem, and D. A. Young. 2019. 
% Deep glacial troughs and stabilizing ridges unveiled beneath the margins of the
% Antarctic ice sheet, Nature Geoscience. https://doi.org/10.1038/s41561-019-0510-8
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools  
% for Matlab. Computers & Geosciences. 104 (2017) pp.151-157. 
% http://dx.doi.org/10.1016/j.cageo.2016.08.003
% 
%% Author Info 
% This function was written by Chad A. Greene of NASA Jet Propulsion 
% Laboratory, December 2019. 
% 
% See also: bedmachine_data and bedmachine_interp. 

%% Shortcut to a simple map: 

if nargin==0
   h(1) = bedmachine('coast','color',[0.5725    0.5843    0.5686]); 
   h(2) = bedmachine('gl','color',[0.5725    0.5843    0.5686]);
   if nargout==0
      clear h 
   end
   return
end

%% Set defaults

outline = true; 
contourplot = false; 

%% Parse inputs

if nargin>0
   if ~ismember(lower(variable),{'gl','hl','coast'})
      outline = false; 
   end
   
   tmp = strncmpi(varargin,'contour',3); 
   if any(tmp)
      contourplot = true; 
      varargin = varargin(~tmp); 
   end
end

%% Get initial conditions: 

da = daspect; 
da = [1 1 da(3)]; 
hld = ishold; 
hold on

mapwasopen = ~isequal(axis,[0 1 0 1]); 
ax = axis; 

%%

if outline
   
   % Load bedmachine-derived outlines: 
   [x,y] = bedmachine_data(variable); 
   
   if mapwasopen
      OutOfBounds = x<ax(1) | x>ax(2) | y<ax(3) | y>ax(4); 
      x(OutOfBounds) = NaN; % trims away everything outside current map extents while keeping the nans that separate different sections of the outline 
      y(OutOfBounds) = NaN; 
   else
      axis([min(x) max(x) min(y) max(y)])
   end
   
   h = plot(x,y,varargin{:}); 
else
   
   if mapwasopen
      [Z,x,y] = bedmachine_data(variable,ax(1:2),ax(3:4)); 
   else
      [Z,x,y] = bedmachine_data(variable); 
   end
   
   if contourplot
      [~,h] = contour(x,y,Z,varargin{:}); 
   else
      
      h = imagesc(x,y,Z); 
      set(h,'alphadata',isfinite(Z)); 
      axis xy; 
   end
   
end

%% Put things back the way we found them: 

daspect(da)
if ~hld
   hold off
end

if mapwasopen
   axis(ax); 
else 
   axis([min(x) max(x) min(y) max(y)]); 
end

%% Clean up: 

if nargout==0 
   clear h
end

end