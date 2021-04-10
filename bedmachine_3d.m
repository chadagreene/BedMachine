function h = bedmachine_3d(lati_or_xi,loni_or_yi,varargin) 
% bedmachine_3d creates a 3D map of BedMachine data. 
% 
%% Requirements 
% This function requires a set of Matlab tools and a Bedmachine dataset, 
% and both will depend on where you're working. Get them here:  
%
% * For Greenland: * 
% Arctic Mapping Tools: https://www.mathworks.com/matlabcentral/fileexchange/63324
% Greenland Bedmachine Data: https://nsidc.org/data/IDBMG4
% 
% * For Antarctica: *
% Antarctic Mapping Tools: https://www.mathworks.com/matlabcentral/fileexchange/47638
% Antarctic Bedmachine Data: https://nsidc.org/data/nsidc-0756
% 
%% Syntax 
% 
%  bedmachine_3d(latlim,lonlim)
%  bedmachine_3d(xlim,ylim)
%  bedmachine_3d(...,'buffer',extrakm)
%  bedmachine_3d(...,'alpha',iceAlpha)
%  bedmachine_3d(...,'greenland')
%  h = bedmachine_3d(...)
% 
%% Description 
%
% bedmachine_3d(latlim,lonlim) creates a 3D map of the region specified by 
% latitude and longitude limits. latlim and lonlin can be two-element arrays 
% to specify just the min and max extents, or latlim,lonlim can be a range 
% of values, and the resulting map will cover the entire range. Plotting might 
% be slow for very large areas, so start small and see what your computer can 
% handle. 
%
% bedmachine_3d(xlim,ylim) as above, but map extents are specified by polar 
% stereographic meters. 
%
% bedmachine_3d(...,'buffer',extrakm) adds a buffer of specified width in 
% kilometers around the region of interest. 
%
% bedmachine_3d(...,'alpha',iceAlpha) sets the transparency of the ice. By default, 
% iceAlpha=1, meaning fully opaque. Use a value between 0 and 1 for semitransparent
% ice. 
%
% bedmachine_3d(...,'greenland') plots Greenland instead of the default Antarctica. 
%
% h = bedmachine_3d(...) returns a structure of all the object graphical handles. 
%
%% Examples 
% For examples with real pretty pictures, type 
% 
%  showdemo bedmachine_3d_documentation 
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
% This function was written by Chad A. Greene, October 2018. 
% 
% See also bedmachine_profile and bedmachine_data. 


%% Input checks: 

assert(nargin>=2,'Plotting the entire continent might be slow, so you must define a region of interest.') 
assert(isequal(size(lati_or_xi),size(loni_or_yi)),'Dimensions of input coordinates must agree.') 
assert(license('test','image_toolbox')==1,'Sorry, this function requires the image processing toolbox.') 

% Set defaults: 
iceAlpha = 1; % transparency of ice. 0=invisible; 1=opaque. 
buf = 0; 
IceSheet = 'antarctica'; 

if nargin>2
   tmp = strncmpi(varargin,'alpha',3); 
   if any(tmp)
      iceAlpha = varargin{find(tmp)+1}; 
      assert(iceAlpha>=0 & iceAlpha<=1,'Ice alpha value must be between 0 and 1.') 
   end
   
   tmp = strncmpi(varargin,'buffer',3); 
   if any(tmp)
      buf = varargin{find(tmp)+1}; 
      assert(isscalar(buf) & buf<5000,'The buffer value must be a scalar in units of kilometers.') 
   end
   
   if any(strncmpi(varargin,'greenland',4))
      IceSheet = 'greenland';    
   end
end
   
%% Load data: 

[bed,x,y] = bedmachine_data('bed',lati_or_xi,loni_or_yi,'buf',buf,IceSheet); 
sfz = bedmachine_data('surface',lati_or_xi,loni_or_yi,'buf',buf,IceSheet); 
base = bedmachine_data('base',lati_or_xi,loni_or_yi,'buf',buf,IceSheet); 
mask = bedmachine_data('mask',lati_or_xi,loni_or_yi,'buf',buf,IceSheet); 

[X,Y] = meshgrid(x,y); 

%% Get the geometry of this 

L = bwlabel(mask~=0); 
P = regionprops(mask~=0);
ind = find([P.Area]==max([P.Area])); 
icesheet = L==ind; 

% Boundaries of the ice sheet: 
B = bwboundaries(icesheet); 
ind = sub2ind(size(X),B{1}(:,1),B{1}(:,2)); 
xi = X(ind); 
yi = Y(ind); 
sfzi = sfz(ind); 
basei = base(ind); 

% Boundaries of the bed: 
B = bwboundaries(true(size(X))); 
ind = sub2ind(size(X),B{1}(:,1),B{1}(:,2)); 
xb = X(ind); 
yb = Y(ind); 
bedb = bed(ind);

% Boundaries of the ocean: 
B = bwboundaries(imdilate(mask==0,strel('disk',1))); 
ind = sub2ind(size(X),B{1}(:,1),B{1}(:,2)); 
xo = X(ind); 
yo = Y(ind); 
bedo = bed(ind); 

base(mask~=3) = nan; 
sfz(mask==0) = nan; % ocean
sfz(mask==1) = nan; % rock

%% Plot the 3D map:

hold on
view(3) 
daspect([1 1 1/30])
clear h 

% Plot bed:
h.bed = surface(x,y,bed); 
shading interp
axis tight off

% Set the colormap for the bed: 
if exist('cmocean.m','file')
   cmocean topo 
   caxis([-1 1]*max(abs(bed(:))))
else
   disp 'I cannot find the cmocean colormap function, which is part of the Climate Data Toolbox. That''s okay; we will use Matlab''s default colormap.'
end

% Plot base:
h.base = surface(x,y,base); 
h.base.EdgeColor = 'none'; 
h.base.FaceColor = [0.8431    1.0000    0.9961]; 
h.base.FaceAlpha = iceAlpha; 

% Plot ice surface: 
h.sfz = surface(x,y,sfz); 
h.sfz.EdgeColor = 'none'; 
h.sfz.FaceColor = [0.8431    1.0000    0.9961]; 
h.sfz.FaceAlpha = iceAlpha; 

% Plot ice side: 
h.iceside = surface([xi xi],[yi yi],[basei sfzi]); 
h.iceside.EdgeColor = 'none'; 
h.iceside.FaceColor = [0.8431    1.0000    0.9961]; 
h.iceside.FaceAlpha = iceAlpha; 

% Plot bed side: 
h.bedside = surface([xb xb],[yb yb],[min(bed(:))*ones(size(bedb)) bedb]); 
h.bedside.EdgeColor = 'none'; 
h.bedside.FaceColor = 'k'; 

% plot ocean surface: 
h.oceansfz = patch(xo,yo,zeros(size(xo)),'b'); 
h.oceansfz.EdgeColor = [0.0118    0.4431    0.6118]; 
h.oceansfz.FaceColor = [0.0118    0.4431    0.6118]; 
h.oceansfz.FaceAlpha = 0.1; 

% Plot ocean side: 
h.oceanside = surface([xo xo],[yo yo],[bedo zeros(size(bedo))]); 
h.oceanside.EdgeColor = 'none'; 
h.oceanside.FaceColor = [0.0118    0.4431    0.6118]; 
h.oceanside.FaceAlpha = 0.1; 
h.oceansfz.SpecularStrength = 1;

% plot ocean corners: 
if bed(1,1)<0
   h.corner(1) = plot3(X(1,1)*[1 1],Y(1,1)*[1 1],[0 bed(1,1)],'color',[0.0118    0.4431    0.6118]); 
end
if bed(end,1)<0
   h.corner(2) = plot3(X(end,1)*[1 1],Y(end,1)*[1 1],[0 bed(end,1)],'color',[0.0118    0.4431    0.6118]); 
end
if bed(1,end)<0
   h.corner(3) = plot3(X(1,end)*[1 1],Y(1,end)*[1 1],[0 bed(1,end)],'color',[0.0118    0.4431    0.6118]); 
end
if bed(end,end)<0
   h.corner(4) = plot3(X(end,end)*[1 1],Y(end,end)*[1 1],[0 bed(end,end)],'color',[0.0118    0.4431    0.6118]); 
end

% Set lighting parameters:
material dull
camlight 
% lighting phong % phong is slow but perhaps slightly prettier 

%% Clean up: 

if nargout==0 
   clear h
end

end
