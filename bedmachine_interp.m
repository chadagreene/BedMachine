function zi = bedmachine_interp(variable,lati_or_xi,loni_or_yi,varargin)
% bedmachine_interp interpolates BedMachine data to any coordinates.
% The data are from Morlighem et al.'s BedMachine datasets.
% 
% http://sites.uci.edu/morlighem/dataproducts/
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
%  zi = bedmachine_interp(variable,lati,loni) 
%  zi = bedmachine_interp(variable,xi,yi)
%  zi = bedmachine_interp(...,IceSheet)
%  zi = bedmachine_interp(...,'datum',datum)
%  zi = bedmachine_interp(...,'method',InterpMethod)
% 
%% Description 
% 
% zi = bedmachine_interp(variable,lati,loni) interpolates a specified 
% BedMachine variable to the given by geo coordinates lati,loni. The 
% variable can be: 
%    * 'mask'      0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice, 4 = non-Greenland land or Vostok
%    * 'surface'   meters relative to EIGEN-EC4 geoid.
%    * 'thickness' meters
%    * 'bed'       meters relative to EIGEN-EC4 geoid.
%    * 'errbed'    meters 
%    * 'firn'      meters firn air thickness (must be added to surface or thickness for true)
%    * 'source'    Greenland: 0 = none, 1 = gimpdem, 2 = Mass conservation, 3 = synthetic, 4 = interpolation, 5 = hydrostatic equilibrium, 6 = kriging, 7 = RTOPO-2, 8 = gravity inversion, 9 = Millan et al. 2021, 10+ = bathymetry data)
%                  Antarctic: 1 = REMA/IBCSOv2, 2 = Mass conservation, 3 =interpolation, 4 = hydrostatic, 5 = streamline diffusion, 6 = gravity inversion, 7=seismic, 10=multibeam  
%    * 'geoid'     meters above WGS84 ellipsoid
%    * 'dataid'    1=GIMPdem or REMA, 2=Radar, 7=seismic, 10=multibeam
%    * 'base'      meters base of the ice sheet (bottom of ice shelves, but same as bed over grounded ice.) 
%    * 'wct'       meters water column thickness (derived, not an official BedMachine product.) 
%    * 'taf'       meters thickness above flotation (derived, not an official BedMachine product.) 
%    * 'flex'      dimensionless coefficient of tidal flexure (can slightly exceed 1; see Vaughan 1995 or Holdsworth 1969; derived, not an official BedMachine product; Requires Image Processing Toolbox. ) 
%    * 'head'      meters freshwater equivalent, static pressure head (derived, not an official BedMachine product.)  
%
% zi = bedmachine_interp(variable,xi,yi) As above, but for polar stereographic coordinates
% xi,yi in meters (ps70 for Greenland; ps71 for Antarctica). The function automatically 
% determines whether input coordinates are geo or polar stereographic via the islatlon function. 
% 
% zi = bedmachine_interp(...,IceSheet) specifies either 'greenland' or 'antarctica' (default). 
% 
% zi = bedmachine_interp(...,'datum',datum) specifies a datum as either 'geoid' (default) 
% or 'ellipsoid' for wgs84. 
% 
% zi = bedmachine_interp(...,'method',InterpMethod) specifies any interpolation method allowed by 
% the interp2 function. 
% 
%% Examples
% For examples, type 
%
%   showdemo bedmachine_interp_documentation  
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
% See also bedmachine_data and bedmachine_profile. 

%% Initial error checks: 

narginchk(3,Inf)
assert(~isnumeric(variable),'Error: variable must be a string, e.g. ''bed'', ''surface'', ''mask'', etc') 
assert(isnumeric(lati_or_xi),'Error: Input coordinates must be numeric.') 
assert(isequal(size(lati_or_xi),size(loni_or_yi)),'Error: Dimensions of input coordinates must agree.') 

%% Set Defaults: 

IceSheet = 'antarctica'; 
datum = 'geoid'; 

% Interpolation method: 
if ismember(variable,{'mask','source','dataid'})
   method = 'nearest'; 
else
   method = 'bilinear'; 
end

%% Parse inputs: 

% Are inputs georeferenced coordinates or polar stereographic?
if islatlon(lati_or_xi,loni_or_yi)
   % Check hemisphere: 
   if any(lati_or_xi(:)>0)
      [xi,yi] = ll2psn(lati_or_xi,loni_or_yi); % The ll2psn function is part of Arctic Mapping Tools package, the lesser known sibling of Antarctic Mapping Tools. 
      IceSheet = 'greenland'; 
   else
      IceSheet = 'antarctica'; % This might be declared explicitly by the user later, but just in case they forget, this will do what they probably want to do. 
      [xi,yi] = ll2ps(lati_or_xi,loni_or_yi); % The ll2ps function is in the Antarctic Mapping Tools package.
   end
else 
   xi = lati_or_xi;
   yi = loni_or_yi;    
end

% Which ice sheet? (Default is already set to Antarctica unless input coordinates have negative latitudes)
if any(strncmpi(varargin,'antarctica',3))
   IceSheet = 'antarctica';  
elseif any(strncmpi(varargin,'greenland',4))
   IceSheet = 'greenland';    
end
   

% Check datum: 
tmp = strncmpi(varargin,'datum',3); 
if any(tmp)
   datum = varargin{find(tmp)+1}; 
end

%% Load data: 

[Z,x,y] = bedmachine_data(variable,xi,yi,...
   'buffer',2,IceSheet,'datum',datum,'xy'); 

% mask out ocean for thickness or surface (to prevent interpolating at the 0 nan boundary) 
if ismember(variable(1:3),{'thi','sur'})
   mask = bedmachine_data('mask',xi,yi,'buffer',2,IceSheet,'xy'); 
   Z(mask==0) = nan; 
end
%% Interpolate: 

% Use FastInterp (a subfunction below) only if 'linear','bilinear', or 'nearest' AND number of query points is fewer than 10^4. 
if ismember(lower(method(1:3)),{'lin','bil','nea'}) & numel(xi)<10^4
   zi = FastInterp(x,y,Z,xi,yi,method); 
else
   zi = interp2(x,y,Z,xi,yi,method);
end

end


function zi = FastInterp(x,y,z,xi,yi,method)
% FastInterp is a subfunction that performs bilinear or nearest-neighbor
% interpolation on 2D datasets. In many cases, FastInterp is faster than
% interp2. However, if the number of query points is large (>10^4), you
% may find that interp2 is faster. Note that the number of *sample points*
% can be large with FastInterp, but the number of *query points* cannot. 
% 
% This function was written by Mathieu Morlighem (mmorligh@uci.edu)
% (inspired by Nathaniel Brahms) and was subsequently slightly modified by 
% Chad A. Greene in October 2018. 

%get data size
[M, N] = size(z);

% Get X and Y library array spacing
ndx = 1/(x(2)-x(1));    
ndy = 1/(y(2)-y(1));

% Begin mapping xi and yi vectors onto index space by subtracting library
% array minima and scaling to index spacing
xi = (xi - x(1))*ndx;       
yi = (yi - y(1))*ndy;

% Fill Zi with NaNs
zi = NaN(size(xi));

switch lower(method(1:3))
   case 'nea' % allows 'near' or 'nearest'

   % Find the nearest point in index space:
   rxi = round(xi)+1;  
   ryi = round(yi)+1;
   
   % Find points that are in X,Y range:
   flag = rxi>0 & rxi<=N & ~isnan(rxi) & ryi>0 & ryi<=M & ~isnan(ryi);
   
   % Map subscripts to indices:
   ind = ryi + M*(rxi-1);
   zi(flag) = z(ind(flag));

   case {'lin','bil'} % allows 'linear', 'lin', or 'bilinear'
   % Transform to unit square
   fxi = floor(xi)+1;  
   fyi = floor(yi)+1; % x_i and y_i
   
   dfxi = xi-fxi+1;    
   dfyi = yi-fyi+1;   % Location in unit square

   % flagIn determines whether the requested location is inside of the data arrays
   flagIn = fxi>0 & fxi<N & ~isnan(fxi) & fyi>0 & fyi<M & ~isnan(fyi);

   %Toss all out-of-bounds variables now to save time
   fxi  = fxi(flagIn);  
   fyi  = fyi(flagIn);
   dfxi = dfxi(flagIn); 
   dfyi = dfyi(flagIn);

   %Find bounding vertices
   ind1 = fyi + M*(fxi-1);     % indices of (  x_i  ,  y_i  )
   ind2 = fyi + M*fxi;         % indices of ( x_i+1 ,  y_i  )
   ind3 = fyi + 1 + M*fxi;     % indices of ( x_i+1 , y_i+1 )
   ind4 = fyi + 1 + M*(fxi-1); % indices of (  x_i  , y_i+1 )

   % Bilinear interpolation
   zi(flagIn) = ...
   z(ind1).*(1-dfxi).*(1-dfyi) + ...
   z(ind2).*dfxi.*(1-dfyi) + ...
   z(ind4).*(1-dfxi).*dfyi + ...
   z(ind3).*dfxi.*dfyi;

   otherwise
      error('Unrecognized interpolation method. Must be ''nearest'' or ''bilinear'' for FastInterp.') 
end

end

