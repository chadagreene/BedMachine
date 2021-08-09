function varargout = bedmachine_data(variable,varargin) 
% bedmachine_data loads data from Morlighem et al.'s BedMachine datasets.
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
%  Z = bedmachine_data(variable)
%  Z = bedmachine_data(variable,lati,loni)
%  Z = bedmachine_data(variable,xi,yi) 
%  Z = bedmachine_data(...,'buffer',extrakm)
%  Z = bedmachine_data(...,IceSheet) 
%  Z = bedmachine_data(...,'datum',datum) 
%  [Z,x,y] = bedmachine_data(...) 
%  [Z,Lat,Lon] = bedmachine_data(...,'geo')
%  [x,y] = bedmachine_data(outline)
%  [lat,lon] = bedmachine_data(outline,'geo')
% 
%% Description 
% 
% Z = bedmachine_data(variable) loads a specified variable for the whole ice sheet, 
% with elevations relative to the EIGEN-EC4 geoid. The variable can be: 
%    * 'mask'      0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice, 4 = non-Greenland land or Vostok
%    * 'surface'   meters relative to EIGEN-EC4 geoid. (for true surface you must add 'firn')
%    * 'thickness' meters
%    * 'bed'       meters relative to EIGEN-EC4 geoid.
%    * 'errbed'    meters 
%    * 'firn'      meters firn air thickness (must be added to surface or thickness for true) 
%    * 'source'    Greenland: 0 = none, 1 = gimpdem, 2 = Mass conservation, 3 = synthetic, 4 = interpolation, 5 = hydrostatic equilibrium, 6 = kriging, 7 = RTOPO-2, 8 = gravity inversion, 10+ = bathymetry data)
%                  Antarctic: 1 = REMA/IBCSO, 2 = Mass conservation, 3 = interpolation, 4 = hydrostatic, 5 = Kriging, 6 = gravity inversion
%    * 'geoid'     meters above WGS84 ellipsoid
%    * 'base'      meters base of the ice sheet (bottom of ice shelves, but same as bed over grounded ice.) 
%    * 'wct'       meters water column thickness (derived, not an official BedMachine product.) 
%    * 'taf'       meters thickness above flotation (derived, not an official BedMachine product.) 
%    * 'flex'      dimensionless coefficient of tidal flexure (can slightly exceed 1; see Vaughan 1995 or Holdsworth 1969; derived, not an official BedMachine product; Requires Image Processing Toolbox. ) 
%    * 'head'      meters freshwater equivalent, static pressure head (derived, not an official BedMachine product.)  
%
% Z = bedmachine_data(variable,lati,loni) returns only enough BedMachine data to fully
% encompass a set of points given by geo coordinates lati,loni. This is a
% good way to save computer memory, only loading, analyzing, and plotting the
% data you need to work in a region of interest. 
% 
% Z = bedmachine_data(variable,xi,yi) As above, but for polar stereographic coordinates
% xi, yi in meters (ps70 for Greenland; ps71 for Antarctica). The function automatically 
% determines whether input coordinates are geo or polar stereographic via the islatlon function. 
% 
% Z = bedmachine_data(...,'buffer',extrakm) as above, but adds a buffer around the input
% coordinates. This option is useful for loading only the data in your region
% of interest, plus a little extra around the sides for good measure. 
% 
% Z = bedmachine_data(...,IceSheet) specifies either 'greenland' or  
% 'antarctica' (default). 
% 
% Z = bedmachine_data(...,'datum',datum) specifies either 'greenland' or 'antarctica' (default). 
% 
% [Z,x,y] = bedmachine_data(...) returns polar stereographic meters (ps70 for
% Greenland; ps71 for Antarctica) as 1d arrays x,y. 
% 
% [Z,Lat,Lon] = bedmachine_data(...,'geo') returns gridded geo coordinates corresponding to 
% each pixel in Z. 
% 
% [x,y] = bedmachine_data(outline) gives the (derived) polar stereographic coordinates of 
% any of these outlines: 
%    * 'gl'      grounding line
%    * 'coast'   coast line
%    * 'hl'      approximate hydrostatic line given by flex=0.99. 
% 
% [lat,lon] = bedmachine_data(outline,'geo') as above, but returns geographic 
% coordinates. 
% 
%% Examples
% For examples, type 
%
%   showdemo bedmachine_data_documentation  
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
% See also bedmachine_interp and bedmachine_profile. 

%% Initial error checks: 

narginchk(1,Inf)
assert(~isnumeric(variable),'Error: variable must be a string, e.g. ''bed'', ''surface'', ''mask'', etc') 

%% Set defaults: 

% NOTE: To update the dataset filename manually, change it here:

GreenlandFilename = 'BedMachineGreenland-2021-04-20.nc';
%AntarcticaFilename = 'BedMachineAntarctica_2019-11-05_v01.nc'; % the OG
AntarcticaFilename = 'BedMachineAntarctica_2020-07-15_v02.nc'; 

subset = false;  % use whole data set (not a regional subset) by default 
extrakm = 0;     % zero buffer by default
xyout = true;    % give xy arrays by default
IceSheet = 'antarctica'; % Change the default ice sheet here if you'd like. Can be 'greenland' or 'antarctica', but ensure to change it in bedmachine_interp too.  
datum = 'geoid'; % default is EIGEN-6C4

rho_ice = 917;   % ice density for thickness above flotation and head calculation
rho_sw = 1027;   % seawater density for thickness above flotation 
rho_fw = 1000;   % freshwater density for head calculation 

E = 0.88e+9;     % elastic modulus from vaughan1995tidal abstract, although marsh2014grounding found 1.4e9, (for ice flexure calculation).
nu = 0.3;        % Poisson's ratio from Table 1 of vaughan1995tidal, (for ice flexure calculation)
gravity = 9.81;  % gravitational acceleration m/s^2, (for ice flexure calculation)  

%% Parse inputs: 

if nargin>1
   
   % Check for subset based on input coordinates: 
   if isnumeric(varargin{1}) 
      subset = true; 
      lati_or_xi = varargin{1}; 
      loni_or_yi = varargin{2}; 
      
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

      % Add a buffer around the edges of the data:
      tmp = strncmpi(varargin,'buffer',3); 
      if any(tmp)
         extrakm = varargin{find(tmp)+1}; 
         assert(numel(extrakm)<3,'Error: buffer must be one or two elements, in kilometers.') 
      end
   end

   % Is the user requesting x and y outputs instead of default lat,lon grid? 
   if any(strcmpi(varargin,'geo')) 
      xyout = false; 
   end

   % Which ice sheet? (Default is already set to unless input coordinates have negative latitudes) 
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
   
end

%% Define filename: 

switch IceSheet 
   case 'greenland'
      filename = GreenlandFilename;
      
   case 'antarctica'
      filename = AntarcticaFilename; 
      
   otherwise
      error('Unrecognized ice sheet and I have no clue how we got here.') 
end

assert(exist(filename,'file')==2,['Error: cannot find ',filename,'. If you have already downloaded the data, make sure Matlab can find it. If you do not have the data, get it here: https://nsidc.org/data/idbmg4']) 

%% Load data 

x = double(ncread(filename,'x')); 
y = double(ncread(filename,'y')); 

if subset
   
   if isscalar(extrakm)
      extrakm = [extrakm extrakm]; 
   end
   
    % A crude manual fix for when a single xi,yi lies between pixels: 
    if isscalar(xi)
          extrakm = [max([extrakm(1) 1]) max([extrakm(2) 1])]; 
    end
    
    % Get xlimits (xl) and ylimits (yl) of input coordinates + buffer:
    xl = [min(xi(:))-extrakm(1)*1000 max(xi(:))+extrakm(1)*1000];
    yl = [min(yi(:))-extrakm(2)*1000 max(yi(:))+extrakm(2)*1000];
    
    % Region of rows and columns of pixels to read: 
    ci=find((y>=yl(1))&(y<=yl(2)));
    ri=find((x>=xl(1))&(x<=xl(2)));
else
    ci = 1:length(y); 
    ri = 1:length(x); 
end

% Load data: 
switch lower(variable)
   case 'base'
      th = double(ncread(filename,'thickness',[ri(1) ci(1)],[length(ri) length(ci)]));
      th(th==-9999) = NaN; 
      sfz = double(ncread(filename,'surface',[ri(1) ci(1)],[length(ri) length(ci)]));
      sfz(sfz==-9999) = NaN; 
      Z = sfz-th; 
      
   case 'wct'
      bed = double(ncread(filename,'bed',[ri(1) ci(1)],[length(ri) length(ci)]));
      bed(bed==-9999) = NaN; 
      th = double(ncread(filename,'thickness',[ri(1) ci(1)],[length(ri) length(ci)]));
      th(th==-9999) = NaN; 
      sfz = double(ncread(filename,'surface',[ri(1) ci(1)],[length(ri) length(ci)]));
      sfz(sfz==-9999) = NaN; 
      Z = sfz-th-bed; 
      Z(abs(Z)<1) = 0; 
      
   case 'taf'
      assert(exist('base2freeboard.m','file')==2,'Cannot find base2freeboard function. Make sure you have Antarctic Mapping Tools in a place where Matlab can find it.')
      
      bed = double(ncread(filename,'bed',[ri(1) ci(1)],[length(ri) length(ci)]));
      bed(bed==-9999) = NaN; 
      th = double(ncread(filename,'thickness',[ri(1) ci(1)],[length(ri) length(ci)]));
      th(th==-9999) = NaN; 
      
      freeboard = base2freeboard(bed,'rhow',rho_sw,'rhoi',rho_ice); % base2freeboard is a function in AMT
      freeboard(bed>=0) = bed(bed>=0); % Accounts for bedrock above sea level.

      Z = th - (freeboard-bed);
      Z(Z<0) = 0; 
      
   case 'flex'
      mask = ncread(filename,'mask',[ri(1) ci(1)],[length(ri) length(ci)]); 
      th = double(ncread(filename,'thickness',[ri(1) ci(1)],[length(ri) length(ci)]));
      th(th==-9999) = 0; 
      
      grounded = mask==1 | mask==2;

      % Calculate distance to nearest grounded grid cell (requires Image Toolbox):
      GridCellSize = median(abs(diff(x))); 
      dst2ground = double(bwdist(grounded)) * GridCellSize;

      % Calculate Eq 4 of Vaughan1995tidal, (originally from Holdsworth 1969):
      beta = (3*rho_sw * gravity*((1-nu^2)./(E.*th.^3))).^(1/4);
      Z = 1 - exp(-beta.*dst2ground) .* (cos(beta.*dst2ground)+sin(beta.*dst2ground));
      Z(mask==1) = 0; % takes care of NaN output
      
   case 'head' 
      bed = double(ncread(filename,'bed',[ri(1) ci(1)],[length(ri) length(ci)]));
      bed(bed==-9999) = NaN; 
      th = double(ncread(filename,'thickness',[ri(1) ci(1)],[length(ri) length(ci)]));
      th(th==-9999) = NaN; 
      mask = ncread(filename,'mask',[ri(1) ci(1)],[length(ri) length(ci)]); 

      % Calculate head: 
      Z = th.*rho_ice./rho_fw + bed; 

      Z(mask==0) = NaN; % ocean 
      Z(mask==3) = NaN; % floating ice 

   case 'gl'
      assert(nargout==2,'Error: For grounding line, outputs must be [lat_or_x,lon_or_y] = bedmachine_data(''gl'')')
      
      % Grounding line was derived by: 
      % g = mask==1 | mask==2 | mask==4; 
      % [C,h] = contour(x,y,double(g),[0.5 0.5],'k');
      % [x1,y1] = C2xyz; 
      % gl = polyshape(x1,y1); 
      
      P = load('BedMachine_Antarctica_v1_outlines.mat','gl'); 
      x = P.gl.Vertices(:,1); 
      y = P.gl.Vertices(:,2); 
      
   case 'coast'
      assert(nargout==2,'Error: For coast line, outputs must be [lat_or_x,lon_or_y] = bedmachine_data(''coast'')')
      
      P = load('BedMachine_Antarctica_v1_outlines.mat','coast'); 
      x = P.coast.Vertices(:,1); 
      y = P.coast.Vertices(:,2); 
      
   case 'hl'
      assert(nargout==2,'Error: For hydrostatic line, outputs must be [lat_or_x,lon_or_y] = bedmachine_data(''hl'')')
      
      % Hydrostatic line was derived as: 
      % [F,x,y] = bedmachine_data('flex');
      % [C3,h3] = contour(x,y,F,[0.99 0.99],'k');
      % [x3,y3] = C2xyz; 
      % g3 = polyshape(x3,y3); 
      
      P = load('BedMachine_Antarctica_v1_outlines.mat','hl'); 
      x = P.hl.Vertices(:,1); 
      y = P.hl.Vertices(:,2); 
      
   otherwise 
      try
         Z = double(ncread(filename,variable,[ri(1) ci(1)],[length(ri) length(ci)]));
      catch
         error('Unable to load data. If the variable is correct, perhaps the spatial extents are the problem?') 
      end
      % Take care of NaNs: 
      Z(Z==-9999) = NaN; 
end

% Convert from geoid to ellipsoid reference if user requested it:
if ismember(variable,{'bed','surface'})
   switch lower(datum(1:3))
      case {'ell','wgs'} % accepts ellipsoid or wgs-84
         Z = Z + double(ncread(filename,'geoid',[ri(1) ci(1)],[length(ri) length(ci)]));
      case {'eig','geo'} % accepts eigen-6c4 or geoid
         % no further work required
      otherwise
         error('Unrecognized datum. Must be ''ellipsoid'' or ''geoid''.')
   end
end
   
%% Final adjustments for the export: 

% Just the grid
if nargout==1
   % Orient it correctly: 
   varargout{1} = flipud(rot90(Z));
end

% Continent outline only: 
if nargout==2
   if xyout
      varargout{1} = x; 
      varargout{2} = y; 
   else
      [varargout{1},varargout{2}] = ps2ll(x,y); 
   end
end

% Grid and coordinates: 
if nargout==3
   % Orient it correctly: 
   varargout{1} = flipud(rot90(Z));
   
   if xyout
      varargout{2} = x(ri); 
      varargout{3} = y(ci); 
   else
      
      % Meshgridding a whole continent of high res data might make the computer stop working, so warn the user for large datasets 
      if (length(ri)*length(ci))>1e7
         answer = questdlg('Warning: Gridding the geo coordinates of an area this large could slow your computer to a crawl. You may prefer to cancel and try again using the ''xy'' option. Do you wish to cancel?',...
            'Memory Warning',...
            'Go for it anyway','Cancel','Cancel'); 
         if strcmp(answer,'Cancel')
            return
         end
      end
         
      % Grid the points so we can get lat,lon coordinates of each grid point:  
      [X,Y] = meshgrid(x(ri),y(ci));
      
      if strcmpi(IceSheet,'greenland')
         [varargout{2},varargout{3}] = psn2ll(X,Y); 
      else
         [varargout{2},varargout{3}] = ps2ll(X,Y); 
      end
   end
end

end