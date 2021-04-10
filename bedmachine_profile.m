function [hice,hbed,hwater] = bedmachine_profile(lati_or_xi,loni_or_yi,varargin)
% bedmachine_profile plots a 2D profile of ice, water, and rock elevations
% along any path in Greenland or Antarctica. The data are from Morlighem et
% al.'s BedMachine datasets.
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
%  bedmachine_profile(lati,loni)
%  bedmachine_profile(xi,yi)
%  bedmachine_profile(...,IceSheet)
%  bedmachine_profile(...,'horiz',HorizontalAxisValues)
%  bedmachine_profile(...,PatchProperty,PatchValue)
%  bedmachine_profile(...,'wgs84')
%  [hice,hbed,hwater] = bedmachine_profile(...)
% 
%% Description
% 
% bedmachine_profile(lati,loni) plots a side-view profile along a path given by
% geo coordinates lat,lon. lat and lon must be 1D arrays of equal length. If only two
% points are entered, an equally-spaced 1000-point line is created between those points. 
% 
% bedmachine_profile(xi,yi) as above, but for polar stereographic meters xi,yi.
% For Greenland, xi,yi are ps70. For Antarctica, it's ps71. 
% 
% bedmachine_profile(...,IceSheet) specifies the ice sheet as 'greenland'  
% or 'antarctica' (default). 
%
% bedmachine_profile(...,'horiz',HorizontalAxisValues) specifies horizontal axis
% values where HorizontalAxisValues is a 1D monotonically-increasing or decreasing
% array of dimensions corresponding to lat and lon. By default,
% HorizontalAxisValues are calculated as pathdistps. If you prefer to
% plot profiles with respect to some other values such as latitude of a north/south
% transect, use bedmachine_profile(lat,lon,'horiz',lat). 
%
% bedmachine_profile(...,PatchProperty,PatchValue) specifies edge line width, face
% color, and edge color of ice, water, or bed. The following properties may
% be specified: 
% 
%     * 'IceFace',ColorSpec
%     * 'IceEdge',ColorSpec
%     * 'IceEdgeWidth',LineWidth
%     * 'WaterFace',ColorSpec
%     * 'WaterEdge',ColorSpec
%     * 'WaterEdgeWidth',LineWidth
%     * 'BedFace',ColorSpec
%     * 'BedEdge',ColorSpec
%     * 'BedEdgeWidth',LineWidth
%     * 'Sky',ColorSpec
% 
% bedmachine_profile(...,'wgs84') plots profiles relative to the WGS84
% ellipsoid. (Profiles are plotted relative to the EIGEN-EC4 geoid by default.) 
% Note: Ice surface, ice thickness, and bed elevations are plotted correctly
% when using the 'wgs84' tag; however, ocean surfaces always appear at zero, 
% which could be innacurate by up to 60 meters or so.
% 
% [hice,hbed,hwater] = bedmachine_profile(...) returns handles of ice, bed,
% and water patch objects. 
% 
%% Examples
% For examples, type 
%
%   showdemo bedmachine_profile_documentation  
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
% See also bedmachine_3d, bedmachine_interp, and bedmachine_data. 

%% Error checks: 

assert(isvector(lati_or_xi)==1,'Profile plotting requires that input coordinates in the form of a vector.'); 
assert(isvector(loni_or_yi)==1,'Profile plotting requires that input coordinates in the form of a vector.'); 
assert(isscalar(lati_or_xi)==0,'You must enter more than one point for a profile.') 
assert(numel(lati_or_xi)==numel(loni_or_yi),'Input lat vector must be the same size as input lon.');

%% Set defaults: 

datum = 'geoid'; 
IceSheet = 'antarctica'; 
distax = true; 

% Colors: 
iceface = [0.8431    1.0000    0.9961]; 
iceedge = [0.4549    0.5922    0.5882]; 
waterface = [0.2745    0.5529    0.6902]; 
wateredge = 'none'; 
bedface = [0 0 0];    
bededge = 'none'; %[0.1137 0.0078      0];  
sky = 'w'; 

% Line widths:
iceedgewidth = 1; 
wateredgewidth = 1; 
bededgewidth = 1; 

%% Parse user inputs: 

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
    
if nargin>2
    
    % Horizontal Axis Values: 
    tmp = strncmpi(varargin,'horizontalaxis',3); 
    if any(tmp)
        horizax = varargin{find(tmp)+1}; 
        distax = false; % don't use track distance as the horizontal axis. 
        
        % Columnate horizontal axis to ensure consistent behavior: 
        horizax = horizax(:); 
        assert(isnumeric(horizax)==1,'horizontal axis values must be numeric.')
        assert(numel(horizax)==numel(lati_or_xi),'If you enter values for a horizontal axis, they must correspond to input lat,lon values. It looks like you have entered a horizontal axis array that does not match the size of lat and lon.') 
    end
    
    % All other inputs are colors or Line widths:
    tmp = strcmpi(varargin,'iceface');
    if any(tmp)
        iceface = varargin{find(tmp)+1}; 
    end
    
    tmp = strcmpi(varargin,'bedface');
    if any(tmp)
        bedface = varargin{find(tmp)+1}; 
    end
    
    tmp = strcmpi(varargin,'waterface');
    if any(tmp)
        waterface = varargin{find(tmp)+1}; 
    end
    
    tmp = strcmpi(varargin,'iceedge');
    if any(tmp)
        iceedge = varargin{find(tmp)+1}; 
    end
    
    tmp = strcmpi(varargin,'bededge');
    if any(tmp)
        bededge = varargin{find(tmp)+1}; 
    end
    
    tmp = strcmpi(varargin,'wateredge');
    if any(tmp)
        wateredge = varargin{find(tmp)+1}; 
    end
    
    tmp = strncmpi(varargin,'iceedgewidth',8);
    if any(tmp)
        iceedgewidth = varargin{find(tmp)+1}; 
        assert(isscalar(iceedgewidth)==1,'iceedgewidth value must be scalar.')
    end
    
    tmp = strncmpi(varargin,'bededgewidth',8);
    if any(tmp)
        bededgewidth = varargin{find(tmp)+1}; 
        assert(isscalar(bededgewidth)==1,'iceedgewidth value must be scalar.')
    end
    
    tmp = strncmpi(varargin,'wateredgewidth',10);
    if any(tmp)
        wateredgewidth = varargin{find(tmp)+1}; 
        assert(isscalar(wateredgewidth)==1,'iceedgewidth value must be scalar.')
    end
    
    tmp = strcmpi(varargin,'sky');
    if any(tmp)
        sky = varargin{find(tmp)+1}; 
    end
    
    % Which ice sheet? (Default is already set to Antarctica unless input coordinates have negative latitudes)
    if any(strncmpi(varargin,'greenland',5))
       IceSheet = 'greenland';    
    end
 
    % Check datum: 
    tmp = strncmpi(varargin,'datum',3); 
    if any(tmp)
       datum = varargin{find(tmp)+1}; 
    end
end

%% Reformat coordinates as needed:  

% Get rid of nans: 
isf = isfinite(xi) & isfinite(yi); 
xi = xi(isf); 
yi = yi(isf); 

% If only two input points, turn them into a 1000 point line: 
if numel(xi)==2
    xi = linspace(xi(1),xi(2),1000); 
    yi = linspace(yi(1),yi(2),1000); 
end

% Columnate for consistency: 
xi = xi(:); 
yi = yi(:); 

if distax
   switch lower(IceSheet(1:3))
      case 'gre'
         horizax = pathdistpsn(xi,yi,'km'); 
      case 'ant'
         horizax = pathdistps(xi,yi,'km'); 
      otherwise
         error(['Unrecognized ice sheet ',IceSheet'])
   end
else
   horizax = horizax(isf); % ensures same elements as input coordinates 
end
   

%% Get data

% Get bedmachine data: 
sfz = bedmachine_interp('surface',xi,yi,'datum',datum,IceSheet); 
bed = bedmachine_interp('bed',xi,yi,'datum',datum,IceSheet); 
thck = bedmachine_interp('thickness',xi,yi,'datum',datum,IceSheet); 

sfz(isnan(sfz))=0; 
thck(isnan(thck)) = 0; 
icebottom = sfz-thck;

% Indices of non-ocean data: 
% nreal = isfinite(icebottom);

% Some vertical limits: 
maxsfz = max(sfz(isfinite(sfz))); 
if isempty(maxsfz)
    maxsfz = 0; 
end
minbed = min(bed(isfinite(bed))); 
padding = (maxsfz-minbed)/20; 

%% Generate plot 

% Draw water:
hwater=fill([horizax(1);horizax(end);horizax(end);horizax(1)],[0;0;minbed-padding;minbed-padding],waterface);
set(hwater,'edgecolor',wateredge,'linewidth',wateredgewidth)
hold on;

% Draw bed:
realbed = find(isfinite(bed)); 
hbed = fill([horizax(realbed);horizax(realbed(length(realbed)));horizax(realbed(1))],[bed(realbed);minbed-padding;minbed-padding],bedface);
set(hbed,'edgecolor',bededge,'linewidth',bededgewidth)

% Draw ice:
hice = patch([horizax;flipud(horizax)],[sfz;flipud(icebottom)],iceface); 
set(hice,'edgecolor',iceedge,'linewidth',iceedgewidth);

% Format axes:
axis tight; 
box off;
if ismember(lower(datum(1:3)),{'geo','eig'})
   ylabel('elevation (m w.r.t EIGEN-EC4)')
else
   ylabel('elevation (m w.r.t WGS84)')
end

% Only label x axis if it's a distance axis: 
if distax
    xlabel('distance along profile (km)')
end

% Sky color
set(gca,'color',sky)
set(gcf,'color',sky)

%% Clean Up

if nargout==0
    clear hice
end
