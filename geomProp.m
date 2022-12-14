function [Ad,A,C,S,E] = geomProp(stemImage,decayImage,imageType,...
    pixelExtent,varargin)
%{
DESCRIPTION
Using PiCUS sonic tomograms or binary images depicting decayed tree parts, 
this function computes several geometric properties of a measured tree 
part, including its cross-sectional area, A (m^2); cicularity, C 
(dimensionless); solidity, S (dimensionless); and eccentricity, 
E (dimensionless). The function also computes the total cross-sectional 
damaged area, Ad (%). 
NOTES
If analyzing PiCUS sonic tomograms, the function requires a stem image
showing the blue trunk boundary outline and a decay image showing the
visualized decay pattern. If analyzing binary images, the two images should 
contain white pixels corresponding to solid, undamaged wood. See the
related article for more details: 
Burcham, D.C., Brazee, N.J., Marra, R.E., Kane, B. 2023. Geometry matters
for sonic tomography of trees. 
[Ad,A,C,S,E] = geomProp(stemImage,decayImage,imageType,pixelExtent);
[Ad,A,C,S,E] = geomProp(stemImage,decayImage,imageType,pixelExtent,...
'colors',value);
REQUIRED INPUTS
stemImage: string - filepath for image file showing outline of the 
measured tree part
decayImage: string - filepath for image file showing internal decay
imageType: string - type of image used for analysis, either a tomogram 
('pit') or binary image ('bin')
pixelExtent: 2x2 matrix - spatial extent of pixels in the two images with
the first and second column entries denoting the spatial extent of pixels 
in the x- and y-directions, respectively, and the first and second row 
entries corresponding to the stem and decay image, respectively
OPTIONAL INPUTS
colors: string - colors used to select damaged parts, either green, violet,
and blue ('GVB'), violet and blue ('VB'), or blue ('B')
PARAMETER                     CLASS       DEFAULT VALUE
---------------------------------------------------------------------------
colors                        string      'GVB'
OUTPUTS
Ad: Percent (%) of total damaged cross sectional area.
A: Cross-sectional area (m^2)
C: Circularity (dimensionless)
S: solidity (dimensionless)
E: eccentricity (dimensionless)
Copyright 2022 Daniel C. Burcham
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
---------------------------------------------------------------------------
%}

format long

% Check function call
if nargin < 4
    error('geomProp requires at least four inputs');
elseif ~isempty(varargin)
    if mod(length(varargin),2)~=0
        if (length(varargin)==1 && isstruct(varargin{1}))
            varargin = reshape([fieldnames(varargin{1})...
                struct2cell(varargin{1})]',1,[]);
        else
            error(strcat('Inputs must be paired: geomProp(stemImage,',...
                'decayImage,imageType,pixelExtent,''PropertyName'',',...
                'PropertyValue,...)'));
        end
    elseif ~any(strcmp(varargin{find(strcmp(varargin,'colors'))+1},...
            {'GVB','VB','B'}))
        error('colors must be set equal to ''GVB'', ''VB'', or ''B''');
    end
end

% Default parameters
Colors = 'GVB';

% User-specified parameters
if ~isempty(varargin)
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'colors'
                Colors = varargin{i+1};
        otherwise
            error([varargin{i} 'is not a valid property for the geomProp function.']);
        end
    end
end
input1 = find(strcmp({'B','VB','GVB'},Colors))-1;

% Import image files (.jpg)
A = imread(stemImage,'jpg');
B = imread(decayImage,'jpg');
sz1 = size(A);
R1 = imref2d(sz1,pixelExtent(1,1),pixelExtent(1,2));
sz2 = size(B);
R2 = imref2d(sz2,pixelExtent(2,1),pixelExtent(2,2));

% Segment image
O = segment(A,B,R1,R2,input1,imageType);

% Search for open cavities and modify outer boundary
for i=1:size(O{2,1},1)-1
    [in,~]=inpolygon(O{2,1}{i+1,1}(:,1),O{2,1}{i+1,1}(:,2),...
        O{1,1}(:,1),O{1,1}(:,2));
    if any(~in)
        N = subtract(polyshape(O{2,1}{1,1}(:,1),O{2,1}{1,1}(:,2),...
            'Simplify',false,'KeepCollinearPoints',true),...
            polyshape(O{2,1}{i+1,1}(:,1),O{2,1}{i+1,1}(:,2),...
            'Simplify',false,'KeepCollinearPoints',true));
        O{2,1}{1,1} = [N.Vertices(:,1),N.Vertices(:,2)];
        O{2,1}{i+1,1}=[];
        clear N
    end
end
if isShapeMultipart(O{2,1}{1,1}(:,1),O{2,1}{1,1}(:,2))
    [X,Y]=polysplit(O{2,1}{1,1}(:,1),O{2,1}{1,1}(:,2));
    U=cell(sum(cellfun('length',X)>=75),1);
    ia=1;
    for i=1:length(X)
        if length(X{i,1}) > 75
            U{ia,1}=[X{i,1} Y{i,1}];
            ia=ia+1;
        end
    end
    O{2,1}{1,1}=U;
else
    U=cell(1,1);
    U{1,1} = O{2,1}{1,1};
    O{2,1}{1,1}=U;
end

% Reshape O to remove empty cells and do isShape
O{2,1} = O{2,1}(~cellfun('isempty',O{2,1}));
if iscell(O{2,1}{1,1})
    O{2,1}{1,1}=O{2,1}{1,1}(~cellfun('isempty',O{2,1}{1,1}));
end

% Geometric properties of solid area

% Area
A = polyarea(O{1,1}(:,1),O{1,1}(:,2));

% Circularity
dp = diff(O{1,1}([1:end 1], :), 1, 1);
L = sum(hypot(dp(:, 1), dp(:, 2)));
C = 4*pi*A/L^2;

% Eccentricity
V = eig(cov([O{1,1}(:,1),O{1,1}(:,2)]));
V = 1.414.*sqrt(V);
E = sqrt(1-(min(V)/max(V))^2);

% Solidity = Area/ConvexArea
idx = convhull(O{1,1}(:,1),O{1,1}(:,2));
Ar = polyarea(O{1,1}(idx,1),O{1,1}(idx,2));
S = A/Ar;

% Decayed area
D = zeros(size(O{2,1},1)+size(O{2,1}{1,1},1),1);
for i = 1:size(O{2,1},1)+size(O{2,1}{1,1},1)
    if i == 1
        D(i) = polyarea(O{1,1}(:,1),O{1,1}(:,2));
    elseif i >= 2 && i <= size(O{2,1}{1,1},1)+1
        D(i) = polyarea(O{2,1}{1,1}{i-1,1}(:,1),O{2,1}{1,1}{i-1,1}(:,2));
    elseif i > size(O{2,1}{1,1},1)+1
        D(i) = polyarea(O{2,1}{i-size(O{2,1}{1,1},1),1}(:,1),...
            O{2,1}{i-size(O{2,1}{1,1},1),1}(:,2));
    end
end
D(size(O{2,1}{1,1},1)+2:end) = -1.*D(size(O{2,1}{1,1},1)+2:end);
Ad = 1-sum(D(2:end)/D(1));
end