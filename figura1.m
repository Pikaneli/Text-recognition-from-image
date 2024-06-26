function varargout = figura1(varargin)
% FIGURA1 MATLAB code for figura1.fig
%      FIGURA1, by itself, creates a new FIGURA1 or raises the existing
%      singleton*.
%
%      H = FIGURA1 returns the handle to a new FIGURA1 or the handle to
%      the existing singleton*.
%
%      FIGURA1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIGURA1.M with the given input arguments.
%
%      FIGURA1('Property','Value',...) creates a new FIGURA1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before figura1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to figura1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help figura1

% Last Modified by GUIDE v2.5 16-May-2024 17:09:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @figura1_OpeningFcn, ...
                   'gui_OutputFcn',  @figura1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before figura1 is made visible.
function figura1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to figura1 (see VARARGIN)

% Choose default command line output for figura1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes figura1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = figura1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ucitaj.
function ucitaj_Callback(hObject, eventdata, handles)
% hObject    handle to ucitaj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp', 'Slike (*.jpg,*.png,*.bmp)'}, 'Odaberi sliku');

% Provera da li je odabrana slika ili je dijalog za izbor datoteke zatvoren bez odabira
if isequal(filename,0) || isequal(pathname,0)
    disp('Odabir slike je otkazan.');
else
    % Učitavanje odabrane slike
    img = imread(fullfile(pathname, filename));
    
    % Prikaz odabrane slike u odgovarajućem axes komponenti GUI-ja
    axes(handles.axes1);
    imshow(img);
    title('Originalna slika');
    handles.img = img;
    axes(handles.axes2);
    title('')
    cla;
    set(handles.text1, 'String', '');
end
guidata(hObject, handles);


% --- Executes on button press in min_filtar.
function kontrast1_Callback(hObject, eventdata, handles)
% hObject    handle to min_filtar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.img;
E = 3; m = 100;   
r=double(img);
J=1./(1+ ( m./(r+eps) ).^E );
axes(handles.axes2);
J=im2uint8(J);
imshow(J);
title('Slika nakon poboljsanog kontrasta');
%set(handles.text1, 'String', '');
handles.slika2=J;
guidata(hObject, handles);

% --- Executes on button press in min_filtar.
function min_filtar_Callback(hObject, eventdata, handles)
% hObject    handle to min_filtar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.img;
img_gray = rgb2gray(img);
J = ordfilt2(img_gray, 1, ones(3,3));
axes(handles.axes2);
imshow(J);
title('Slika nakon min filtra');
%set(handles.text1, 'String', '');
handles.slika2=J;
guidata(hObject, handles);


% --- Executes on button press in sobel.
function sobel_Callback(hObject, eventdata, handles)
% hObject    handle to sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.img;
gx = imfilter(img, [-1 0 1; -2 0 2; -1 0 1], 'symmetric'); 

gy = imfilter(img, [-1 -2 -1; 0 0 0; 1 2 1], 'symmetric');  

g = abs(gx)+abs(gy);

axes(handles.axes2);
imshow(g,[]);
hold on
BW = edge(rgb2gray(img),'sobel');
imshow(BW);
title('Sobelov gradijent');
hold off
%set(handles.text1, 'String', '');
handles.slika2=g;
guidata(hObject, handles);

% --- Executes on button press in laplasijan1.
function laplasijan1_Callback(hObject, eventdata, handles)
% hObject    handle to laplasijan1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.img;
img=im2double(img);
w = [0 1 0; 1 -4 1; 0 1 0]; 

J = imfilter(img, w, 'replicate');  

J1 = img - J;        
J2 = mat2gray(J1);
low_high = stretchlim(J2);                              
J3 = imadjust(J2, low_high);
axes(handles.axes2);
imshow(J3);
title('Izostrena slika (centar -4)');
%set(handles.text1, 'String', '');
handles.slika2=J3;
guidata(hObject, handles);

% --- Executes on button press in laplasijan2.
function laplasijan2_Callback(hObject, eventdata, handles)  %furijeova transformacija
% hObject    handle to laplasijan2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.img;
img=rgb2gray(img);
[M N] = size(img);
[x y] = meshgrid(1:N,1:M);   
I = (-1).^(x + y);                                   
FT = fft2(double(img).*I);
axes(handles.axes2);
imshow(log(abs(FT)),[]);
title('Furijeova transformacija');
%set(handles.text1, 'String', '');
handles.slika2=FT;
guidata(hObject, handles);

% --- Executes on button press in max_filtar.
function max_filtar_Callback(hObject, eventdata, handles)   %morfoloski bliska slika
% hObject    handle to max_filtar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.img;
img=rgb2gray(img);

x=imbinarize(img);

x=~x;
g=strel('disk',5);
x=imclose(x,g);
axes(handles.axes2);
imshow(x);
title('Morfoloski bliska slika');
%set(handles.text1, 'String', '');
handles.slika2=x;
guidata(hObject, handles);


% --- Executes on button press in generisi.
function generisi_Callback(hObject, eventdata, handles)
% hObject    handle to generisi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%slika2 = handles.slika2;
slika1 = handles.img;

I = im2gray(slika1);

% Detect MSER regions.
[mserRegions, mserConnComp] = detectMSERFeatures(I, ... 
    "RegionAreaRange",[100 10000],"ThresholdDelta",4);

mserStats = regionprops(mserConnComp, "BoundingBox", "Eccentricity", ...
    "Solidity", "Extent", "Euler", "Image");

% Compute the aspect ratio using bounding box data.
bbox = vertcat(mserStats.BoundingBox);
w = bbox(:,3);
h = bbox(:,4);
aspectRatio = w./h;

% Threshold the data to determine which regions to remove. These thresholds
% may need to be tuned for other images.
filterIdx = aspectRatio' > 1.5; 
filterIdx = filterIdx | [mserStats.Eccentricity] > .8 ;
filterIdx = filterIdx | [mserStats.Solidity] < .5;
filterIdx = filterIdx | [mserStats.Extent] < 0.4 | [mserStats.Extent] > 0.7;
filterIdx = filterIdx | [mserStats.EulerNumber] < -3;

% Remove regions
mserStats(filterIdx) = [];
mserRegions(filterIdx) = [];

regionImage = mserStats(6).Image;
regionImage = padarray(regionImage, [1 1]);

% Compute the stroke width image.
distanceImage = bwdist(~regionImage); 
skeletonImage = bwmorph(regionImage, "thin", inf);

strokeWidthImage = distanceImage;
strokeWidthImage(~skeletonImage) = 0;

strokeWidthValues = distanceImage(skeletonImage);   
strokeWidthMetric = std(strokeWidthValues)/mean(strokeWidthValues);

strokeWidthThreshold = 0.4;
strokeWidthFilterIdx = strokeWidthMetric > strokeWidthThreshold;

for j = 1:numel(mserStats)
    
    regionImage = mserStats(j).Image;
    regionImage = padarray(regionImage, [1 1], 0);
    
    distanceImage = bwdist(~regionImage);
    skeletonImage = bwmorph(regionImage, "thin", inf);
    
    strokeWidthValues = distanceImage(skeletonImage);
    
    strokeWidthMetric = std(strokeWidthValues)/mean(strokeWidthValues);
    
    strokeWidthFilterIdx(j) = strokeWidthMetric > strokeWidthThreshold;
    
end

% Remove regions based on the stroke width variation
mserRegions(strokeWidthFilterIdx) = [];
mserStats(strokeWidthFilterIdx) = [];

bboxes = vertcat(mserStats.BoundingBox);

% Convert from the [x y width height] bounding box format to the [xmin ymin
% xmax ymax] format for convenience.
xmin = bboxes(:,1);
ymin = bboxes(:,2);
xmax = xmin + bboxes(:,3) - 1;
ymax = ymin + bboxes(:,4) - 1;

% Expand the bounding boxes by a small amount.
expansionAmount = 0.02;
xmin = (1-expansionAmount) * xmin;
ymin = (1-expansionAmount) * ymin;
xmax = (1+expansionAmount) * xmax;
ymax = (1+expansionAmount) * ymax;

% Clip the bounding boxes to be within the image bounds
xmin = max(xmin, 1);
ymin = max(ymin, 1);
xmax = min(xmax, size(I,2));
ymax = min(ymax, size(I,1));

% Show the expanded bounding boxes
expandedBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];

% Compute the overlap ratio
overlapRatio = bboxOverlapRatio(expandedBBoxes, expandedBBoxes);

% Set the overlap ratio between a bounding box and itself to zero to
% simplify the graph representation.
n = size(overlapRatio,1); 
overlapRatio(1:n+1:n^2) = 0;

% Create the graph
g = graph(overlapRatio);

% Find the connected text regions within the graph
componentIndices = conncomp(g);

% Merge the boxes based on the minimum and maximum dimensions.
xmin = accumarray(componentIndices', xmin, [], @min);
ymin = accumarray(componentIndices', ymin, [], @min);
xmax = accumarray(componentIndices', xmax, [], @max);
ymax = accumarray(componentIndices', ymax, [], @max);

% Compose the merged bounding boxes using the [x y width height] format.
textBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];

% Remove bounding boxes that only contain one text region
numRegionsInGroup = histcounts(componentIndices);
textBBoxes(numRegionsInGroup == 1, :) = [];

% Show the final text detection result.
ITextRegion = insertShape(slika1, "rectangle", textBBoxes,"LineWidth",3);
allText = '';
set(handles.text1, 'String', '');

% Perform OCR on each detected text region and accumulate results.
for i = 1:size(textBBoxes, 1)
    bbox = textBBoxes(i, :);
    region = imcrop(I, bbox);
    ocrResult = ocr(region);
    allText = sprintf('%s%s', allText, ocrResult.Text);
    set(handles.text1, 'String', allText);
end

axes(handles.axes1);
imshow(ITextRegion);
title('Detektovan text');
guidata(hObject, handles);

% --- Executes on button press in sacuvaj.
function sacuvaj_Callback(hObject, eventdata, handles)
% hObject    handle to sacuvaj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = questdlg('Želite li da sačuvate editovanu sliku?');
if strcmp(answer, 'Yes') == 1
    if isempty(get(handles.axes2, 'Children'))
        msgbox('Prvo izvršite promene nad slikom.', 'Upozorenje', 'warn');
    else
        % Korisnik izabire ime i lokaciju za čuvanje slike
        [filename, pathname] = uiputfile({'*.jpg', 'JPEG slika (*.jpg)'}, 'Sačuvaj sliku kao');
        
        % Provera da li je korisnik odabrao ime i lokaciju za čuvanje
        if isequal(filename,0) || isequal(pathname,0)
            disp('Čuvanje slike je otkazano.');
        else
            filepath = fullfile(pathname, filename);
            slika = handles.slika2;
            imwrite(slika, filepath, 'jpg');
            msgbox('Podaci su sačuvani.', 'Sačuvano');
        end
    end
end
guidata(hObject, handles);

% --- Executes when selected object is changed in binarizacija.
function binarizacija_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in binarizacija 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.img;
radio_2_5 = get(handles.radiobutton1, 'Value');
radio_3 = get(handles.radiobutton2, 'Value');
radio_3_5 = get(handles.radiobutton3, 'Value');
radio_4 = get(handles.radiobutton4, 'Value');
radio_4_5 = get(handles.radiobutton5, 'Value');
radio_5 = get(handles.radiobutton6, 'Value');
radio_5_5 = get(handles.radiobutton7, 'Value');
radio_6 = get(handles.radiobutton8, 'Value');
prag=0;
if radio_2_5==1
    prag=0.25;
end
if radio_3==1
    prag=0.3;
end
if radio_3_5==1
    prag=0.35;
end
if radio_4==1
    prag=0.4;
end
if radio_4_5==1
    prag=0.45;
end
if radio_5==1
    prag=0.5;
end
if radio_5_5==1
    prag=0.55;
end
if radio_6==1
    prag=0.6;
end
BW2=im2bw(img,prag);
axes(handles.axes2);
imshow(BW2);
title('Slika nakon binarizacije');
%set(handles.text1, 'String', '');
handles.slika2=BW2;
guidata(hObject, handles);


% --- Executes on button press in iseci.
function iseci_Callback(hObject, eventdata, handles)
% hObject    handle to iseci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.img;
[y, x] = ginput(2);
img=img(min(round(x)):max(round(x)),min(round(y)):max(round(y)),:);
axes(handles.axes1);
imshow(img);
title('Isecena slika');
%set(handles.text1, 'String', '');
handles.slika2=img;
handles.img=img;
guidata(hObject, handles);
