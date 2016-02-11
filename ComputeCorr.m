function varargout = ComputeCorr(varargin)
%COMPUTECORR M-file for ComputeCorr.fig
%      COMPUTECORR, by itself, creates a new COMPUTECORR or raises the existing
%      singleton*.
%
%      H = COMPUTECORR returns the handle to a new COMPUTECORR or the handle to
%      the existing singleton*.
%
%      COMPUTECORR('Property','Value',...) creates a new COMPUTECORR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ComputeCorr_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      COMPUTECORR('CALLBACK') and COMPUTECORR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in COMPUTECORR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ComputeCorr

% Last Modified by GUIDE v2.5 11-Feb-2016 01:17:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ComputeCorr_OpeningFcn, ...
                   'gui_OutputFcn',  @ComputeCorr_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before ComputeCorr is made visible.
function ComputeCorr_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ComputeCorr
handles.output = hObject;

handles.rawim=[];
handles.segim=[];
handles.radius=0;
handles.final_im=[];
% Update handles structure
if strcmp(get(hObject,'Visible'),'off')
    axes(handles.imcell);
    set(gca,'xtick',[],'ytick',[])
    axes(handles.imbacksub);
    set(gca,'xtick',[],'ytick',[])
    cla;
end
guidata(hObject, handles);
% UIWAIT makes ComputeCorr wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ComputeCorr_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function ComputeCorr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ComputeCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function nclass_Callback(hObject, eventdata, handles)
% hObject    handle to nclass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nclass as text
%        str2double(get(hObject,'String')) returns contents of nclass as a double


% --- Executes during object creation, after setting all properties.
function nclass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nclass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function loadspot_Callback(hObject, eventdata, handles)
% hObject    handle to loadspot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rnafile, pathname]= uigetfile({'*.loc;*locx;*txt','Quantification file'; '*.*','All (*.*)'}, 'Pick the mRNA file');
% data should be in the following format:  x y intensity
if ~ (isequal(rnafile,0) || isequal(pathname,0))
    rnaf = fullfile(pathname, rnafile);
    handles.rna = load(rnaf);
end

guidata(hObject, handles)

% --------------------------------------------------------------------
function raw_Callback(hObject, eventdata, handles)
% hObject    handle to raw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, pathname]= uigetfile({'*.tif;*.tiff','Microscopy image file'; '*.*','All (*.*)'}, 'Pick the Cell file');
if ~ (isequal(file,0) || isequal(pathname,0))
    rawfile = fullfile(pathname, file);
    handles.rawim = imread(rawfile);
    curview = strcmp(get(get(handles.view,'SelectedObject'),'Tag'), 'rawview');
    if(curview == 1)
        axes(handles.imcell);
        imshow(handles.rawim,[]);
    end

    handles.final_im = handles.rawim;
    axes(handles.imbacksub);
    imshow(handles.final_im,[]);
end
guidata(hObject, handles);


% --- Executes on slider movement.
function ballsize_Callback(hObject, eventdata, handles)
% hObject    handle to ballsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

coeff = str2double(get(hObject, 'String'));
if isnan(coeff)
    coeff = handles.radius;
end
curview = strcmp(get(get(handles.backsub,'SelectedObject'),'Tag'), 'rollball');
if(curview==1)
    background=imopen(handles.rawim, strel('disk', round(coeff)));
    handles.final_im=imsubtract(handles.rawim,background);
end
axes(handles.imbacksub);
imshow(handles.final_im,[]);
%imshow(mat2gray(imsubtract(handles.rawim, mean(handles.rawim(:)))));
%image( find(abs(image-background) <= threshold) ) = 0;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ballsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ballsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when selected object is changed in backsub.
function backsub_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in backsub 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
stype = get(hObject, 'Tag');

if strcmp(stype, 'outmean') &&  ~isempty(handles.segim) &&  ~isempty(handles.rawim)
    background = mean(handles.rawim(handles.segim==0));
    handles.final_im = imsubtract(handles.rawim,background);
    axes(handles.imbacksub);
    imshow(handles.final_im,[]);
end

guidata(hObject, handles);


function intensity_list = getDensity(handles)
cell_list = sort(unique(handles.segim(:)))';
cell_list = cell_list(cell_list~=0);
intensity_list = double(zeros(length(cell_list),7));
% format :  cell, cell_int, area, total_int, min_int, max_int, mean_int, density
for i=1:length(cell_list)
    bw = handles.segim==cell_list(i);
    area = regionprops(bw, 'Area');
    im = double(handles.final_im);
    total_int = double(sum(im(bw), 'double'));
    max_int = max(max(im(bw)));
    min_int = min(min(im(bw)));
    density = double(total_int / area.Area);
    intensity_list(i,:) = [i, double(cell_list(i)), area.Area, min_int, max_int, density, total_int];
end

% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    intensity_list = getDensity(handles);
    
    nclass = str2double(get(handles.nclass , 'String'));
    values = intensity_list(:,end-1);
    cluster = ones(size(values));
    no_err = true;
    if isnan(nclass)
        nclass = 2;
    end
    try
        cluster = ckmeans(values, nclass);
    catch
        error = errordlg('Try again!  if you got this error again, Please decrease the number of clusters', 'Error'); 
        no_err = false;
        uiwait(error)
    end

    if no_err
        fprintf('Number of clusters : %d\n', nclass);
        intensity_list(:, end+1) = cluster;
        cluster_list = zeros(nclass, 6);
        for j=1:nclass
            cluster_list(j,1:2) = [j, mean(values(cluster==j))];
        end
        
        trans_infos = zeros(length(intensity_list(:,2)), 3);
        int_vals = intensity_list(:, 2);
        for i=1:numel(int_vals)
            int_val = intensity_list(i);
            trans_infos(i, 1:2) = sum(handles.rna((handles.rna(:,4)==int_val & handles.rna(:,5)>0), 5:6), 1);
            trans_infos(i,3) = nnz(handles.rna(:,4)==int_val);
        end
        intensity_list(:, end+1:end+3) = trans_infos;
        
        for j=1:nclass
            cluster_list(j, 3) = nnz(intensity_list(:, 8)==j);
            cluster_list(j, 4:6) =  sum(intensity_list(intensity_list(:, 8)==j, 9:11));
        end
        
        % get cluster mean density  
        [filename, pathname] = uiputfile('intensity_by_nucleus.txt', 'Save Intensity');
        header={'Cell_num', 'Cell_int', 'Area', 'minInt','maxInt', 'densityInt', 'totalInt', 'cluster', 'TransSite', 'Nascents', 'totalSpots'};
        fid = fopen(fullfile(pathname, filename), 'w');
        if fid == -1; error('Cannot open file: %s', outfile); end
        fprintf(fid, '%s\t', header{:});
        fprintf(fid, '\n');
        fclose(fid);
        dlmwrite(fullfile(pathname, filename), intensity_list,'delimiter', '\t', '-append', 'precision','%.2f');
        
        % save information per cluster now
        [filename, pathname2] = uiputfile(fullfile(pathname, 'intensity_by_cluster.txt'), 'Save clusters');
        fid = fopen(fullfile(pathname2, filename), 'w');
        if fid == -1; error('Cannot open file: %s', outfile); end
        header={'Cluster', 'meanDensity', 'Nucleus', 'Transites','Nascents', 'totalSpots'};
        fid = fopen(fullfile(pathname2, filename), 'w');
        if fid == -1; error('Cannot open file: %s', outfile); end
        fprintf(fid, '%s\t', header{:});
        fprintf(fid, '\n');
        fclose(fid);
        dlmwrite(fullfile(pathname2, filename), cluster_list,'delimiter', '\t', '-append');
       
    end
 catch err
    wdlg = warndlg('Error, missing data, please load image or spots', 'Incomplete informations'); 
    disp(err.message)
    uiwait(wdlg)
 
 end


% --- Executes when selected object is changed in view.
function view_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in view 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'rawview'
        axes(handles.imcell);
        imshow(handles.rawim,[]);
    otherwise
        axes(handles.imcell);
        imshow(handles.segim,[]);
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function segmented_Callback(hObject, eventdata, handles)
% hObject    handle to segmented (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, pathname]= uigetfile({'*.tif;*.tiff','Microscopy image file'; '*.*','All (*.*)'}, 'Pick the Segmented file');
if ~ (isequal(file,0) || isequal(pathname,0))
    segfile = fullfile(pathname, file);
    finfo = imfinfo(segfile);
    mask = imread(segfile);
    if ~ (strcmp(finfo.ColorType, 'grayscale') && min(unique(mask))==0)
        mask = bwlabel(im2bw(mat2gray(mask),0),4);
    end
    handles.segim = mask;
    m= mode(handles.segim(:)');

    bw = handles.segim==m;
    majoraxis= regionprops(bw, 'MajorAxisLength');
    %minoraxis= regionprops(bw, 'MinorAxisLength');
    handles.radius=(majoraxis(1).MajorAxisLength)/2;
    curview = strcmp(get(get(handles.view,'SelectedObject'),'Tag'), 'rawview');
    if curview == 0
        axes(handles.imcell);
        imshow(handles.segim, []);
    end
end
guidata(hObject, handles);


% --- Executes on button press in guess.
function guess_Callback(hObject, eventdata, handles)
% hObject    handle to guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

intensity_list = getDensity(handles);
values = intensity_list(:,end-1);
total_case = numel(unique(values));
valid_index(values, total_case);


function valid_index(data, N)
clusters = {};
NC = 2:N;
Rd = 'Euclidean';
Sil = zeros(1,N);
DB = zeros(1,N);
CH = zeros(1,N);
KL = zeros(1,N);
Ha = zeros(1,N);
Re = strcmp(Rd, 'Euclidean');
% (2) Internal validity indices when true labels are unknown
for i = NC
    clusters{i} = ckmeans(data, i);
    R = silhouette(data, clusters{i}, Rd);
    Sil(i) = mean(R);        % average Silhouette
    % Davies-Bouldin, Calinski-Harabasz, Krzanowski-Lai
    [DB(i), CH(i), KL(i), Ha(i), ST] = valid_internal_deviation(data, clusters{i}, Re);
    S = ind2cluster(clusters{i});
    %[Hom(i), Sep(i), wtertra(i)] = valid_internal_intra(Dist, S, Re, dmax);
end

kl = KL(NC);
ha = Ha(NC);
nl = length(NC);
S = trace(ST);
kl = [S kl];
ha = [S ha];
R = abs(kl(1:nl)-kl(2:nl+1));
S = [R(2: end) R(end)];
kl = R./S;
kl(nl) = kl(nl-1);
R = ha(1:nl)./ha(2:nl+1);
ha = (R-1).*(1-[NC(1)-1 NC(1:nl-1)]-1); 
KL(NC) = kl;
Ha(NC) = ha;

% (3) plotting indices
SR = [Sil; DB; CH; KL; Ha];
kfind = [20 20 20 20 2 1 2 2 5 2 2]; 
FR = {'Silhouette (Sil) ==> Highest', 'Davies-Bouldin (DB) ==> Lowest', 'Calinski-Harabasz (CH) ==> Highest', 'Krzanowski-Lai (KL) ==> Highest', ...
     'Hartigan' };

valid_index_plot(SR(:,NC), NC, kfind, FR);

