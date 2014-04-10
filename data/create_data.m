function varargout = create_data(varargin)
% CREATE_DATA M-file for create_data.fig
%      CREATE_DATA, by itself, creates a new CREATE_DATA or raises the existing
%      singleton*.
%
%      H = CREATE_DATA returns the handle to a new CREATE_DATA or the handle to
%      the existing singleton*.
%
%      CREATE_DATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CREATE_DATA.M with the given input arguments.
%
%      CREATE_DATA('Property','Value',...) creates a new CREATE_DATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before create_data_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to create_data_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help create_data

% Last Modified by GUIDE v2.5 26-Apr-2007 01:11:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @create_data_OpeningFcn, ...
    'gui_OutputFcn',  @create_data_OutputFcn, ...
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


% --- Executes just before create_data is made visible.
function create_data_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to create_data (see VARARGIN)

% Choose default command line output for create_data

clc

cd /home/bar747/matlab/gamma_i/data

%z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/data']; eval(cmd)

set(0,'DefaultFigureColor',[0.961 0.988 0.965])

handles.output = hObject;

handles.grdx = '1:1:nx';

handles.grdy = '1:1:ny';

handles.grdz = '1:1:nz';

handles.select = 'ocean==5 & lats>=10';

load handles_last, handles.open_path = plast, clear plast

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes create_data wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = create_data_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function grid_res_longs_edit_Callback(hObject, eventdata, handles)
% hObject    handle to grid_res_longs_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_res_longs_edit as text
%        str2double(get(hObject,'String')) returns contents of grid_res_longs_edit as a double

handles.grdx = get(hObject,'String');

% Update handles structure
guidata(hObject, handles)



% --- Executes during object creation, after setting all properties.
function grid_res_longs_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grid_res_longs_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

handles.grdy = get(hObject,'String')

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to select_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of select_edit as text
%        str2double(get(hObject,'String')) returns contents of select_edit as a double

handles.select = get(hObject,'String');

%handles.select = [handles.select, ' & n>2']

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function select_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in select_ok_pushbutton.
function select_ok_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to select_ok_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global s ct t p g longs lats ocean n longss dist

if length(size(g))~=length(size(s)) | size(g)~=size(s)
    g = nan*ones(size(s));
end

[nz,ny,nx] = size(s);

longs0 = longs; lats0 = lats; longss = nanmean(longs), latss = nanmean(lats');

inds = find(isnan(s(1,:,:))); longs0(inds) = nan; lats0(inds) = nan;

whos

cmd = ['inds_o = find(', handles.select, ');'], eval(cmd)

%cmd = ['[longs_min,inds_o] = min(abs(longss-', handles.select, ')); longss(inds_o)'], eval(cmd)       %%%%%%%%   the inclusion criteria

if length(inds_o)==0
    error('****    no data matching selection criteria    ****')
end

inds = setdiff(1:nx*ny,inds_o); %[length(inds),nx*ny];

s(:,inds) = nan; t(:,inds) = nan; p(:,inds) = nan; g(:,inds) = nan;

longs0(inds) = nan; lats0(inds) = nan; ocean(inds) = nan; n(inds) = nan;

longs_min = nanmin(longs0(:)); longs_max = nanmax(longs0(:));

lats_min = nanmin(lats0(:)); lats_max = nanmax(lats0(:));

indsx = find(longs_min<=longss&longss<=longs_max);

indsy = find(lats_min<=latss&latss<=lats_max);

nx = length(indsx); ny = length(indsy);

s = s(:,indsy,indsx); t = t(:,indsy,indsx);

p = p(:,indsy,indsx);

if exist('g'), g = g(:,indsy,indsx); end
if exist('ocean'), ocean = ocean(indsy,indsx); end
if exist('n'), n = n(indsy,indsx); end
if exist('longs'), longs = longs(indsy,indsx); end
if exist('lats'), lats = lats(indsy,indsx); end

%%      squeeze a section if selected

s = squeeze(s); t = squeeze(t); p = squeeze(p); g = squeeze(g);

longs = squeeze(longs); lats = squeeze(lats);

figure

if length(size(s))==2
    longs_range = max(longs)-min(longs);
    lats_range = max(lats)-min(lats);
    set(gcf,'Name','ocean depth index','NumberTitle','off','Color',[0.961 0.988 0.965])
    if longs_range<lats_range
        longs = median(longs)*ones(size(lats));
        fpcolor(lats,(1:nz),t(:,:))
        xlabel('latitudes')
    elseif longs_range>lats_range
        lats = median(lats)*ones(size(longs));
        fpcolor(longs,(1:nz),t(:,:))
        xlabel('longitudes')
    else
        error('****    ERROR 1 in create_data.m    ****')
    end
    set(gca,'ydir','reverse','color','k')
elseif length(size(s))==3
    longss = nanmean(longs); latss = nanmean(lats'); dj_pltmp(longss,latss,squeeze(t(1,:,:)))
else
    error('****    ERROR 2 in create_data.m')
end



% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear global, clc, dj_tic

global s t p g ocean n longs lats

cmd = ['cd ', handles.open_path]; eval(cmd)

[FileName,PathName] = uigetfile('*.mat');

cmd = ['load ', PathName, FileName]; eval(cmd)

handles.open_path = PathName;

cd /home/bar747/matlab/gamma_i/data

%z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/data']; eval(cmd)

% update handles structure

guidata(hObject, handles);

% save file location

plast = handles.open_path; save handles_last plast

% adjust pressure if needs be

dims = size(p);
if length(dims)==2 & dims(2)==1
    if length(size(s))==2
        [nz,nx] = size(s);
        p = reshape(p(:)*ones(1,nx),nz,nx);
    else
        [nz,ny,nx] = size(s);
        p = reshape(p(:)*ones(1,nx*ny),nz,ny,nx);
    end
end

% adjust longs and lats if needs be

[nz,ny,nx] = size(s); dims = size(longs);

if length(size(s))==2
    
else
    if dims(1)~=ny | dims(2)~=nx
        longs = ones(ny,1)*longs(:)'; lats = lats(:)*ones(1,nx);
    end
end

% adjust ocean and n if needs be

% ocean = []

if exist('ocean')~=1 | length(ocean)==0 | exist('n')~=1 | length(n)==0
    disp( 'computing oceans and depth ...')
    dj_tic, oceans_and_n, dj_toc, dj_pause(1)
end

% figure, dj_pltmp(nanmean(longs),nanmean(lats')',ocean)

ok = '... loaded data'

return



% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cd ../webstuff
open create_data_css.html
cd ../data



% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s ct t p g ocean n longs lats

%global s1 t1

%s = s1; t = t1;

cmd = ['cd ', handles.open_path]; eval(cmd)

[file,path] = uiputfile('*.mat','save data as');

cmd = ['save ', path, file, ' s ct t p g ocean n longs lats '], eval(cmd)

plast = handles.open_path;

cd /home/bar747/matlab/gamma_i/data
% z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/data']; eval(cmd)

save handles_last plast


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double



% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z1_Callback(hObject, eventdata, handles)
% hObject    handle to grid_res_longs_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_res_longs_edit as text
%        str2double(get(hObject,'String')) returns contents of grid_res_longs_edit as a double



function select_edit_Callback(hObject, eventdata, handles)
% hObject    handle to select_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of select_edit as text
%        str2double(get(hObject,'String')) returns contents of select_edit as a double


handles.select = get(hObject,'String')

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in create_data_pushbutton.
function create_data_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to create_data_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s ct t p g longs lats ocean n

check_data

% replicate a section

longs_mean = nanmean(longs(:))

dims = size(t);
if length(dims)==2 && dims(2)>1
    longs_range = max(longs)-min(longs);
    lats_range = max(lats)-min(lats);
    if longs_range<lats_range
        xx = lats
    elseif longs_range>lats_range
        xx = longs
    else
        error(' ERROR 1 in create_data_pushbutton')
    end
    pp = nanmean(p(:,:)');
    figure, whos
    fpcolor(xx,pp,squeeze(t(:,:)))
    set(gca,'ydir','reverse','color',[0 0 0]), %figure(gcf)
    [nz,ny] = size(s);
    s = reshape(repmat(s,1,3),[nz,ny,3] ); t = reshape(repmat(t,1,3),[nz,ny,3] );
    p = reshape(repmat(p,1,3),[nz,ny,3] );
    longs_mean = mean(longs(:)); longs = [longs_mean-0.1; longs_mean; longs_mean+0.1];
    ocean = reshape(repmat(ocean(:),1,3),[ny,3] ); n = reshape(repmat(n(:),1,3),[ny,3] );
end

% adjust longs and lats if needs be

dims = size(longs)
if sum(size(longs)==size(n))~=2
    [longs,lats] = meshgrid(longs,lats);
end

% compute ct if needs be

if length(ct)==0 | sum(size(s)==size(t))~=sum(size(ct)==size(t))
    dj_disp('computing ct ...')
    inds = find(isfinite(t));
    ss = s(inds); tt = t(inds); pp = p(inds);
    ctt = ct_from_t(ss,tt,pp);
    ct = nan*ones(size(t)); ct(inds) = ctt;
end

whos('global'),  dj_toc



function grid_res_depth_edit_Callback(hObject, eventdata, handles)
% hObject    handle to grid_res_depth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_res_depth_edit as text
%        str2double(get(hObject,'String')) returns contents of grid_res_depth_edit as a double



% --- Executes on button press in grid_res_ok_pushbutton.
function grid_res_ok_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to grid_res_ok_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global s t p g ocean n longs lats

[nz,ny,nx] = size(s)

eval(['ix = ', handles.grdx, ',']); eval(['jy = ', handles.grdy, ',']);

s = s(:,jy,ix); t = t(:,jy,ix); p = p(:,jy,ix);

if exist('g'), g = g(:,jy,ix); end
if exist('ocean'), ocean = ocean(jy,ix); end
if exist('n'), n = n(jy,ix); end
if exist('longs'), longs = longs(jy,ix); end
if exist('lats'), lats = lats(jy,ix); end

whos

