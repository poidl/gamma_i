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

%cd d:/neutrals/ness8/data
default_dir = which('create_data.m');
cd ([default_dir(1:end-14)])
%z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/data']; eval(cmd) 

set(0,'DefaultFigureColor',[0.961 0.988 0.965])

handles.output = hObject;
handles.grdx = '1:1:nx';
handles.grdy = '1:1:ny';
handles.grdz = '1:1:nz';
handles.select = 'ocean==5 & lat>=10';

try
    load handles_last
    handles.open_path = plast
    clear plast
catch
    handles.open_path = default_dir(1:end-14);
end

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


function grid_res_long_edit_Callback(hObject, eventdata, handles)
% hObject    handle to grid_res_long_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_res_long_edit as text
%        str2double(get(hObject,'String')) returns contents of grid_res_long_edit as a double

handles.grdx = get(hObject,'String');

% Update handles structure
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function grid_res_long_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grid_res_long_edit (see GCBO)
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


global SA CT t p g long lat ocean n longs lats dist

if length(size(g))~=length(size(SA)) | size(g)~=size(SA)
    g = nans(size(SA));
end

[nz,ny,nx] = size(SA);

long0 = long;
lat0 = lat;
longs = nanmean(long);
lats = nanmean(lat');

Inan = find(isnan(SA(1,:,:)));
long0(Inan) = nan;
lat0(Inan) = nan;

whos

cmd = ['inds_o = find(', handles.select, ');'], eval(cmd)   

%cmd = ['[long_min,inds_o] = min(abs(longs-', handles.select, ')); longs(inds_o)'], eval(cmd)       %%%%%%%%   the inclusion criteria

if length(inds_o)==0
    error('****    no data matching selection criteria    ****')
end

inds = setdiff(1:nx*ny,inds_o); %[length(inds),nx*ny];

SA(:,inds) = nan;
CT(:,inds) = nan;
%t(:,inds) = nan;
p(:,inds) = nan;
g(:,inds) = nan; 

long0(inds) = nan;
lat0(inds) = nan;
ocean(inds) = nan;
n(inds) = nan;

long_min = nanmin(long0(:));
long_max = nanmax(long0(:));

lat_min = nanmin(lat0(:));
lat_max = nanmax(lat0(:));

indsx = find(long_min<=longs & longs<=long_max);
indsy = find(lat_min<=lats & lats<=lat_max);

nx = length(indsx);
ny = length(indsy);

SA = SA(:,indsy,indsx);
CT = CT(:,indsy,indsx);
%t = t(:,indsy,indsx);
p = p(:,indsy,indsx); 

if exist('g','var'), g = g(:,indsy,indsx); end
if exist('ocean','var'), ocean = ocean(indsy,indsx); end
if exist('n','var'), n = n(indsy,indsx); end
if exist('long','var'), long = long(indsy,indsx); end
if exist('lat','var'), lat = lat(indsy,indsx); end
    
%      squeeze a section if selected

SA = squeeze(SA);
CT = squeeze(CT);
%t = squeeze(t);
p = squeeze(p);
g = squeeze(g);
long = squeeze(long);
lat = squeeze(lat); 

figure

if length(size(SA))==2
    long_range = max(long) - min(long);
    lat_range = max(lat) - min(lat);
    set(gcf,'Name','ocean depth index','NumberTitle','off','Color',[0.961 0.988 0.965])
    if long_range < lat_range
        long = median(long)*ones(size(lat));
        pcolor(lat,(1:nz),CT(:,:))
        xlabel('Latitude')
    elseif long_range>lat_range
        lat = median(lat)*ones(size(long));
        pcolor(long,(1:nz),CT(:,:))
        xlabel('Longitude')
    else
        error('****    ERROR 1 in create_data.m    ****')
    end
    set(gca,'ydir','reverse','color','k')
elseif length(size(SA))==3
	longs = nanmean(long);
    lats = nanmean(lat');
    dj_pltmp(longs,lats,squeeze(CT(1,:,:)))    
else
    error('****    ERROR 2 in create_data.m')
end


% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%clear global, clc, dj_tic

global SA CT t p g ocean n long lat
try
    %cmd = ['cd ', handles.open_path]; eval(cmd)
    cd ([handles.open_path])
catch
    default_dir = which('create_data.m');
    cd ([default_dir(1:end-14)])
end

[FileName,PathName] = uigetfile('*.mat');

cmd = ['load ', PathName, FileName]; eval(cmd)

handles.open_path = PathName;

default_dir = which('create_data.m');
cd ([default_dir(1:end-14)])

%cd d:/neutrals/ness8/data

%z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/data']; eval(cmd) 

% update handles structure

guidata(hObject, handles);

% save file location

plast = handles.open_path;
save handles_last plast 

% adjust pressure if needs be

dims = size(p);
if length(dims)==2 & dims(2)==1
  if length(size(SA))==2
    [nz,nx] = size(SA);
    p = reshape(p(:)*ones(1,nx),nz,nx); 
  else
    [nz,ny,nx] = size(SA);
    p = reshape(p(:)*ones(1,nx*ny),nz,ny,nx); 
  end
end
 
% adjust long and lat if needs be

[nz,ny,nx] = size(SA);
dims = size(long);

if length(size(SA))==2

else
  if dims(1)~=ny | dims(2)~=nx
    long = ones(ny,1)*long(:)';
    lat = lat(:)*ones(1,nx);
  end
end

% adjust ocean and n if needs be

% ocean = []

if exist('ocean','var')~=1 | length(ocean)==0 | exist('n','var')~=1 | length(n)==0
    disp( 'computing oceans and depth ...')
    dj_tic
    oceans_and_n
    dj_toc
    dj_pause(1)
end

% figure, dj_pltmp(nanmean(long),nanmean(lat')',ocean)

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

global SA CT t p g ocean n long lat

%global s1 t1 

%s = s1; t = t1;

cmd = ['cd ', handles.open_path]; eval(cmd)

[file,path] = uiputfile('*.mat','save data as');

cmd = ['save ', path, file, ' SA CT t p g ocean n long lat '], eval(cmd)

plast = handles.open_path;

%cd d:/neutrals/ness8/data
% z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/data']; eval(cmd) 
default_dir = which('create_data.m');
cd ([default_dir(1:end-14)])

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
% hObject    handle to grid_res_long_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_res_long_edit as text
%        str2double(get(hObject,'String')) returns contents of grid_res_long_edit as a double



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

global SA CT t p g long lat ocean n

%check_data
%--------------------------------------------------------------------------
% delete casts with internal missing data and less than three levels
%--------------------------------------------------------------------------                           
Icast = find(isfinite(SA(1,:)));
[nz, nxy] = size(SA);
ss = SA(:,Icast);
z = flipud(ss);
[nzz, nxx] = size(z);
zz = nans(1,nxx);
for i = 1:nxx
    indss = find(isfinite(z(:,i)));
    zz(i) = nz - indss(1);
end

% (i) delete shallow
Idelete = find(n(Icast) < zz);
%number_of_shallow_casts = length(Idelete);
%salt = SA(:,inds(Idelete));
SA(:,Icast(Idelete)) = nan;
CT(:,Icast(Idelete)) = nan;
%t(:,Icast(Idelete)) = nan;
  p(:,Icast(Idelete)) = nan;
g(:,Icast(Idelete)) = nan;
ocean(Icast(Idelete)) = nan;
n(Icast(Idelete)) = nan;

%(ii) delete less than 3 bottles
Idelete = find(n < 3); 
%number_of_casts_less_than_3_bottles = length(Idelete);
%salt = s(:,inds_delete)
SA(:,Idelete) = nan;
CT(:,Idelete) = nan;
%t(:,Idelete) = nan;
  p(:,Icast(Idelete)) = nan;
g(:,Idelete) = nan;
ocean(Idelete) = nan;
n(Idelete) = nan;

%--------------------------------------------------------------------------


% replicate a section
long_mean = nanmean(long(:))

dims = size(CT);
if (length(dims) == 2) && (dims(2) > 1)
    long_range = max(long) - min(long); 
    lat_range = max(lat) - min(lat);
    if long_range < lat_range
        xx = lat
    elseif long_range > lat_range
        xx = long
    else
        error(' ERROR 1 in create_data_pushbutton')
    end
    pp = nanmean(p(:,:)');
   
    whos
    
    figure
    fpcolor(xx,pp,squeeze(CT(:,:)))
    set(gca,'ydir','reverse','color',[0 0 0]), %figure(gcf)
    [nz,ny] = size(SA);
    SA = reshape(repmat(SA,1,3),[nz,ny,3] );
    CT = reshape(repmat(CT,1,3),[nz,ny,3] ); 
    %t = reshape(repmat(t,1,3),[nz,ny,3] ); 
    p = reshape(repmat(p,1,3),[nz,ny,3] ); 
    long_mean = mean(long(:));
    long = [long_mean-0.1; long_mean; long_mean+0.1];
    ocean = reshape(repmat(ocean(:),1,3),[ny,3] );
    n = reshape(repmat(n(:),1,3),[ny,3] );
end

% adjust long and lat if needs be
dims = size(long)
if sum(size(long) == size(n))~=2
    [long,lat] = meshgrid(long,lat);
end

% compute CT if needs be
% if exist('CT','var')==0 | length(CT)==0 | sum(size(SA)==size(t))~=sum(size(CT)==size(t))
%     dj_disp('computing CT ...')
%     CT = nan(size(SA));
%     data = SA.*t.*p;
%     Idata = find(isfinite(data)); 
% %     ss = s(inds);
% %     tt = t(inds); 
% %     pp = p(inds);
%     CT(Idata) = gsw_CT_from_t(SA(Idata),t(Idata),p(Idata));
% %     ct = nans(size(t));
% %     ct(inds) = ctt;
% end

whos('global')
dj_toc


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
global SA CT t p g ocean n long lat

[nz,ny,nx] = size(SA)

% eval(['ix = ', handles.grdx, ',']);
% eval(['jy = ', handles.grdy, ',']);
keyboard
ix = handles.grdx;
eval(['jy = ', handles.grdy, ',']);

SA = SA(:,jy,ix);
CT = CT(:,jy,ix);
%t = t(:,jy,ix);
p = p(:,jy,ix); 

if exist('g','var'), g = g(:,jy,ix); end
if exist('ocean','var'), ocean = ocean(jy,ix); end
if exist('n','var'), n = n(jy,ix); end
if exist('long','var'), long = long(jy,ix); end
if exist('lat','var'), lat = lat(jy,ix); end

whos


