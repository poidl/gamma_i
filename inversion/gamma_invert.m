function varargout = gamma_invert(varargin)
% GAMMA_INVERT M-file for gamma_invert.fig
%      GAMMA_INVERT, by itself, creates a new GAMMA_INVERT or raises the existing
%      singleton*.
%
%      H = GAMMA_INVERT returns the handle to a new GAMMA_INVERT or the handle to
%      the existing singleton*.
%
%      GAMMA_INVERT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAMMA_INVERT.M with the given input arguments.
%
%      GAMMA_INVERT('Property','Value',...) creates a new GAMMA_INVERT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gamma_invert_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gamma_invert_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gamma_invert

% Last Modified by GUIDE v2.5 14-Sep-2007 12:04:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gamma_invert_OpeningFcn, ...
                   'gui_OutputFcn',  @gamma_invert_OutputFcn, ...
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


% --- Executes just before gamma_invert is made visible.
function gamma_invert_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gamma_invert (see VARARGIN)

% Choose default command line output for gamma_invert

clc

gamma_invert_dir = which('gamma_invert.m');
cd ([gamma_invert_dir(1:end-15)])
% z = matlabroot; 
% cmd = ['cd ', z(1), ':/neutrals/ness8/inversion']; 
% eval(cmd) 

set(0,'DefaultFigureColor',[0.961 0.988 0.965])

handles.output = hObject;
handles.inv_meth = 1;
handles.eos = 'eos80_t'; 
handles.zbc = 0;
handles.quad = 0;
handles.linb = 0;
handles.h2vwt = 1;
handles.helwt = 1e-6;
handles.bdrywt = 1e-5;

load handles_last 
 handles.open_path = plast; 
 clear plast

handles.maxits = 100;

% Update handles structure
guidata(hObject, handles);

% my_computer = set_computer
% 
% if strcmp(my_computer,'laptop')
%    set(gcf,'Position',[1,1,626,438])
% elseif strcmp(my_computer,'desktop')
%    set(gcf,'Position',[190,43,53,26])
% end


% UIWAIT makes gamma_invert wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gamma_invert_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in init_pushbutton.
function init_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to init_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s t ct p z g longs lats ocean n 

global k_east r_east k_west r_west k_north r_north k_south r_south

clc, dj_tic

%compute_east = 0; compute_west = 0; compute_north = 0; compute_south = 0;

compute_east = 1; compute_west = 1; compute_north = 1; compute_south = 1;
 
compute_helicities = 0;

compute_vertical = 0; 

compute_vertical_xzb = 0;

compute_gvals = 1; 

bc = handles.zbc;

if strcmp(handles.eos,'eos80_t')    
	eosno = 1      
elseif strcmp(handles.eos,'eos05_th')  
	eosno = 2      
elseif strcmp(handles.eos,'eos05_ct')
	eosno = 3    
end

if size(ct) == 0
    dj_disp('computing conservative temperature ... to fix for eos80_t')
    inds = find(isfinite(s(1,:)));
    ss = s(:,inds);
    tt = t(:,inds);
    pp = p(:,inds);
    ctt = ct_from_t(ss,tt,pp);
    ct = nan*ones(size(s)); 
    ct(:,inds) = ctt;
end

figure

if compute_east==1
    dj_disp('computing easterly intersections ...')
    handles

    if strcmp(handles.eos,'eos80_t')
        [k_east,r_east] = intersections_east(s,t,p,ocean,n,longs,lats,handles);
    elseif strcmp(handles.eos,'eos05_th')
        error('not yet implemented in gamma_invert')
        %[k_east,r_east] = intersections_east(s,ct,p,ocean,n,longs,lats,eos);
    elseif strcmp(handles.eos,'eos05_ct')
        [k_east,r_east] = intersections_east(s,ct,p,ocean,n,longs,lats,eosno,bc);
    end
    save intersections/intersections_east k_east r_east
    dj_toc
else
    load intersections/intersections_east
end
       
if compute_west==1
    dj_disp('computing westerly intersections ...')
    if strcmp(handles.eos,'eos80_t')
        [k_west,r_west] = intersections_west(s,t,p,ocean,n,longs,lats,handles);
    elseif strcmp(handles.eos,'eos05_th')
        error('not yet implemented in gamma_invert')
        %[k_east,r_east] = intersections_east(s,ct,p,ocean,n,longs,lats,eos);
    elseif strcmp(handles.eos,'eos05_ct')
        [k_west,r_west] = intersections_west(s,ct,p,ocean,n,longs,lats,eosno);
    end
    save intersections/intersections_west k_west r_west
	dj_toc
else
    load intersections/intersections_west
end
    
if compute_north==1
    dj_disp('computing northerly intersections ...')
    if strcmp(handles.eos,'eos80_t')
        [k_north,r_north] = intersections_north(s,t,p,ocean,n,longs,lats,handles);
    elseif strcmp(handles.eos,'eos05_th')
        error('not yet implemented in gamma_invert')
        %[k_east,r_east] = intersections_east(s,ct,p,ocean,n,longs,lats,eos);
    elseif strcmp(handles.eos,'eos05_ct')
        [k_north,r_north] = intersections_north(s,ct,p,ocean,n,longs,lats,eosno);
    end
    save intersections/intersections_north k_north r_north
	dj_toc
else
    load intersections/intersections_north
end

if compute_south == 1
    dj_disp('computing southerly intersections ...')
    if strcmp(handles.eos,'eos80_t')
        [k_south,r_south] = intersections_south(s,t,p,ocean,n,longs,lats,handles);
    elseif strcmp(handles.eos,'eos05_th')
        error('not yet implemented in gamma_invert')
        %[k_east,r_east] = intersections_east(s,ct,p,ocean,n,longs,lats,eos);
    elseif strcmp(handles.eos,'eos05_ct')
        [k_south,r_south] = intersections_south(s,ct,p,ocean,n,longs,lats,eosno);
    end
    save intersections/intersections_south k_south r_south
	dj_toc
else
    load intersections/intersections_south
end

if compute_helicities == 4
    dj_disp('computing helicity ...')
    [h,h_east,h_north,h_west,h_south] = helicities_4(s,t,ct,p,ocean,n,longs,lats,bc);
 	save intersections/helicities h h_east h_north h_west h_south
	dj_toc
else
    load intersections/helicities
end

if compute_vertical == 1
    dj_disp('computing vertical vitals ... to fix for eos80_t')
    [k_vert,r_vert,gprod_vert,no_eqs] = vertical_vitals(s,ct,p,g,ocean,n,longs,lats,handles);
 	save intersections/vertical_vitals k_vert r_vert gprod_vert no_eqs
	dj_toc
else
    load intersections/vertical_vitals
end

if compute_vertical_xzb == 1;
    [gc1,gc2,gc3,gc4,gc5,nxz_beqs] = bxz_equations(s,ct,p,g,ocean,n,longs,lats,handles);
    save intersections/bxz_equations gc1 gc2 gc3 gc4 gc5 nxz_beqs
else
    load intersections/bxz_equations
end

if compute_gvals == 1
    dj_disp('computings deciles  ... ')
    [inds_bg, g_bdry] = boundary_gammas(g,long,lat,n);
    gamma_invert_dir = which('gamma_invert.m');
 	save ([[gamma_invert_dir(1:end-15)],'/intersections/gamma_boundary inds_bg g_bdry'])
	dj_toc
else
    load intersections/gamma_boundary
end

% --- Executes on button press in zonalbc_quadratic_radiobutton.
function zonalbc_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to zonalbc_quadratic_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zonalbc_quadratic_radiobutton

handles.zbc = 1-handles.zbc

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in quadratic_radiobutton.
function quadratic_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to quadratic_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of quadratic_radiobutton

handles.quad = 1-handles.quad;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in linb_radiobutton.
function linb_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to linb_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of linb_radiobutton

handles.linb = 1-handles.linb

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in direct_radiobutton.
function direct_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to direct_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of direct_radiobutton

handles.inv_meth = 1

% Update handles structure
guidata(hObject, handles);





% --- Executes on button press in iterative_radiobutton.
function iterative_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to iterative_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iterative_radiobutton


handles.inv_meth = 2

% Update handles structure
guidata(hObject, handles);





function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


handles.h2vwt = str2double(get(hObject,'String'))

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function density_weight_Callback(hObject, eventdata, handles)
% hObject    handle to density_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density_weight as text
%        str2double(get(hObject,'String')) returns contents of density_weight as a double


handles.helwt = str2double(get(hObject,'String'))

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function density_weight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function boundary_weight_Callback(hObject, eventdata, handles)
% hObject    handle to boundary_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of boundary_weight as text
%        str2double(get(hObject,'String')) returns contents of boundary_weight as a double


handles.bdrywt = str2double(get(hObject,'String'))

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function boundary_weight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to boundary_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear global, clc, dj_tic

global s t ct p z g longs lats ocean n

%z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/inversion']; eval(cmd)

try
    handles.open_path
    cmd = ['cd ', handles.open_path]; eval(cmd)
catch
    gamma_invert_dir = which('gamma_invert.m');
    handles.open_path = gamma_invert_dir(1:end-15);
end
[FileName,PathName] = uigetfile('*.mat');

cmd = ['load ', PathName, FileName]; eval(cmd)

handles.open_path = PathName

% Update handles structure
guidata(hObject, handles);

%z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/inversion']; eval(cmd) 
gamma_invert_dir = which('gamma_invert.m');
cd ([gamma_invert_dir(1:end-15)])

if length(size(p))<3
    [nz,ny,nx] = size(s);
    p = reshape(p(:)*ones(1,nx*ny),nz,ny,nx); 
end

% unknown
% indss = find(~isfinite(s(1,:)));
% s(:,indss) = nan;
% t(:,indss) = nan;
% p(:,indss) = nan;
% g(:,indss) = nan; 

if exist('ocean')~=1 | length(ocean)==0
    disp( 'computing oceans and depth ...')
    cd ../data, oceans_and_n, cd ../inversion, dj_pause(1)
end

ok = 'loaded data'
whos('global')
dj_toc

plast = handles.open_path;
save handles_last plast 



% --- Executes on button press in invert_pushbutton.
function invert_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to invert_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s t ct p z g longs lats ocean n

global gn limit_panel2

gamma_invert_dir = which('gamma_invert.m');
cd ([gamma_invert_dir(1:end-15)])

% z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/inversion']; eval(cmd) 

gn = g;

invert



% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s t ct p z g longs lats ocean n 

%z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/inversion']; eval(cmd) 
try
    cmd = ['cd ', handles.open_path]; eval(cmd)
catch
    gamma_invert_dir = which('gamma_invert.m');
    cd ([gamma_invert_dir(1:end-15)])
end

[file,path] = uiputfile('*.mat','save data as');
ok = 'here'

cmd = ['save ', path, file, ' s t ct p z g longs lats ocean n'], eval(cmd)

gamma_invert_dir = which('gamma_invert.m');
cd ([gamma_invert_dir(1:end-15)])

% z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/inversion']; eval(cmd) 

plast = handles.open_path
save handles_last plast 





% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

open local_help.html





function iterations_Callback(hObject, eventdata, handles)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterations as text
%        str2double(get(hObject,'String')) returns contents of iterations as a double


handles.maxits = str2double(get(hObject,'String'))

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function iterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


