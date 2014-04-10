function varargout = fitting_gui(varargin)
% FITTING_GUI M-file for fitting_gui.fig
%      FITTING_GUI, by itself, creates a new FITTING_GUI or raises the existing
%      singleton*.
%
%      H = FITTING_GUI returns the handle to a new FITTING_GUI or the handle to
%      the existing singleton*.
%
%      FITTING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITTING_GUI.M with the given input arguments.
%
%      FITTING_GUI('Property','Value',...) creates a new FITTING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fitting_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fitting_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help fitting_gui

% Last Modified by GUIDE v2.5 19-Oct-2007 15:57:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fitting_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @fitting_gui_OutputFcn, ...
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


% --- Executes just before fitting_gui is made visible.
function fitting_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fitting_gui (see VARARGIN)

% Choose default command line output for fitting_gui

dj_tic, clc

global h_numerator h_denominator h_normalise

cd d:/neutrals/ness8/an_equation

handles.output = hObject;

handles.ct = 1;

handles.numerator = '15';   h_numerator = handles.numerator; 
handles.denominator = '0';  h_denominator = handles.denominator;

handles.lin = 0;
handles.nonlin = 1;

handles.directory = 'data_new/R';

handles.normalise = 0;      h_normalise = handles.normalise;

handles.boundary = 1;

load handles_last
handles.open_path = plast
clear plast
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fitting_gui wait for user response (see UIRESUME)
% uiwait(handles.fitting_gui);


% --- Outputs from this function are returned to the command line.
function varargout = fitting_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_ct.
function popupmenu_ct_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_ct contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ct


% Choose default command line output for fitting_gui
%handles.output = hObject;

%handles.popup_ct = get(hObject,'String'); % returns popupmenu_ct contents as cell array
%handles.ct = handles.popup_ct{get(hObject,'Value')};

handles.ct = get(hObject,'Value')


% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_ct_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_ct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in view_pushbutton.
function view_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to view_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s ct t p g ocean n longs lats 

view_code


% --- Executes on selection change in numerator_popupmenu.
function numerator_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to numerator_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns numerator_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from numerator_popupmenu

%contents = get(hObject,'String');
%handles.npars = contents{get(hObject,'Value')}

contents = get(hObject,'String');

handles.numerator = contents{get(hObject,'Value')}

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function numerator_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numerator_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in contours_pushbutton.
function contours_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to contours_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contour_code



% --- Executes during object creation, after setting all properties.
function ok_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contours_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in gamma_model.
function gamma_model_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns gamma_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gamma_model

contents = get(hObject,'String');
handles.gamma_model = contents{get(hObject,'Value')}


% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function gamma_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function popupmenu_area_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in ok_pushbutton.
function ok_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ok_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function contours_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contours_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on button press in save_pushbutton.
function save_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to save_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cmd = ['copyfile(''gamma_p', handles.numerator, '.dat'', ''', handles.directory, ''')']

eval(cmd)


% --- Executes on button press in save_pushbutton.
function restore_button_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to save_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in linear_linear_radiobutton.
function linear_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to linear_linear_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of linear_linear_radiobutton


handles.lin = get(hObject,'Value');

handles.nonlin = 1-handles.lin

% Update handles structure
guidata(hObject, handles)



% --- Executes on button press in view_pushbutton.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to view_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc

if handles.ct==1
    uh_oh = 'not yet implemented!'
else
    conservative_temperature = 0
end

data_set =  handles.data

limit_calculations = handles.area-1

set_computer

file_formats

if ~exist('s'), fgui_data, end

twigs_code





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

global s t p g ocean n longs lats

cmd = ['cd ', handles.open_path]; eval(cmd)

[FileName,PathName] = uigetfile('*.mat');

cmd = ['load ', PathName, FileName]; eval(cmd)

handles.open_path = PathName

% Update handles structure
guidata(hObject, handles);

z = matlabroot; cmd = ['cd ', z(1), ':/neutrals/ness8/an_equation']; eval(cmd) 

if length(size(p))<3
    [nz,ny,nx] = size(s);
    p = reshape(p(:)*ones(1,nx*ny),nz,ny,nx); 
end

ok = 'loaded data'

if exist('ocean')~=1 | length(ocean)==0
    disp( 'computing oceans and depth ...')
    cd ../data, oceans_and_n, cd ../an_equation, dj_pause(1)
end

plast = handles.open_path; pwd, save handles_last plast 

whos('global'),  dj_toc


% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in denominator_popupmenu.
function denominator_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to denominator_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns denominator_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from denominator_popupmenu


contents = get(hObject,'String');

handles.denominator = contents{get(hObject,'Value')}

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in restore_button.
function restore_button_Callback(hObject, eventdata, handles)
% hObject    handle to restore_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cmd = ['copyfile(''./', handles.directory, '/gamma_p', handles.numerator, '.dat'', ''',pwd, ''')']

eval(cmd)



function directory_Callback(hObject, eventdata, handles)
% hObject    handle to directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of directory as text
%        str2double(get(hObject,'String')) returns contents of directory as a double


 handles.directory = get(hObject,'String')
 
 
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function directory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in global_rfunc_pushbutton.
function global_rfunc_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to global_rfunc_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cmd = ['copyfile(''rfunc_', handles.numerator, '_', handles.denominator, '_global.dat'',', ... 
                    '''rfunc_', handles.numerator, '_', handles.denominator, '.dat'')']

eval(cmd)


% --- Executes on button press in fit_pushbutton.
function fit_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to fit_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s ct t p g ocean n longs lats 

global h_numerator h_denominator h_normalise

clc

h_numerator = handles.numerator; h_denominator = handles.denominator;

h_normalise = handles.normalise;

cd d:/neutrals/ness8/an_equation

if handles.lin
    fitl_rfunc
elseif handles.nonlin
    nm
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


handles.normalise = get(hObject,'Value')

% Update handles structure
guidata(hObject, handles)




% --- Executes on selection change in popupmenu11.
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11

global h_numerator

contents = get(hObject,'String');

handles.numerator = contents{get(hObject,'Value')}

h_numerator = handles.numerator;

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu12.
function popupmenu12_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu12

global h_denominator

contents = get(hObject,'String');

handles.denominator = contents{get(hObject,'Value')}

% Update handles structure
guidata(hObject, handles);

h_denominator = handles.denominator;


% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in normalise_variables_checkbox.
function normalise_variables_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to normalise_variables_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalise_variables_checkbox

global h_normalise

handles.normalise = get(hObject,'Value')

h_normalise = handles.normalise;

% Update handles structure
guidata(hObject, handles)



% --- Executes on button press in boundary_equations_checkbox.
function boundary_equations_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to boundary_equations_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of boundary_equations_checkbox


handles.boundary = get(hObject,'Value')

% Update handles structure
guidata(hObject, handles)



% --- Executes on button press in linear_radiobutton.
function radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to linear_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of linear_radiobutton




% --- Executes on button press in non_linear_radiobutton.
function non_linear_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to non_linear_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of non_linear_radiobutton



handles.nonlin = get(hObject,'Value'); 

handles.lin = 1-handles.nonlin

% Update handles structure
guidata(hObject, handles);


