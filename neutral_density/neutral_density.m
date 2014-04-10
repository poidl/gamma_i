function varargout = neutral_density(varargin)
% NEUTRAL_DENSITY M-file for neutral_density.fig
%      NEUTRAL_DENSITY, by itself, creates a new NEUTRAL_DENSITY or raises the existing
%      singleton*.
%
%      H = NEUTRAL_DENSITY returns the handle to a new NEUTRAL_DENSITY or the handle to
%      the existing singleton*.
%
%      NEUTRAL_DENSITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEUTRAL_DENSITY.M with the given input arguments.
%
%      NEUTRAL_DENSITY('Property','Value',...) creates a new NEUTRAL_DENSITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before neutral_density_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to neutral_density_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools help.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help neutral_density

% Last Modified by GUIDE v2.5 26-Mar-2008 11:23:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @neutral_density_OpeningFcn, ...
                   'gui_OutputFcn',  @neutral_density_OutputFcn, ...
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


% --- Executes just before neutral_density is made visible.
function neutral_density_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to neutral_density (see VARARGIN)

% Choose default command line output for neutral_density

handles.output = hObject;

handles.ntp = 1;
handles.ns = 1;
handles.gfunc = 0;
handles.sigp = 0;
handles.pr0 = 2000;
handles.glevels = '[  [20 : 1: 26]  [27 : 0.2 : 27.8]  [28 : 0.1 : 28.5]  ]';
handles.eos = 'eos05_ct';
handles.helicity = 1;
handles.epsilon = 1;

%pwd, load handles_last, plast, dj_pause(0)
load handles_last
handles.open_path = plast;
clear plast

% set(0,'Units','pixels')
% set(gcf,'Units','pixels','Position',[390,75,62,30])

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = neutral_density_OutputFcn(hObject, eventdata, handles) 
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


% Choose default command line output for neutral_density
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


% --- Executes on button press in gamma_n_pushbutton.
function gamma_n_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_n_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s ct t p g longs lats 

longss = nanmean(longs(:,:)); latss = nanmean(lats(:,:)');

%figure,dj_pltmp(longss,latss,squeeze(g(21,:,:))), dj_pause(1)

dj_disp('labeling missing gamma_n values with gamma_rf ...')

method = 2; % compile_gl_t

% s = change(s,'==',nan,-99); t = change(t,'==',nan,-99);

% if length(size(s))==2
%         [g,dgl,dgu] = quick_section_glabel(s,t,p,longs,lats);
% elseif length(size(s))==3
%     if length(size(p))==2 
%         [nz,ny,nx] = size(s);
%         p = reshape(p(:)*ones(1,nx*ny),nz,ny,nx);
%     end
% %    longs = longs(:); lats = lats(:); [longs2,lats2] = meshgrid(longs,lats); 
%     dj_tic
%     if method == 1
%         [g,gl,gh] = quick_global_glabel(s,t,p,longs,lats);
%     elseif method == 2
%         inds = find(isfinite(squeeze(s(1,:,:)))&-80<=lats&lats<=64);
%         ss = s(:,inds); tt = t(:,inds); pp = p(:,inds);
%         longss = longs(inds); latss = lats(inds);                          
% %    [gg,dggl,dggh] = gamma_n(ss,tt,pp,longss,latss);
% 
% %         [gg,dggl,dggh] = quick_section_glabel(ss,tt,pp,longss,latss);
% 
%         dj_tic
%         disp('called gamma_n ...'), clear g, whos
%         [gg,dggl,dggh] = gamma_n(ss,tt,pp,longss,latss);
%         dj_toc
%         g = nan*ones(size(s)); g(:,inds) = gg;
%         dj_pause(0), whos
%     else
%         error('****    impossible in function gamma_n_pushbutton    ****')
%     end
%    dj_toc
%     s = change(s,'==',-99,nan); t = change(t,'==',-99,nan);
%     g = change(g,'<',0,nan); g = reshape(g,size(s));
%     
% end

% 
% g = change(g,'==',0,nan);

inds = find(isnan(g)&isfinite(s)); length(inds)

ctt = ct_from_t(s(inds),t(inds),p(inds));figure

cd gfunction

g(inds) = gfunc(s(inds),ctt);

cd ..

figure,dj_pltmp(longss,latss,squeeze(g(21,:,:))), dj_pause(1)


if is_a_section(s,t) == 1
    g(:,:,1) = g(:,:,2); g(:,:,3) = g(:,:,2);
    ok = 'section fixed'
end

ok = 'done', dj_toc




% --- Executes on button press in gamma_rf_pushbutton.
function gamma_rf_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_rf_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc, dj_tic

global s ct p g longs lats ocean n

dj_disp('labeling with gamma_rf ...')

gamma_rf
    
dj_disp('                      ... done')


% --- Executes during object creation, after setting all properties.
function gamma_uipanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_uipanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function gamma_n_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_n_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function gamma_rf_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_rf_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function settings_uipanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to settings_uipanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function neutral_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to neutral_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

clc

global h_glevels

z = matlabroot; cmd = ['cd ', z(1), 'home/bar747/matlab/gamma_i/neutral_density']; eval(cmd) 

handles.ntp = 0;
handles.ns = 1;
handles.gfunc = 0;
handles.sigp = 0;

handles.glevels = '[  [20 : 1: 26]  [27 : 0.2 : 27.8]  [28 : 0.1 : 28.5]  ]'

h_glevels = handles.glevels;

load handles_last
handles.open_path = plast;
clear plast

%my_computer = set_computer

% Update handles structure
guidata(hObject, handles);

% if strcmp(my_computer,'laptop')
%    set(gcf,'Position',[172,53,52,30])
% elseif strcmp(my_computer,'desktop')
% set(0,'Units','pixels')
% set(gcf,'Units','pixels','Position',[390,75,62,30])
% end


% --- Executes on button press in igamma_nversions_pushbutton.
function inversions_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to igamma_nversions_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s t p g longs lats

gamma_n_inversions


% --- Executes on button press in stability_pushbutton.
function stability_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to stability_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s t p g longs lats

buoyancy_frequencies


% --- Executes on button press in options_pushbutton.
function options_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to options_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s t p g longs lats

clc

z = matlabroot; cmd = ['cd ', z(1), 'home/bar747/matlab/gamma_i/neutral_density']; eval(cmd) 

output = tangent_plane_gui;

handles.ntp = output(1);

handles.ns = output(2);

handles.sigp = output(3);

handles.pr0 = output(4)

% Update handles structure
guidata(hObject, handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

global h_glevels

handles.glevels = get(hObject,'String')

h_glevels = handles.glevels

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


% --- Executes on button press in ok1_pushbutton.
function ok1_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ok1_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s ct p g longs lats glevels

global sns tns ctns pns eos

clc

dj_disp('computing density surfaces ...')


cmd = ['glevels = ', handles.glevels, ';']; eval(cmd)

inds = find(isfinite(g(1,:)));

ss = s(:,inds); ctt = ct(:,inds); pp = p(:,inds); gg = g(:,inds);

nlevels = length(glevels); [nz,ny,nx] = size(s);

sns = nan*ones(nlevels,ny,nx); ctns = sns; pns = sns;

% if strcmp(eos,'eos80_t')    
%     [snss,tnss,pnss,dsnss,dtnss,dpnss] = nsfcs(ss,tt,pp,gg,glevels);
%     tns = sns;
%     snss = change(snss,'<=',-99,nan); tnss = change(tnss,'<=',-99,nan); pnss = change(pnss,'<=',-99,nan); 
%     sns(:,inds) = snss; tns(:,inds) = tnss; pns(:,inds) = pnss;
% elseif strcmp(eos,'eos05_ct')  
%     ctt = ct_from_t(ss,tt,pp);
    
[snss,ctnss,pnss,dsnss,dctnss,dpnss] = neutral_surfaces(ss,ctt,pp,gg,glevels);

sns(:,inds) = snss; ctns(:,inds) = ctnss; pns(:,inds) = pnss;

dj_disp('                           ... done')

%     tnss = t_from_ct(snss,ctnss,pnss);
%     tns = nan*ones(size(ctns)); tns(:,inds) = tnss; 
    
%    inds_trpls = find(finite(pnss) & pnss~=0); no_triples = length(inds_trpls)
%    snss(inds_trpls) = nan; tnss(inds_trpls) = nan; pnss(inds_trpls) = nan;
% end

% indss = find(finite(pnss)); length(indss)
% 
% ctmean = nanmean(pnss(indss)); ctmedian = median(pnss(indss));
% 
% stats = [ctmean, ctmedian]
% 
% figure, dj_pltmp(longs(1,:),lats(:,1),squeeze(ctns(1,:,:)))
% 
% title(['mean ', num2str(ctmean,'%8.2f'), '     median ', num2str(ctmedian,'%8.2f')])


% --- Executes on button press in Veronis_pushbutton.
function Veronis_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Veronis_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global s t p g longs lats

global epsilon_x epsilon_y

global called_ggrads

clc

z = matlabroot; cmd = ['cd ', z(1), 'home/bar747/matlab/gamma_i/neutral_density/']; eval(cmd) 

%ok = 'tangent plane gradients ...'

called_ggrads = 0;

[ave,percent2,epsilon_x,epsilon_y] = get_gradients(handles);



% --- Executes on button press in eos05_ct_radiobutton.
function eos05_ct_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to eos05_ct_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eos05_ct_radiobutton

handles.eos = 'eos05_ct';

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in eos80_t_radiobutton.
function eos80_t_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to eos80_t_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eos80_t_radiobutton

handles.eos = 'eos80_t';

% Update handles structure
guidata(hObject, handles);




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

z = matlabroot; cmd = ['cd ', z(1), 'home/bar747/matlab/gamma_i/neutral_density/']; eval(cmd) 

if length(size(p))<3
    [nz,ny,nx] = size(s);
    p = reshape(p(:)*ones(1,nx*ny),nz,ny,nx); 
end

ok = 'loaded data'

if exist('ocean')~=1 | length(ocean)==0
    disp( 'computing oceans and depth ...')
    cd ../data, oceans_and_n, cd ../neutral_density, dj_pause(1)
end

plast = handles.open_path; pwd, save handles_last plast 

whos('global'),  dj_toc




% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s ct t p  g ocean n longs lats

%global s1 t1 

%s = s1; t = t1;

cmd = ['cd ', handles.open_path]; eval(cmd)


[file,path] = uiputfile('*.mat','save data as');

cmd = ['save ', path, file, ' s ct t p g ocean n longs lats '], eval(cmd)

plast = handles.open_path;

z = matlabroot; cmd = ['cd ', z(1), 'home/bar747/matlab/gamma_i/neutral_density']; eval(cmd) 

save handles_last plast 




% --- Executes on button press in gamma_na_pushbutton.
function gamma_na_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_na_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cd ../webstuff
    open neutral_density_css.html
cd ../neutral_density



% --------------------------------------------------------------------
function Save_tmp_Callback(hObject, eventdata, handles)
% hObject    handle to Save_tmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

saveas(gcf,'ndg','bmp')



% --- Executes during object creation, after setting all properties.
function gamma_pr_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_rf_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in gamma_pr_pushbutton.
function gamma_pr_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_pr_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global s ct t p g longs lats ocean n
    
indss = find(finite(s));

ss = s(indss); tt = t(indss); pp = p(indss); 
    
if size(ct)==0
    ctt = ct_from_t(ss,tt,pp);
    ct = nan*ones(size(t)); ct(indss) = ctt;
else
    ctt = ct(indss);
end

cd ../an_equation
    copyfile('rfunc_ct7_9.dat', 'rfunc_ct.dat')
    load rfunc_ct.dat
    g(indss) = rfunc_ct7_9(ss,ctt,rfunc_ct);
    
cd ../neutral_density    

ok = '... labelled'



% --- Executes during object creation, after setting all properties.
function gamma_density_uipanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_density_uipanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in tangent_planes_potential_density_checkbox.
function tangent_planes_potential_density_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to tangent_planes_potential_density_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tangent_planes_potential_density_checkbox

handles.ntp = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in potential_density_checkbox.
function potential_density_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to potential_density_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of potential_density_checkbox

handles.sigp = 1-handles.sigp;

if handles.sigp==1
    handles.ns = 0
else
    handles.ns = 1
end

% Update handles structure
guidata(hObject, handles);


function pr0_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pr0_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pr0_edit as text
%        str2double(get(hObject,'String')) returns contents of pr0_edit as a double

handles.pr0 = str2double(get(hObject,'String'))

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pr0_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pr0_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tangent_planes_checkbox.
function tangent_planes_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to tangent_planes_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tangent_planes_checkbox


% --- Executes on button press in gamma_c_pp_pushbutton.
function gamma_c_pp_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_c_pp_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



gamma_pp_contours



% --- Executes on button press in pressure_pushbutton.
function pressure_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pressure_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s ct p g longs lats ocean n

p_on_g_surfaces


% --- Executes on button press in location_pushbutton.
function location_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to location_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

view_sns_ctns


% --- Executes on button press in gamma_c_r_pushbutton.
function gamma_c_r_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_c_r_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


gamma_r_contours


% --- Executes on button press in gamma_p_pushbutton.
function gamma_p_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_p_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s ct p g longs lats ocean n
    
global h_numerator h_denominator h_normalise

which gamma_p, gamma_p


% --- Executes on button press in extension_pushbutton.
function extension_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to extension_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gamma_surface_extension


% --- Executes on button press in gamma_pp_pushbutton.
function gamma_pp_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_pp_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s t p ct g longs lats ocean n

ok = 1

cd ../an_equation

g = gamma_pp(s,ct,p,g,ocean,longs,lats,15);

cd ../neutral_density


% --- Executes during object creation, after setting all properties.
function gamma_p_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_p_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function gamma_pp_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_pp_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


