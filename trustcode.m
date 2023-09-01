% close all
function varargout = trustcode(varargin)
% TRUSTCODE M-file for trustcode.fig
%      TRUSTCODE, by itself, creates a new TRUSTCODE or raises the existing
%      singleton*.
%
%      H = TRUSTCODE returns the handle to a new TRUSTCODE or the handle to
%      the existing singleton*.
%
%      TRUSTCODE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRUSTCODE.M with the given input arguments.
%
%      TRUSTCODE('Property','Value',...) creates a new TRUSTCODE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trustcode_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trustcode_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trustcode

% Last Modified 23-Aug-2011 12:58:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trustcode_OpeningFcn, ...
                   'gui_OutputFcn',  @trustcode_OutputFcn, ...
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


% --- Executes just before trustcode is made visible.
function trustcode_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trustcode (see VARARGIN)

% Choose default command line output for trustcode
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trustcode wait for user response (see UIRESUME)
% uiwait(handles.Trustcode);


% --- Outputs from this function are returned to the command line.
function varargout = trustcode_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in selectafile.
function selectafile_Callback(hObject, eventdata, handles)
% hObject    handle to selectafile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns selectafile contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectafile
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function selectafile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectafile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[c, b] =uigetfile({'*.IMA';'*.DCM';'*.REC'},'Select .REC or the 1st .IMA file');
d=get(handles.selectafile,'string');
% e=get(handles.FSG,'string')
% C{end+1}=fullfile(b,c);
% d{1:end}=C{end:1};
d=fullfile(b,c);
set(handles.selectafile,'string',d);
% set(handles.FSG,'string',e);
guidata(hObject,handles);


function rownum_Callback(hObject, eventdata, handles)
% hObject    handle to rownum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rownum as text
%        str2double(get(hObject,'String')) returns contents of rownum as a double
input=str2num(get(hObject,'String'));
if(isempty(input))
    set(hObject,'String','0');
else
    input=num2str(input);
set(handles.rownum,'String',input);
end;

guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function rownum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rownum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function columnum_Callback(hObject, eventdata, handles)
% hObject    handle to columnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of columnum as text
%        str2double(get(hObject,'String')) returns contents of columnum as a double
input=str2num(get(hObject,'String'));
if(isempty(input))
    set(hObject,'String','0');
else
    input=num2str(input);
set(handles.columnum,'String',input);
end;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function columnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to columnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dynamic_Callback(hObject, eventdata, handles)
% hObject    handle to dynamic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dynamic as text
%        str2double(get(hObject,'String')) returns contents of dynamic as a double
input=str2num(get(hObject,'String'));
if(isempty(input))
    set(hObject,'String','0');
else
    input=num2str(input);
set(handles.dynamic,'String',input);
end;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function dynamic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dynamic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in eTe.
function eTe_Callback(hObject, eventdata, handles)
% hObject    handle to eTe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eTe
input=get(handles.eTe,'Value');
if(input==1)
    qdiag = questdlg('What is the inner loop?', 'Order','eTE','Rep', 'eTE');
    eTEdiag = questdlg('What is minimal eTE?', 'Min eTE (ms)', '1', '0.44', '0.44');
    min_eTE = str2double(eTEdiag);
    if strcmpi(qdiag, 'eTE')
        seq=[min_eTE min_eTE 40 40 80 80 160 160 min_eTE min_eTE 40 40 80 80 160 160 min_eTE min_eTE 40 40 80 80 160 160];
    elseif strcmpi(qdiag, 'Rep')
        seq=[min_eTE min_eTE min_eTE min_eTE min_eTE min_eTE 40 40 40 40 40 40 80 80 80 80 80 80 160 160 160 160 160 160];
    end
    seq1=sprintf('%g  ', seq);
    set(handles.sequence,'Max',size(seq,1));
    set(handles.sequence,'String',seq1);
else
    set(handles.sequence,'String',' ');
end;
guidata(hObject,handles);



function sequence_Callback(hObject, eventdata, handles)
% hObject    handle to sequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sequence as text
%        str2double(get(hObject,'String')) returns contents of sequence as a double
seq=(get(handles.sequence,'String'));
set(handles.sequence,'String',seq);
 seq=str2num(get(handles.sequence,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Hct_Callback(hObject, eventdata, handles)
% hObject    handle to Hct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Hct as text
%        str2double(get(hObject,'String')) returns contents of Hct as a double
udef=(get(handles.Hct,'String'));
set(handles.Hct,'String',udef);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Hct_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hct.
function hct_Callback(hObject, eventdata, handles)
% hObject    handle to hct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hct
input=get(handles.hct,'Value');
if(input==1)
    % set default hct based on the gender of the subject
    qdiag = questdlg('What is the gender?','gender','Male','Female','Unknown','Unknown');
    if strcmpi(qdiag, 'Male')
        def=0.42;
    elseif strcmpi(qdiag, 'Female')
        def=0.40;
    elseif strcmpi(qdiag, 'Unknown')
        def=0.41;
    end
    def=num2str(def);
    set(handles.Hct,'String',def);
else
    set(handles.Hct,'String',' ');
end;
guidata(hObject,handles);


% --- Executes on button press in deltar2.
function process_Callback(hObject, eventdata, handles)
% hObject    handle to deltar2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath([pwd '\utility']);
if get(handles.checkbox5,'Value')==1
    bIsNeonate = 1;
else
    bIsNeonate = 0;
end

Threshold_epsilon = str2num(get(handles.edit18,'String'));
i=str2num(get(handles.dynamic,'string'));
j=length(str2num(get(handles.sequence,'String')));
if(i~=j||mod(i,2)~=0);
h=errordlg('The length of the sequence and the dynamic number doesnot match or odd dynamic number entered');
set(handles.T2, 'String',' ');
set(handles.R2,'String',' ');
set(handles.stderrR2,'String',' ');
set(handles.R2,'String',' ');
set(handles.DeltaR2,'String',' ');
set(handles.Yv,'String',' ');
set(handles.StderrofYv,'String',' ');
return;
end;
filename1=get(handles.selectafile,'String');
matrix(1)=str2num(get(handles.rownum,'String'));
matrix(2)=str2num(get(handles.columnum,'String'));
matrix(3)=str2num(get(handles.dynamic,'String'));
seq=str2num(get(handles.sequence,'String'));
hct=str2num(get(handles.Hct,'String'));
 set(handles.Trustcode, 'HandleVisibility', 'off');
get(0,'currentfig');

% [t2,ci,Yv]=trust51combine_SS_4te_6dyn_GUI_linux(filename1,matrix,seq,hct);
% Yifan Gou, 8/23/2023: Add ARTS to TRUST process
[t2,ci,Yv, LabelEff] = trustcode_pl_GUI_ARTS(filename1,matrix,seq,hct,bIsNeonate,Threshold_epsilon);
% figure(findall(0,'Tag','MyGUI'));
t2
ci
Yv
LabelEff
t2s=num2str(t2);
set(handles.T2, 'String', t2s);
% ishandle(handles.R2)
r2=1000/(t2);
r2s=num2str(r2);
set(handles.R2,'String',r2s);% S refers to string value of variable
deltar2=1000/ci(2)-1000/ci(1)
deltaR2s=num2str(deltar2);
set(handles.DeltaR2,'String',deltaR2s);
stderr=deltar2/(2*1.96);
stderrs=num2str(stderr);
set(handles.stderrR2,'String',stderrs);
Yvs=num2str(Yv);
set(handles.Yv,'String',Yvs);
r2up=r2+stderr;
r2lw=r2-stderr;
t2up=1000/r2lw;
t2lw=1000/r2up;
Yv1(1) = 100*neoT2toY(t2up,hct);
Yv1(2) = 100*neoT2toY(t2lw,hct);
stderrYv=(Yv1(1)-Yv1(2))/2/1.96;
stderrYvs=num2str(stderrYv);
set(handles.StderrofYv,'String',stderrYvs);
output = [t2, Yv, deltar2, stderrYv, hct]; % just make it easier to copy to excel, not necessary
[fname_path fname_body fname_ext] = fileparts(filename1);
filename = [fname_path filesep fname_body];
fresult=[fname_path filesep 'trust_',fname_body,'.txt'];
    f=fopen(fresult, 'wt');
    fprintf(f,'%s\r\n',filename);
    fprintf(f,'\r\n');
    fprintf(f,'1). T2 (ms):');
    fprintf(f,'%3.3f\r\n',t2);
    fprintf(f,'\r\n');
    fprintf(f,'2). R2 (1/s):');
    fprintf(f,'%3.3f\r\n',r2);
    fprintf(f,'\r\n');
    fprintf(f,'%s\r\n','3). 95% confidence'); 
    fprintf(f,'\n');
    fprintf(f,'interval of R2 (1/s):');
    fprintf(f,'%3.3f\r\n',deltar2);
    fprintf(f,'\r\n');
    fprintf(f,'4). standard error of R2 (1/s):');
    fprintf(f,'%3.3f\r\n',stderr);
    fprintf(f,'\r\n');
    fprintf(f,'%s','5). Yv (%):'); 
    fprintf(f,'%3.3f\r\n',Yv);
    fprintf(f,'\r\n');
    fprintf(f,'%s','6). standard error of Yv (%):'); 
    fprintf(f,'%3.3f\r\n',stderrYv);
    fprintf(f,'\r\n');
    fprintf(f,'%s','7). Labelling Efficiency (%):'); 
    fprintf(f,'%3.3f%%\r\n',LabelEff*100);
    fprintf(f,'\r\n');
    fclose(f);
guidata(hObject,handles)

function T2_Callback(hObject, eventdata, handles)
% hObject    handle to T2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T2 as text
%        str2double(get(hObject,'String')) returns contents of T2 as a double
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function T2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R2_Callback(hObject, eventdata, handles)
% hObject    handle to R2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R2 as text
%        str2double(get(hObject,'String')) returns contents of R2 as a double
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function R2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DeltaR2_Callback(hObject, eventdata, handles)
% hObject    handle to DeltaR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DeltaR2 as text
%        str2double(get(hObject,'String')) returns contents of DeltaR2 as a double
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function DeltaR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeltaR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stderrR2_Callback(hObject, eventdata, handles)
% hObject    handle to stderrR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stderrR2 as text
%        str2double(get(hObject,'String')) returns contents of stderrR2 as a double
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function stderrR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stderrR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Yv_Callback(hObject, eventdata, handles)
% hObject    handle to Yv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Yv as text
%        str2double(get(hObject,'String')) returns contents of Yv as a double
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Yv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Yv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StderrofYv_Callback(hObject, eventdata, handles)
% hObject    handle to StderrofYv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StderrofYv as text
%        str2double(get(hObject,'String')) returns contents of StderrofYv as a double


% --- Executes during object creation, after setting all properties.
function StderrofYv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StderrofYv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Defaultimageparameters.
function Defaultimageparameters_Callback(hObject, eventdata, handles)
% hObject    handle to Defaultimageparameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Defaultimageparameters
input=get(handles.Defaultimageparameters,'Value');
if(input==1)
    matrix(1)=64;
    matrix(2)=64;
    matrix(3)=24;
    m1=num2str(matrix(1));
    m2=num2str(matrix(2));
    m3=num2str(matrix(3));
    set(handles.rownum,'String',m1);
    set(handles.columnum,'String',m2);
    set(handles.dynamic,'String',m3);
else
set(handles.rownum,'String',' ');
set(handles.columnum,'String',' ');
set(handles.dynamic,'String',' ');
end;
    
guidata(hObject,handles);


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox5,'Value')==0
    set(handles.checkbox5,'Value',1);
end
set(handles.checkbox6,'Value',0);
% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox6,'Value')==0
    set(handles.checkbox6,'Value',1);
end
set(handles.checkbox5,'Value',0);
% Hint: get(hObject,'Value') returns toggle state of checkbox6
