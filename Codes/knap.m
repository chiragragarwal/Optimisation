function knap()
close all
clear all
clc

global n;
global W;
global g_modify;

prompt = {'Enter the total number of items : '};
dlg_title = 'Enter Data';
num_lines = 1;
def = {'5'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
n = str2num(answer{1});

% Input Benefits and Weights for each Item
dat{n,2} = [];

 for i=1:n
    for j=1:2
        dat{i,j} = '';
    end
 end

columnname = {};
rowname = {};

columnname{1} = 'Weights';
columnformat{1} = 'numeric';


columnname{2} = 'Benefits';
columnformat{2} = 'numeric';


for i=1:n
   rowname{i} = strcat('Item ',num2str(i)); 
end

P = figure('Name','Fractional Knapsack Algorithm','NumberTitle','off','Position',[50 50 400 400] );

t = uitable('Parent',P,...
            'Position',[0 0 400 320],'Data', dat,... 
            'ColumnName', columnname,...
            'ColumnEditable', true,...
            'ColumnFormat', columnformat,...
            'RowName',rowname, ...
            'FontSize',12,'ForegroundColor','k', ...
            'FontName','Comic Sans MS');

% ADDING A BUTTON AT THE BOTTOM

b1 = uicontrol('Parent',P,...
            'Style','Pushbutton',... 
            'Units','points', ...
            'Callback',@Pushbutton1_Callback,...
         'Position',[10 250 83.1724 30.4138], ...
         'String','Click To Solve', ...
         'Tag','checkbox1' ); 
     
b3 = uicontrol('Parent',P,...
            'Style','edit',... 
            'String','','Position',[270 335 120 30],...
            'FontSize',12,'ForegroundColor','k', ...
            'FontName','Comic Sans MS');
        
b4 = uicontrol('Parent',P,...
            'Style','text',... 
            'FontSize',13,'ForegroundColor','k', ...
            'FontName','Times New Roman', ...
             'Position',[150 335 120 30], ...
            'String','Max. Weight = '); 
        
        
 P;
 function Pushbutton1_Callback(hObject,eventdata)
if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take appropriate action
    g = get(t,'Data');
    g_modify = g;
    
    W = get(b3,'String');
    W = str2num(W);
    
    t = zeros(n,2);
    for i=1:n
        for j=1:2
            t(i,j) = str2num(g{i,j});
        end
    end 
    
    
    close all;
    knap_solve(t,n,W)     % Fucntion Call         
end
end  % end of Pushbutton Callback    

function knap_solve(t,n,W)
c = 0;
b = zeros(1,n) ;
w = b ;
x = w ;
v = x ;
index = 0;
wsum = 0;
for i=1:n
   w(i) = t(i,1);
   b(i) = t(i,2);
   wsum = wsum + w(i);
end  

% Value Array
for j = 1:n
	v(j) = b(j)/w(j);
end

k = 1;
while c<min(W,wsum)
	v_max = max(v);
    v_max_index = find(v == max(v));
    x(k) = min(w(v_max_index),W-c);
    c = c + x(k);
    index(k) = v_max_index;
    k = k+1;
    v(v_max_index) = -Inf;   
end   

dat1 = {};
rowname1 = {};
for i=1:k-1
       dat1{i,1} = num2str(index(i));
       dat1{i,2} = num2str(x(i)); 
       rowname1{i} = num2str(i);
end    

columnname1 = {};


columnname1{1} = 'Item Number';
columnformat1{1} = 'numeric';


columnname1{2} = 'Weight';
columnformat1{2} = 'numeric';

% Final Printing
T = figure('Name','Fractional Knapsack Algorithm','NumberTitle','off','Position',[50 50 350 400],...
           'DeleteFcn',@Figure_Close_Callback);

t1 = uitable('Parent',T,...
            'Position',[0 0 350 320],'Data', dat1,... 
            'ColumnName', columnname1,...
            'ColumnEditable', false,...
            'RowName',rowname1, ...
            'ColumnFormat',columnformat1,...
            'FontSize',12,'ForegroundColor','k', ...
            'FontName','Comic Sans MS');

% ADDING A BUTTON AT THE BOTTOM

b2 = uicontrol('Parent',T,...
            'Style','text',... 
            'FontSize',13,'ForegroundColor','k', ...
            'FontName','Times New Roman', ...
             'Position',[20 340 200 35], ...
            'String','FINAL SELECTION'); 
        
 T;

end

function Figure_Close_Callback(hObject,eventdata)
        %          Construct a questdlg with two options
            choice = questdlg('Would you like to solve a new problem ?', ...
             'Solution Achieved. Try another ?', ...
             'Yes','No','View/Modify Input Data','No');
            % Handle response
            switch choice
                case 'Yes'
                    close all; clear all;
                    start();
                    
                case 'No'
                    return  
                    
                case 'View/Modify Input Data'   
                   knap_modify(g_modify,n,W);
            end       
end

end
