function graph_adj()
clear all; close all;
clc;

prompt = {'Enter the number of Vertices : '};
dlg_title = 'No. of Vertices';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);
n = str2num(answer{1});

dist{n,n} = [];
global dist_modify;

 for i=1:n
    for j=1:n
        dist{i,j} = '';
        dist{i,i} = '0';       
    end
 end
 
P = figure('Name','Enter the Adjacency Matrix (No Value Entered will be cosidered as 0)','NumberTitle','off','Position',[0 0 1000 600] );
t = uitable('Parent',P,...
            'Position',[0 100 1000 430],'Data', dist,...
            'FontSize',10,'ForegroundColor','k', ...
            'FontName','Comic Sans MS', ...
            'ColumnEditable', true );


b = uicontrol('Parent',P,...
            'Style','Pushbutton',... 
            'Units','points', ...
            'Callback',@Pushbutton1_Callback,...
         'Position',[110 410 83.1724 30.4138], ...
         'String','Click To Solve', ...
         'Tag','checkbox1' );
 P;
 
function Pushbutton1_Callback(hObject,eventdata)
if (get(hObject,'Value') == get(hObject,'Max'))

g = get(t,'Data');
new = {};
for i=1:n
    for j=1:n      
        new{i,j} = str2num(g{i,j});
    end
end

for i=1:n
    for j=1:n
        sz = size(new{i,j});
        if sz(2)==0
           new{i,j} = 0; 
        end    
    end
end    

for i=1:n
    for j=1:n
        temp = new{i,j};         
            if new{i,j} ~= new{j,i}
               str = strcat('Edge ',num2str(i),',',num2str(j));
               waitfor(msgbox('Values mismatch for same pair of vertices in Adjacency Matrix ! Re-Enter the Values.',str,'warn')); 

               prompt = {'Enter the new value : '};
               dlg_title = str;
               num_lines = 1;
               def = {'0'};
               answer = inputdlg(prompt,dlg_title,num_lines,def);

               temp = answer{1};
               new{i,j} = str2num(temp);
               new{j,i} = str2num(temp);
            else
                new{i,j} = temp;
                new{j,i} = temp;
            end    
    end
end

dist = new;
dist_modify = {};
for i=1:n
    for j=1:n
        dist_modify{i,j} = num2str(dist{i,j});
    end
end    

plot_graph(dist,n)

end  
end

%% Final Plot Function
function plot_graph(dist,n)
close all;    
theta = 2*pi/n;

vertices = [];
r = 450;
c = 1;

for x=0:theta:(2*pi)
    xcord = r*cos(x);
    w1 = xcord/power(10,-10);
    if abs(w1)<1
         xcord=0;  
    end   
    ycord = r*sin(x);
    w2 = ycord/power(10,-10);
    if abs(w2)<1
         ycord=0;  
    end 
    vertices(c) = (xcord) + 1i*(ycord);
   c = c+1;
end     

for i=1:n
    X(i) = real(vertices(i));
    Y(i) = imag(vertices(i));
end

T = figure('Name','Graph','NumberTitle','off','Position',[50 50 1200 600],...
            'DeleteFcn',@Figure_Close_Callback);
plot(real(vertices),imag(vertices),'ob','LineWidth',4,...
                'MarkerSize',10)        

for i=1:n
    if X(i)>=0 && Y(i)>=0
        text(X(i)+20,Y(i)+20,num2str(i),'FontSize',18,'Color','r'); 
    elseif X(i)<=0 && Y(i)>=0
        text(X(i)-30,Y(i)+30,num2str(i),'FontSize',18,'Color','r');
    elseif X(i)<=0 && Y(i)<=0
        text(X(i)-30,Y(i)-30,num2str(i),'FontSize',18,'Color','r');
    elseif X(i)>=0 && Y(i)<=0    
        text(X(i)+20,Y(i)-20,num2str(i),'FontSize',18,'Color','r');
    end    
end     

for i=1:n
    for j=1:n
          g = [X(i) X(j)];
          h = [Y(i) Y(j)];
          if(dist{i,j} == 1)
             line(g,h,'LineWidth',1.5,'Color','k')
          end   
    end
end

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
                    graph_adj_modify(dist_modify,n);
            end       
end
end