function graph_complete()
clear all;
close all;
clc;

prompt = {'Enter the number of Vertices : '};
dlg_title = 'No. of Vertices';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);
n = str2num(answer{1});

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

T = figure('Name','Complete Graph','NumberTitle','off','Position',[50 50 1200 600],...
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
          line(g,h,'LineWidth',1.5,'Color','k')
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
                    graph_complete();
            end       
end

end

