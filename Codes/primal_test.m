function primal_test()
clear all;
close all;
clc;
prompt = {'Enter the number of variables : ','Enter the number of constraints :'};
dlg_title = 'Constraints/Variables';
num_lines = 1;
def = {'2','2'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
v = str2num(answer{1});
cons = str2num(answer{2});

global cons_org;
global v_org;
v_org = v;
global g_modify;

global U_indices;
global U_indices_new;
global U_id;
global f1;
f1 = 0;
global l_check;
global l_bound;
l_check = 0;
l_bound = 0;

dat{cons+4,v+2} = [];

daddy_dat = [];         %% Final solution table
dat1 = [];              %% each iteration tables
daddy_row = [];         %% Final rows

for i=1:cons+4
    for j=1:v+2
        dat{i,j} = '';
    end
end

for i=1:v
    dat{cons+2,i} = '0';
    dat{cons+3,i} = 'Inf';
    dat{cons+4,i} = 'n';
end

columnname = {};
columneditable = {};
columnformat = {};

rowname = {};

for i=1:v
    str1 = strcat('x',num2str(i));
    columnname{i} = str1 ;
    columnformat{i} = 'numeric';
    columneditable{i} = true;
end

columnname{i+1} = 'Enter <= , >= or = ';
columnformat{i+1} = {'<=' '>=' '='};
columneditable{i+1} = true;
columneditable{i+2} = true;
columnname{i+2}  = 'R.H.S.';
columnformat{i+2}= 'numeric';

rowname = {};
rowname{1} = 'Objective Function';
global rowname1;
rowname1 = [];

for i=2:cons+1
    rowname{i} = strcat('Constraint ',num2str(i-1));
    
end

rowname{i+1} = 'Lower Bound';
rowname{i+2} = 'Upper Bound';

rowname{i+3} = 'Unrestricted ?(y/n)';

P = figure('Name','All Slack Primal Simplex Table Values','NumberTitle','off','Position',[50 50 800 400] );

t = uitable('Parent',P,...
    'Position',[0 0 800 300],'Data', dat,...
    'ColumnName', columnname,...
    'ColumnFormat', columnformat,...
    'ColumnEditable', true,...    
    'FontSize',10,'ForegroundColor','k', ...
    'FontName','Comic Sans MS', ...
    'RowName',rowname );

% ADDING A BUTTON AT THE BOTTOM

b = uicontrol('Parent',P,...
    'Style','Pushbutton',...
    'Units','points', ...
    'Callback',@Pushbutton1_Callback,...
    'Position', [230 260 83.1724 30.4138], ...
    'String','Click To Solve', ...
    'Tag','checkbox1' );

% Create the button group.
h = uibuttongroup('Parent',P,'visible','on','Position',[0 0 83.1724 30.4138]);
% Create two radio buttons in the button group.
u0 = uicontrol('Style','Radio','String','Maximize',...
    'pos',[100 350 100 30],'parent',h,'HandleVisibility','on');
u1 = uicontrol('Style','Radio','String','Minimize',...
    'pos',[200 350 100 30],'parent',h,'HandleVisibility','on');
set(h,'Visible','on');

P;
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
%                     close all; clear all;
                    return  
                    
                case 'View/Modify Input Data'   
                    primal_modify(g_modify,v_org,cons_org);
            end       
    end

    function Pushbutton1_Callback(hObject,eventdata)
        if (get(hObject,'Value') == get(hObject,'Max'))
            % Checkbox is checked-take appropriate action
            g = get(t,'Data');
            g_modify = g;
            for i=1:cons+1
                for j=1:v
                    sz = size(g{i,j});
                    if sz(2)==0
                       g{i,j} = '0'; 
                    end    
                end
            end
      
             for i=2:cons+1
                    sz = size(g{i,v+2});
                    if sz(2)==0
                       g{i,v+2} = '0'; 
                    end      
             end
            
            for i=2:cons+1
                    sz = size(g{i,v+1});
                    if sz(2)==0
                       g{i,v+1} = '<='; 
                    end      
            end
            
            A = zeros(cons,v);
            for i=2:cons+1
                for j=1:v
                    A(i-1,j) = str2num(g{i,j});
                end
            end
           
            
            b = zeros(1,cons);
            for i=2:cons+1
                b(i-1) = str2num(g{i,v+2});
            end
      
            
            c = zeros(1,v);
            for j=1:v
                c(j) = str2num(g{1,j});
            end
         
            
            inq = zeros(1,cons);
            k=1;
            
            for i=2:(cons+1)
                if strcmp('<=',g{i,v+1}) == 1
                  
                    inq(k) = -1;
                    k = k+1;
                end
                
                if strcmp('=',g{i,v+1}) == 1
                   
                    inq(k) = 0;
                    k = k+1;
                end
                
                if strcmp('>=',g{i,v+1}) == 1
                 
                    inq(k) = 1;
                    k = k+1;
                end
            end
            
            % Check if lower bound is greater than upper bound
            for er=1:v
                temp1 = g{cons+2,er};
                temp2 = g{cons+3,er};
                   [temp1 temp2] = check(temp1,temp2,g,er);
                g{cons+2,er} = temp1;
                g{cons+3,er} = temp2;
            end
            
            cons_org = cons;
            
           % Check for upper bounds
            j1=0;
            for i=1:v
               if str2num(g{cons+3,i}) ~= Inf
                  temp = zeros(1,v);
                  temp(i) = 1;
                  A = vertcat(A,temp);
                  inq = horzcat(inq,-1);
                  b = horzcat(b,str2num(g{cons+3,i}));
                  j1 = j1+1;
               end    
            end    
          
            % Check for lower bounds
            j2=0;
            for i=1:v
               if str2num(g{cons+2,i}) ~= 0
                  b = b';
%                   b = b + A(:,i)*str2num(g{cons+2,i});
                  b = b - A(:,i)*str2num(g{cons+2,i});                  
                  f1 = c(i)*str2num(g{cons+2,i});
                  j2 = j2+1;
                  b = b';
                  l_check(i) = 5;
                  l_bound(i) = str2num(g{cons+2,i});
               else   
                  l_check(i) = 0;
                  l_bound(i) = 0;
               end    
            end    
            
            cons = cons + j1 + j2;
            
            % Check for Unrestricted
            j = 0;
            U_indices = 0;
            U_indices_new = 0;
            U_id = 0;
            
            g_modify = g;
            
            for i=1:v
               if g{cons_org+4,i} == 'y'
                    if str2num(g{cons_org+3,i}) ~= Inf || str2num(g{cons_org+2,i}) ~= 0
                        waitfor(msgbox('A variable cannot be both unrestricted and have bounds simultaneously !','ERROR !','warn'));
                        primal_modify(g_modify,v_org,cons_org)
                    end    
                   Atemp = -A(:,i+j);
                   Apart1 = A(:,(1:i+j));
                   Apart2 = A(:,(i+j+1:v+j));
                   
                   Ctemp = -c(:,i+j);
                   Cpart1 = c(:,(1:i+j));
                   Cpart2 = c(:,(i+j+1:v+j));
                   
                   lctemp = 0;
                   lcpart1 = l_check(:,(1:i+j));
                   lcpart2 = l_check(:,(i+j+1:v+j));
                   
                   lbtemp = 0;
                   lbpart1 = l_bound(:,(1:i+j));
                   lbpart2 = l_bound(:,(i+j+1:v+j));
                   
                   l_check = horzcat(lcpart1,lctemp,lcpart2);
                   l_bound = horzcat(lbpart1,lbtemp,lbpart2);
                   
                   U_indices(i) = 1;
                   U_indices_new(i+j) = 1;
                   U_indices_new(i+j+1) = 2;
                   U_id(i+j) = i;
                   U_id(i+j+1) = i;
                   
                   A = horzcat(Apart1,Atemp,Apart2);
                   c = horzcat(Cpart1,Ctemp,Cpart2);
                                    
                   j = j+1;
               else
                   U_indices(i) = 0;
                   U_indices_new(i+j) = 0;
                   U_id(i+j) = i;
               end    
            end    
            U_indices = horzcat(U_indices,zeros(1,j));
            v = v+j;

            for i=v+1:v+cons
               U_id(i) = i; 
            end    
            
            %% Check for maximize or minimize
            
            var = get(u0,'Value');
            if var == 1
                minimize = 0;
            else
                new = get(u1,'Value');
                if new == 1
                    minimize = 1;
                end
            end
            close all;
            primal_solve(c,b,A,inq,minimize,U_indices,U_indices_new,U_id,f1,l_check,l_bound)     % Fucntion Call
        end
    end  % end of Pushbutton Callback

%% Start of Primal Simplex
    function primal_solve(c,b,A,inq,minimize,U_indices,U_indices_new,U_id,f1,l_check,l_bound)
  
        [m,n] = size(A);
        b = b' ;
        global u;
        u = 1;
       
        % Check for negative RHS with >= sign , then reverse the inequality
        for i=1:m
            if b(i)<0 && inq(i)>=0
                inq(i) = -inq(i);
                b(i) = -b(i);
                for j=1:n
                    A(i,j) = -A(i,j);
                end
            end
        end
        
        % Check for negative RHS
        for i=1:m
            if b(i)<0
                waitfor(msgbox('Primal simplex cannot have negative RHS !','Error','warn'));
                primal_modify(g_modify,v_org,cons_org);
                return
            end
        end
        
        % Check for inequality other than <=
        for i=1:m
            if inq(i)>=0
                waitfor(msgbox('Primal simplex can have only <= inequality !','Error','warn'));
                primal_modify(g_modify,v_org,cons_org);
                return
            end
        end
        
        A1 = horzcat(A,eye(m));
        global A2;
        A2 = A1;
        c1 = horzcat(c,zeros(1,m));
        iter = 0;
        
        B_indices = (n+1:n+m) ;         %% indices of Basic Variables
        V_indices = (1:n) ;             %% indices of Non-Basic Variables
        
        if minimize==1               %% For minization problem
            c1 = -c1 ;
        end
        global c2;
        %% Beginning While loop
        while (1==1)
            iter = iter + 1;
            
            B  = A1(:,B_indices);
            Aj = A1(:,V_indices);
            Cb = c1(:,B_indices);
            Cj = c1(:,V_indices);
            
            bstar = B\b ;           %% B\b = inv(B)*b
            
            c2 = c1;
            
            f = Cb*(B\b) + f1;                                        %% Calculate Final Optimal Solution

            if minimize == 1
                f = -f;
            end
            
            cols = n + m;
            
            for i=1:cols
                size_B = size(find(B_indices == i));
                if size_B(2) ~= 0              %% IF BASIC VARIABLE
                    index_B = find(B_indices ==  i);
                    A2(:,i) = 0;
                    A2(index_B,i) = 1;
                    c2(1,i) = 0;
                else                           %% IF NON-BASIC VARIABLE
                    A2(:,i) = B\A1(:,i) ;
                    index_V = find(V_indices == i);
                    
                    c2(i) = Cb*(B\Aj(:,index_V)) - Cj(index_V) ;
                end
            end

            %****************************** Storing values in table *******************************************

            if minimize == 1
                rowname1{iter}{1} = 'z(min)';
            else
                rowname1{iter}{1} = 'z(max)';
            end
            
            for i=2:m+1
                if B_indices(i-1) <= v
                    if U_indices_new(B_indices(i-1)) == 0
                        rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))));
                    elseif U_indices_new(B_indices(i-1)) == 1
                        rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))),'+'); 
                    elseif U_indices_new(B_indices(i-1)) == 2
                        rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))),'-');     
                    end    
                else
                    rowname1{iter}{i} = strcat('sx',num2str(B_indices(i-1)));
                end
            end
            
            
            dat{iter} = [];
            dat1{iter}{1,1} = '1';
            for i=2:m+1
                dat1{iter}{i,1} = '0';
            end
            
            for i=2:m+n+1
                dat1{iter}{1,i} = num2str(c2(i-1));
            end
            
            dat1{iter}{1,m+n+2} = num2str(f);
            
            for i=1:m+n               %% Printing A2
                for j=1:m
                    dat1{iter}{j+1,i+1} = num2str(A2(j,i)) ;
                end
            end
            
            for i=1:n                   %% Check for lower bound and change the bstar values
                l_temp = size(find(B_indices==i));
                if l_check(i) == 5 && l_temp(2)
                    index_B = find(B_indices ==  i);
                    bstar(index_B) = bstar(index_B) + l_bound(i);
                end    
            end
            
            for i=1:m                 %% Printing bstar
                dat1{iter}{i+1,m+n+2} = num2str(bstar(i));
            end
            
            
            %****************************** Done storing values in table *************************************************
            
            if u == -999
                return
            end
        
            ent_var_mat = Cb*(B\Aj) - Cj;
            
            if min(ent_var_mat) >= 0                                      %% Optimality condition satisfied
                
                [B_indices B bstar iter A2 c2 dat1 rowname1 f U_indices U_indices_new U_id l_check l_bound] = alternate_check(A1,B, bstar, B_indices, V_indices, m,ent_var_mat,u,n,A2,minimize,iter,rowname1,dat1,c2,f,b,c1,U_indices,U_indices_new,U_id,l_check,l_bound);            %% Check for ALTERNATE SOLUTIONS
                
                final_print(f,minimize,A2,c2,bstar,m,n,v,cons,B_indices,iter,rowname1,U_indices,U_indices_new,U_id);                                   %% Call FINAL PRINT function
                return
            else
                
                for i=1:m+n+2                 %% Printing BLANK ROW
                    dat1{iter}{m+2,i} = '' ;
                end
                
                rowname1{iter}{m+2} = strcat('******** Iteration No. ',num2str(iter+1),'  *********');
                               
                ent_var_index = find(ent_var_mat == min(ent_var_mat));
                ent_var_index = V_indices(ent_var_index(1));                %% Entering Variable
                
                [A1 B bstar ent_var_index B_indices V_indices m u rowname1 iter U_indices U_indices_new U_id] = subroutine(A1,B,bstar,ent_var_index,B_indices,V_indices,m, u,n,A2,f,minimize,c2,cons,iter,rowname1,U_indices,U_indices_new,U_id);    %% Calling the SUBROUTINE function
            end
        end         %% End of while loop
    end         %% End of primal_solve

%% Final PRINT function
    function final_print(f,minimize,A2,c2,bstar,m,n,v,cons,B_indices,iter,rowname1,U_indices,U_indices_new,U_id)
        
        % ************** PRINTING OPTIMAL TABLE BEGINS HERE **********************************************************************
        T = figure('Name','Simplex Table Values','NumberTitle','off','Position',[50 50 1200 600],...
            'DeleteFcn',@Figure_Close_Callback);
        columnname = {};
        columnformat = {};       
        columnname{1} = 'z';    %% First column name is 'z'
        
        j=2;
        for i=2:v+1               %% Columnames for decision variables
            if U_indices(i-1) == 0
                str1 = strcat('x',num2str(i-1));
                columnname{j} = str1 ;
                j=j+1;
            else
                str1 = strcat('x',num2str(i-1),'+');
                str2 = strcat('x',num2str(i-1),'-');
                columnname{j} = str1 ;
                columnname{j+1} = str2 ;
                j = j+2;
            end    
        end
        
        for i=v+2:m+n+1           %% Columnnames for slack variables
            str1 = strcat('sx',num2str(i-1));
            columnname{i} = str1 ;
        end
        
        columnname{m+n+2} = 'Solution';
        
        if mod(iter,2) == 0     %% If even iterations
            for i=1:2:iter
                xyz = vertcat(dat1{i},dat1{i+1});
                daddy_dat = vertcat(daddy_dat,xyz);
                
                lol = horzcat(rowname1{i},rowname1{i+1});
                daddy_row = horzcat(daddy_row,lol);
            end
            
        else
            if iter~=1
                for i=1:2:iter-1
                    xyz = vertcat(dat1{i},dat1{i+1});
                    daddy_dat = vertcat(daddy_dat,xyz);

                    lol = horzcat(rowname1{i},rowname1{i+1});
                    daddy_row = horzcat(daddy_row,lol);
                end

                daddy_dat = vertcat(daddy_dat,dat1{iter});
                daddy_row = horzcat(daddy_row,rowname1{iter});
            else
                daddy_dat = dat1{1};
                daddy_row = rowname1{1};    
            end    
        end
        
        
        dat1;
        q = uitable('Parent',T,...
            'Position',[0 0 1200 600],'Data', daddy_dat,...
            'ColumnName', columnname,...
            'ColumnEditable', false,...
            'FontSize',10,'ForegroundColor','k', ...
            'FontName','Comic Sans MS', ...
            'RowName',daddy_row );
        
        T;
         %*************** PRINTING OPTIMAL TABLE ENDS HERE **********************************************************************
        
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
%                     close all; clear all;
                    return  
                    
                case 'View/Modify Input Data'   
                    primal_modify(g_modify,v_org,cons_org);
            end       
    end

%% ALTERNATE SOLUTION function
    function [B_indices B bstar iter A2 c2 dat1 rowname1 f U_indices U_indices_new U_id l_check l_bound] = alternate_check(A1,B, bstar, B_indices, V_indices, m,ent_var_mat,u,n,A2,minimize,iter,rowname1,dat1,c2,f,b,c1,U_indices,U_indices_new,U_id,l_check,l_bound)
        var = size(ent_var_mat);
        length = var(2);
        check = 0;
        
        for temp=1:length
            if ent_var_mat(temp)==0 
                x = U_id(V_indices(temp));
                q = find(U_id(B_indices)==x);
                q = size(q);
                if q(2) == 0
                ent_var_index = V_indices(temp);
                waitfor(msgbox('Optimal Solution Found. ALTERNATE SOLUTION Exists ! Click OK to proceed '));
                          
                for i=1:m+n+2                 %% Printing BLANK ROW
                    dat1{iter}{m+2,i} = '' ;
                end
                
                rowname1{iter}{m+2} = strcat('ALTERNATE OPTIMAL SOLUTION ');
                
                iter = iter + 1
               
                [A1 B bstar ent_var_index B_indices V_indices m u rowname1 iter U_indices U_indices_new U_id] = subroutine(A1,B,bstar,ent_var_index,B_indices,V_indices,m,u,n,A2,f,minimize,c2,cons,iter,rowname1,U_indices,U_indices_new,U_id);
                
                B  = A1(:,B_indices);
                Aj = A1(:,V_indices);
                Cb = c1(:,B_indices);
                Cj = c1(:,V_indices);
                
                bstar = B\b ;           %% B\b = inv(B)*b
                
                c2 = c1;
                
                f = Cb*(B\b);                                             %% Calculate Final Alternate Optimal Solution
                if minimize == 1
                    f = -f;
                end
                
                cols = n + m;
                
                for i=1:cols
                    size_B = size(find(B_indices == i));
                    if size_B(2) ~= 0              %% IF BASIC VARIABLE
                        index_B = find(B_indices ==  i);
                        A2(:,i) = 0;
                        A2(index_B,i) = 1;
                        c2(1,i) = 0;
                    else                           %% IF NON-BASIC VARIABLE
                        A2(:,i) = B\A1(:,i) ;
                        index_V = find(V_indices == i);
                        
                        c2(i) = Cb*(B\Aj(:,index_V)) - Cj(index_V) ;
                    end
                end

                %****************************** Storing values in table *******************************************

                if minimize == 1
                    rowname1{iter}{1} = 'z(min)';
                else
                    rowname1{iter}{1} = 'z(max)';
                end
                
                for i=2:m+1
                    if B_indices(i-1) <= v
                        if U_indices_new(B_indices(i-1)) == 0
                              rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))));
                        elseif U_indices_new(B_indices(i-1)) == 1
                            rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))),'+'); 
                        elseif U_indices_new(B_indices(i-1)) == 2
                            rowname1{iter}{i} = strcat('x',num2str(U_id(B_indices(i-1))),'-');     
                        end    
                    else
                        rowname1{iter}{i} = strcat('sx',num2str(B_indices(i-1)));
                    end
                end
               
                dat1{iter}{1,1} = '1';
                for i=2:m+1
                    dat1{iter}{i,1} = '0';
                end
               
                for i=2:m+n+1
                    dat1{iter}{1,i} = num2str(c2(i-1));
                end
                              
                dat1{iter}{1,m+n+2} = num2str(f);
                
                for i=1:m+n               %% Printing A2
                    for j=1:m
                        dat1{iter}{j+1,i+1} = num2str(A2(j,i)) ;
                    end
                end
                
                for i=1:n                   %% Check for lower bound and change the bstar values
                    l_temp = size(find(B_indices==i));
                    if l_check(i) == 5 && l_temp(2)
                        bstar(i) = bstar(i) + l_bound(i)
                    end    
                end
            
                
                for i=1:m                 %% Printing bstar
                    dat1{iter}{i+1,m+n+2} = num2str(bstar(i));
                end
                
                
                %****************************** Done storing values in table*************************************************
                return
            else   
                check = check + 1;
                end    
            end
        end
        
        if check==length
            waitfor(msgbox('Optimal Solution Found ! Click OK to proceed '));
            return
        end
        return
    end

%% SUBROUTINE function
    function [A1 B bstar ent_var_index B_indices V_indices m u rowname1 iter U_indices U_indices_new U_id] = subroutine(A1,B,bstar,ent_var_index,B_indices,V_indices,m,u,n,A2,f,minimize,c2,cons,iter,rowname1,U_indices,U_indices_new,U_id)
        piv_col = B\A1(:,ent_var_index);
        
        ratio = zeros(m,1);
        for i=1:m
            if piv_col(i)>0                              %% Denominator positive
                if bstar(i)>=0                           %% Numerator non-negative
                    ratio(i) = bstar(i)/piv_col(i);      %% Finding minimum ratio
                else
                    ratio(i) = inf;
                end
            else
                ratio(i) = inf;
            end
        end

        k = 0;
        for i=1:m
            if (ratio(i) == inf)       %% If cannot find ratio
                k = k+1;
            end
        end
        
        if k == m
            rowname1{iter}{m+2} = 'UNBOUNDED SOLUTION';
            waitfor(msgbox('The given problem has UNBOUNDED SOLUTION ! Click OK to proceed '));
            final_print(f,minimize,A2,c2,bstar,m,n,v,cons,B_indices,iter,rowname1,U_indices,U_indices_new,U_id)             %% Call FINAL PRINT function
         
            u = -999;
              return
        end
        
        %   If ratio exists, find the minimum ratio
        min_ratio = min(ratio);
        min_ratio_index = find(ratio == min(ratio)) ;       %% Leaving Variable Index
        leave_var_index = B_indices(min_ratio_index(1));
        
        temp = find(B_indices == leave_var_index);
        B_indices(temp) = ent_var_index;
        temp1 = find(V_indices == ent_var_index);
        V_indices(temp1) = leave_var_index;
        return
    end
%% Check Function For Bounds
 function [temp1 temp2] = check(temp1,temp2,g,er)
                if str2num(temp1) > str2num(temp2)
                        str = strcat('x',num2str(er));
                        waitfor(msgbox('The Lower bound is greater than Upper Bound ! Re-Enter the Values.',str,'warn'));
                        
                        prompt = {'Enter the Lower Bound : ','Enter the Upper Bound :'};
                        dlg_title = str;
                        num_lines = 1;
                        def = {'0','Inf'};
                        answer = inputdlg(prompt,dlg_title,num_lines,def);

                        temp1 = answer{1};
                        temp2 = answer{2};
                        
                        [temp1 temp2] = check(temp1,temp2,g,er);
                end            
 end     
end