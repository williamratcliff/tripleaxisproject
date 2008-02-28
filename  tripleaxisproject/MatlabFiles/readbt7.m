function data=readbt7()


%%%%read in files

clear global

mon0=4859*5;


%read in bt7 files


filedir='c:\ZnMn2O4\jan2007\';
filehead='znmn250kscan';
fileend='.bt7';
%scans=[42588];
scans=[42613];
data=[];
for i=1:length(scans)

    file=[filedir filehead num2str(scans(i)) fileend]

    if exist([file])
        fid=fopen([file],'r');
    %    disp('hello')
    end
    
%    [datum]=fscanf(fid,'%f %f',[2 inf]);fclose(fid);
%    data(i).datum=datum;
%end
%end

data(i).qx=[];
data(i).qy=[];
data(i).I=[];
data(i).mon=[];
for j=1:28
            line = fgetl(fid);if ~isstr(line), break, end  
end
 
while 1           
            line = fgetl(fid);if ~isstr(line), break, end
%            line
%            i
            if line(1)=='#', break, end
            remain=line;
            A={};
            while true
                [tok, remain] = strtok(remain);
                A=[A;tok];
                if isempty(tok),  break;  end
            end
%            sscanf(char(A(1)),'%f')

            %A = strread(line,'%f');
            %A2=strread(line,'%d','delimiter', ',');
            data(i).qx=[data(i).qx sscanf(char(A(1)),'%f')];
            data(i).qy=[data(i).qy sscanf(char(A(2)),'%f')];
            data(i).I=[data(i).I sscanf(char(A(9)),'%f')];
            data(i).mon=[data(i).mon sscanf(char(A(8)),'%f')];
end
fclose(fid);
end

%size(data)

I=data(1).I;
qx=data(1).qx;
qy=data(1).qy;
Ierr=sqrt(I);

return
%size(I)


%errorbar(qy,I,Ierr,'or','color','red'); 

%xmin=26;%40;
%xmax=35;%50;
%ymin=0;
%ymax=800;%1000;
%axis([xmin xmax ymin ymax])



