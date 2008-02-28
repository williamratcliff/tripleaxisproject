function data=readbt7(files)


%%%%read in files

clear global

%mon0=4859*5;
%read in bt7 files
%nargin
if nargin==0
filedir='e:\ZnMn2O4\jan2007\';
filehead='znmn250kscan';
fileend='.bt7';
%scans=[42588];
scans=[42613];
files={};
for i=1:length(scans)
    currfile=[filedir filehead num2str(scans(i)) fileend];
    files=[files; currfile];
end


end

%files
data=[];
for i=1:length(files)

    
    file=char(files(i));
    if exist([file])
        fid=fopen([file],'r');
    %    disp('hello')
    end
    
data(i).qx=[];
data(i).qy=[];
data(i).qz=[];
data(i).I=[];
data(i).mon=[];
data(i).E=[];
for j=1:28
            line = fgetl(fid);if ~isstr(line), break, end  
end
 
while 1           
            line = fgetl(fid);if ~isstr(line), break, end
            if line(1)=='#', break, end
            remain=line;
            A={};
            while true
                [tok, remain] = strtok(remain);
                A={A{:},tok};
                if isempty(tok),  break;  end
            end
%            sscanf(char(A(1)),'%f')

            %A = strread(line,'%f');
            %A2=strread(line,'%d','delimiter', ',');
            data(i).qx=[data(i).qx sscanf(char(A(1)),'%f')];
            data(i).qy=[data(i).qy sscanf(char(A(2)),'%f')];
            data(i).qz=[data(i).qz sscanf(char(A(3)),'%f')]; 
            data(i).E=[data(i).E sscanf(char(A(4)),'%f')];
            data(i).I=[data(i).I sscanf(char(A(9)),'%f')];
            data(i).Ierr=sqrt(data(i).I);
            data(i).mon=[data(i).mon sscanf(char(A(8)),'%f')];
end
data(i).qx=data(i).qx';
data(i).qy=data(i).qy';
data(i).qz=data(i).qz';
data(i).E=data(i).E';
data(i).I=data(i).I';
data(i).Ierr=data(i).Ierr';
data(i).mon=data(i).mon';
fclose(fid);
end


%I=data(1).I;
%qx=data(1).qx;
%qy=data(1).qy;
%Ierr=sqrt(I);
%errorbar(qy,I,Ierr,'or','color','red'); 
return





%xmin=26;%40;
%xmax=35;%50;
%ymin=0;
%ymax=800;%1000;
%axis([xmin xmax ymin ymax])



