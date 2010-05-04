function data=readbt7file(file)

%filedir is the file directory:
%ex
%filedir='c:\14213\data\';
%file heading for all of the files:
%filehead='split68';
%file ending:
%fileend='.bt7';
%scan range
%scans=[497,498];
file
%data=[];
%for i=1:length(scans)

    %file=[filedir filehead num2str(scans(i)) fileend]

    if exist([file])
        fid=fopen([file],'r');
   
    end
    

line = fgetl(fid);
remain=line;
           

            
            
            

while 1           
            line = fgetl(fid);if ~isstr(line), break, end

            remain=line;
             A={};
            j=1;
            while true
                [tok, remain] = strtok(remain); 
                if isempty(tok),  break;  end
                A{j}=tok;
                j=j+1;
            end
            if strcmp(A{1},'#Columns'), break; end
            
           
end


for j=2:length(A)
data.(lower(A{j}))=[];    
    
end

fields=fieldnames(data);





while 1           
            line = fgetl(fid);if ~isstr(line), break, end

            remain=line;
            A={};
            j=1;
            while true
                [tok, remain] = strtok(remain); 
                if isempty(tok),  break;  end
                A{j}=tok;
                j=j+1;
            end
            
    
            for j=1:length(fields)
               currvalue=str2num(A{j});
               currfield=char(fields(j));
               
               if length(currvalue)==1
                  data.(currfield)=[data.(currfield) currvalue];
               else
                   data.(currfield)=[data.(currfield) A(j)];
               end
                   
            end
            
end

%data(1).a4
%data(1).detector
fclose(fid);

%plot(data(i).a4,data(i).detector,'s'); hold on
%end



return