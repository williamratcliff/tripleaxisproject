function [data, Q, mon]=nistload_general_noprompt(file);
% This is a basic load routine for MFIT, illustrating the required
% syntax. The routine takes the name of a data file (including path) as a
% parameter and returns the column vectors x, y, err.
% On error the routine must exit with datafile=0.

%------------- Open data file---------------------
x=[]; y= []; err=[]; ylab=''; xlab='';
fid=fopen(file,'r');				
if (fid<0) 
	datafile=0;
	return
end
%-------------- Initialize arrays-----------------
data=[];Q=[];
%============= Load data ===========================================
%------Read header and first row of data ---------------------------
Q=zeros(1,3);
while 1
   line = fgetl(fid);
   if findstr('RAW',line)
      mon=str2num(line(39:50))*str2num(line(52:55));
   end
   if (findstr('E center  Delta E',line))
      line=fgetl(fid);
      Scan=str2num(line);
      Q(1)=Scan(1);Q(2)=Scan(2);Q(3)=Scan(3);Q(4)=Scan(4);Q=Q';
   end
   if (findstr('Q (hkl scan center)',line)) 
      line=fgetl(fid);
      if findstr('  1  ',line)
         while 1 
            line=fgetl(fid);
            if findstr(' Mot:',line) 
               line=fgetl(fid);
				  head=line;
               break 
            end 
         end
      end
      head=line;
      break
   elseif (findstr('Mot:',line))
      line=fgetl(fid);
      head=line;
      break
   elseif (findstr('Motor no.',line))
      head=['Angle' line(13:16) 'Counts'];
      %disp(head);
      break      
   end
   if ~isstr(line), 
      disp('Wrong format: File cannot be loaded');
      break, 
   end
end
%----Process Header to extract labels ---------------------
datastr = [];p = 1;
lh = length(head);
while p < lh
	letpos = find(isletter(head(p:lh)));
	if ~isempty(letpos)
		istart = letpos(1)+p-1;
		spapos = find(isspace(head(istart:lh)));
		if ~isempty(spapos)
			iend = spapos(1)+istart-2;
         toadd = head(istart:iend);
			if ~isempty(toadd) & find(isletter(toadd))
				if isempty(datastr)
					datastr = toadd;
            else
					datastr = str2mat(datastr,toadd);
				end
				p = iend+1;
			end
		else
         iend = lh;
         toadd = head(istart:iend);
         datastr = str2mat(datastr,toadd);
         p = lh;
  		end
	else
		p = lh;
	end
end
%------ Read data and reshape into matrix ----------------------------
data=fscanf(fid,'%f');            % Read data into vector (for speed)
[ncol dumb]=size(datastr);
nrow=length(data)/ncol;
if (nrow*ncol~=length(data))                 
	error('Bad data format');
end
data=data';data=[reshape(data,ncol,nrow)'];
fclose(fid);                                 % close input file