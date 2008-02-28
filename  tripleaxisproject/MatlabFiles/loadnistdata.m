function [data, header, column_labels, flip]=loadnistdata(filename)

%LOADNISTDATA loads NIST data file [NEUTRON]
%   [data, header, column_labels, flip]=loadnistdata(filename)
%
%Loads a data set and header read from neutron data files produced by the ICP 
%program at NIST.  Also produces column heading labels and information
%about the polarization flipper states if available.
%
%Modified: 05/03/01 by Patick K. Mang
%
%%%%% Header Structure %%%%%%%%%%%
%
% filename
% date
% monitor.level
%    "   .neutrons
%    "   .time
% polarization.used
%    "        .flipper1
%    "        .flipper2
% description
% lattconst.a
%    "     .b
%    "     .c
% energy.energy
%    "  .KIfixed
%    "  .KFfixed
% center
% temp.avg
%   " .std
%
%In the case of loading an fp or fpt file the header structure is 
%modified to be the gaussian fitting parameters produced by the NIST
%computer in the order of [Amplitude Center FWHM Background]


%-------------- Initialize Arays ----------------
data = [];
column_labels = [];
	temp = [];
header = [];



%-------------- Open Data File -----------------
[fid, message] = fopen(filename, 'r');
if (fid < 0)
   fprintf(1, 'ERROR on %s oepn: %s\n', filename, message);
   return;
elseif strncmp(filename, 'fp', 2)
   %run the alternate load routine if the file is of the type created
   %when scanning one motor, not executing a buffer.  These are distinguished
   %by starting with 'fp'.  ex: fpx03002.bt9
   
   [data,header]=alternateloadroutine(fid);
   return;
   
 elseif ~(strncmp(filename(end-2:end-1), 'bt',2) | ...
         (strncmp(filename(end-2:end-1), 'ng',2)))
 %run alternate load routine two if the files have already been ...
 %alterated during analysis from the raw data format.  This is ...
 %indicated by their not ending in bt# or ng5.

   data=alternateloadroutine2(fid);
   
   return;
   
   
   
end

%-------------- Read Column Headers ----------------------
text = fgetl(fid);
fpos=ffind(filename, 'Q(x)');
fseek(fid, fpos, 'bof');
r=fgetl(fid);
temp = [];
while r	
	[a r] = strtok(r);  %parse into regions separted by "'"
	temp = strvcat(temp, a);
end
column_labels = temp;


%-------------- Read Data ------------------------
r=fgetl(fid);
while length(r) > 2
	a = sscanf(r, '%f');
	data = [data ; a' ];	
	r = fgetl(fid);
end

%------------ Build Scan Header ----------------
	fseek(fid, 0, -1);   %return to beginning of file

r=fgetl(fid); 	     %get first line of header
if r(1) == '#'
	r=r(2:length(r)); %this line is necessary to remove # signs sometimes
end			  %inserted to comment out the header information for 
temp = [];		  %the program C-Plot.

while r	
	[a r] = strtok(r,'\''');  %parse into regions separted by "'"
       temp = strvcat(temp, a);
end

header.filename = deblank(temp(1,:));	%filename
header.date = deblank(temp(3,:));	%time stamp
[mon pref] = strtok(temp(6,:), '.');    %set monitor
pref = pref(2:length(pref));
header.monitor.level = str2num(mon) * str2num(pref);
header.monitor.neutrons=0; header.monitor.time=0;

if (deblank(temp(7,:))=='NEUT')
	header.monitor.neutrons = 1;
else
	header.monitor.time = 1;
end


r=fgetl(fid);	     %skip over line 2


r=fgetl(fid);        %line 3
i=findstr('F1:', r);  %check if polarization mode used

if ~isempty(i)       %if polarized
	flipper1 = deblank(r(i+5,7));    %format is F1: ON or F1: OFF
	i=findstr('F2:', r);
	flipper2 = deblank(r(i+5,7));
	header.polarization.used = 1;
	flip = 1;
	if (flipper1=='ON')
		header.polarization.flipper1 = 1;
	else
		header.polarization.flipper1 = 0;
	end
	if (flipper2 == 'ON')
		header.polarization.flipper2 = 1;
	else
		header.polarization.flipper2 = 0;
	end
	r = strtok(r, 'F1:');	%strip off polarization information
else
	header.polarization.used = 0;
	flip = 0;
end

header.description = deblank(r);	%set scan description



r=fgetl(fid);	%line 4
r=fgetl(fid);	%line 5

r=fgetl(fid);        %get lattice constants (units:angstrom)
if (r(1) == '#')
	r=r(2:length(r));
end
a = sscanf(r, '%f');
%header.lattconst = sprintf('a=%0.3f  b=%0.3f  c=%0.3f', a(1), a(2), a(3));
header.lattconst.a = a(1);
header.lattconst.b = a(2);
header.lattconst.c = a(3);

r=fgetl(fid);	%line 7

r=fgetl(fid);        %get energy (units:meV)
if(r(1) == '#')
	r=r(2:length(r));
end
a = sscanf(r, '%f');
header.energy.energy = a(3);

r=fgetl(fid);	%determine if initial or final energy is fixed
header.energy.KIfixed = 0; header.energy.KFfixed = 0;
if findstr('EM fixed', r) header.energy.KIfixed = 1; end
if findstr('EA fixed', r) header.energy.KFfixed = 1; end

r=fgetl(fid);	     %get scan center in reciprocal space
if(r(1) == '#')
	r=r(2:length(r));
end
a=sscanf(r, '%f');
header.center = sprintf('(%0.3f  %0.3f  %0.3f)', a(1), a(2), a(3));

header.temp.avg = mean(data(:,5));
header.temp.std = std(data(:,5));

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the alternate load routine for find peak files
function [data,header]=alternateloadroutine(fid)
data=[];
header=[];
temp =[];

line=fgetl(fid);
while line	
  [a line] = strtok(line);  %parse by whitespace
  temp = strvcat(temp, a);
end
num_param = size(temp);
if (num_param(1)==9)
  temp=str2num(temp([3 5 7 9],:));
  header=temp([3 1 4 2]);  %header is now the gaussian fit parameters
%                          %produced by the NIST computers
  for i=1:2
     line=fgetl(fid);      %if a header exists then you need to skip 2
  end                      %lines to get to the data

else
  header=[0 1 1 1];
end

line=fgetl(fid);
while (ischar(line))
   a = sscanf(line, '%f');
   data = [data; a'];
   line=fgetl(fid);
end
fclose(fid);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%alternate load routine two for modified data files
function data=alternateloadroutine2(fid)
data=[];

r=fgetl(fid);
while length(r) > 2
	a = sscanf(r, '%f');
	data = [data ; a' ];	
	r = fgetl(fid);
end
fclose(fid);









