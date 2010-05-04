function summed_scan=addscans(filename_head, filename_tail, scans, normalize)
%ADDSCANS  add together two NIST neutron scans   [NEUTRON]
%
%        [summed_Scan]=addscans(filename_head, filename_tail, scans, normalize)
%            filename_head and filename_tail specify the format of the
%            files to be added together. e.g. 'ncco4' and 'bt2'. Scans is
%            a matrix of the scan numbers to be added. If normalize is
%            specified, all points will be normalized to the same number of
%            minutes and the last column of the return data set will
%            contain appropriately calculated errorbars for each point.
%
%Modified 3/24/03 by Patrick K. Mang

summed_scan=[];
temp=[];
if (nargin <4)
	normalize=0;
end


for i=1:length(scans)
	filename=strcat(filename_head,sprintf('%03d',scans(i)),'.',filename_tail);
	temp=loadnistdata(filename);
	summed_scan=[summed_scan; temp];
end


scanned_var=std(summed_scan(:,1:end));  %find which columns are changing
scanned_var=find(scanned_var > 1e-5);  %index those by column number

if ~isempty(intersect([1 2 3 4], scanned_var))  %if HKLW is being scanned

	%error check that scans are at the same temperature
	if any(abs(summed_scan(:,5)-mean(summed_scan(:,5)))>10)
		error('Temperature variation among scans is greater than 10 degrees!')
	end
	
	scanpts=unique(summed_scan(:,scanned_var(1)));  %identify unique pts
	
	temp=[];  %reinitialize temp
	for i=1:length(scanpts)
		j=find(summed_scan(:,scanned_var(1))==scanpts(i));
		temp=[temp; summed_scan(j(1),1:4) mean(summed_scan(j,5)) sum(summed_scan(j,6:7),1)]; 
		%sum points identical in HKL
	end

	summed_scan=sortrows(temp,scanned_var(1)); 
end


%%%%%%%%%%%%%Normalization%%%%%%%%%%%%%%%%%%%%%%
if (normalize~=0)
	summed_scan(:,8)=sqrt(summed_scan(:,7));
	summed_scan(:,8)=diag(normalize./summed_scan(:,6)*summed_scan(:,8)');
	summed_scan(:,7)=diag(normalize./summed_scan(:,6)*summed_scan(:,7)');
	summed_scan(:,6)=normalize*ones(length(summed_scan),1);
end






