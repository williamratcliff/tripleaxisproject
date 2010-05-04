function files=generatefiles(filedir,filehead,fileend,scans,iceflag)

if nargin < 5
filedir='c:\ZnMn2O4\jan2007\';
filehead='znmnescan424';
fileend='.bt7';
%scans=[42588];
scans=[58,65,70,76,77,85];
iceflag=0;
end

if iceflag==1
files={};
for i=1:length(scans)
    currfile=[filedir filehead num2str(scans(i)) fileend];
    files={files{:}, currfile};
end

else
    fileszero='00';
    fileszero2='0';
    fileszero3='';
   
    files={};
    scans1=scans(find(scans<10));
    scans2=scans(find((scans>=10 & scans <100)));
    scans3=scans(find(scans>=100));
    for i=1:length(scans1)
    currfile=[filedir filehead fileszero num2str(scans1(i)) fileend];
    files={files{:},currfile};
    end
    for i=1:length(scans2)
    currfile=[filedir filehead fileszero2 num2str(scans2(i)) fileend]
    files={files{:},currfile};
    end
    for i=1:length(scans3)
    currfile=[filedir filehead num2str(scans2(i)) fileend]
    files={files{:},currfile};
    end

end



return