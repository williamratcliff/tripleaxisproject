function data=readfiles(files,filetype,scantype)

if nargin < 3
filedir='c:\ZnMn2O4\jan2007\';
filehead='znmnescan424';
fileend='.bt7';
%scans=[42588];
scans=[58,65,70,76,77,85];
iceflag=1;
files=generatefiles(filedir,filehead,fileend,scans,iceflag)
filetype='icp'
%files={'c:\ZnMn2O4\Jan2006\mag10014.bt9'}; %a3scan 
files={'c:\ZnMn2O4\bt9\Oct13_2005\mp20k012.bt9'};
scantype='q'
end

if filetype=='ice'
   data=readbt7generalmy(files); 
end
if filetype=='icp'
        data=[];
        if scantype=='a3'
                for i=1:length(files)
                filename=char(files(i));
                datum=[];
                [datum,Qm, mon]=nistload_general_noprompt(filename);
                %datum
                data(i).T=mean(datum(:,2));
                data(i).Qm=Qm;
                data(i).monitor=mon;
                data(i).I=datum(:,4);
                data(i).Ierr=sqrt(datum(:,4));
                data(i).a3=datum(:,1);
                end
        end
        if scantype=='q'
                for i=1:length(files)
                filename=char(files(i));
                [datum, header, column_labels, flip]=loadnistdata(filename);
                %data(i).data=datum;
                data(i).qx=datum(:,1);
                data(i).qy=datum(:,2);
                data(i).qz=datum(:,3);
                data(i).E=datum(:,4);
                data(i).I=datum(:,6);
                data(i).Ierr=sqrt(datum(:,6));
                data(i).T=header.temp.avg; %mean(datum(:,5));
                data(i).center=sscanf(header.center(2:end-1),'%f %f %f');
                data(i).monitor=header.monitor;
                end
        end
end

return
