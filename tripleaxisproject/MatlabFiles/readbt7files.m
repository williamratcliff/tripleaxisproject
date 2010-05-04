function data=readbt7files(filedir,filehead,fileend,scans)


data=[];
for i=1:length(scans)
    file=[filedir filehead num2str(scans(i)) fileend]
    data=[data readbt7file(file)];
end

return
    
    