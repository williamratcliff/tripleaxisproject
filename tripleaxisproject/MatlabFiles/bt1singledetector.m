function bt1singledetector(detector);
if nargin~=1
    detector=10;
end

clear global;
global T
global M
global Merr


currtemp=[];
currdetector=[];
files={'C:\ZnMn2O4\5Kznm003.bt1';'C:\ZnMn2O4\5Kznm004.bt1';'C:\ZnMn2O4\5Kznm005.bt1'};
files={'C:\ZnMn2O4\5Kznm003.bt1';'C:\ZnMn2O4\5Kznm004.bt1'};




for j=1:length(files)
fid=fopen([char(files(j))],'r');
 for i=1:17
                line = fgetl(fid);if ~isstr(line), break, end
  
 end
 
while 1
           
            line = fgetl(fid);if ~isstr(line), break, end
            if findstr('T=',line)
               token = strtok(line,'C');
               [token2,rem] = strtok(token,'=');
               [token2,rem2] = strtok(rem,' ');
               currtemp=[currtemp str2num(char(rem2))];
            end
            if findstr('$',line)
            line = fgetl(fid);if ~isstr(line), break, end
            end
            %now comes the detectors
            line = fgetl(fid);if ~isstr(line), break, end
            A1 = strread(line,'%d','delimiter', ',');
            line = fgetl(fid);if ~isstr(line), break, end
            A2=strread(line,'%d','delimiter', ',');
            A=[A1' A2'];
            currdetector=[currdetector A(detector)];
            
end

end



start=40;
finish=length(currtemp)-5;
opdata=[currtemp(start:finish);currdetector(start:finish);sqrt(currdetector(start:finish))]';
%plot(opdata(1,:),opdata(2,:))
%return

clf reset;

p=[392 62.3 1.1 0.125 550];

errorbar(opdata(:,1),opdata(:,2),opdata(:,3), 'ko');
hold;
%guess = orderparameter(p,opdata(:,1));
%plot(opdata(:,1), guess, 'r-');
%[BETA,R,Jacobian] = nlinfit(opdata(:,1),opdata(:,2),'orderparameter',p)
%[BETA,CHISQUARE,ERRORS,FITRESULT] = nlfit(opdata(:,1),opdata(:,2),'orderparameter',p,[1 1 1 0 1], opdata(:,3));
[BETA,CHISQUARE,ERRORS,FITRESULT] = nlfit(opdata(:,1),opdata(:,2),'orderparameter',p,[1 1 1 1 1], opdata(:,3));
fit = orderparameter(BETA,2:1:83);
plot(2:1:83, fit);
BETA
ERRORS
text(110, 200, 'Fit Results:');
text(110, 175, 'Tn = 101.3(1.1)');
text(110, 150, 'dTn = 21.3(2.4)');
text(110, 125, 'Beta = 0.25 fixed');

return

%plot(currtemp,currdetector,'o',currtemp,currdetector);
%errorbar(currtemp,currdetector,sqrt(currdetector),'or');

Ttem=currtemp;
Mtem=(currdetector-currdetector(length(currdetector)))/(currdetector(1)-currdetector(length(currdetector)));
Merrtem=sqrt(currdetector)/(currdetector(1)-currdetector(length(currdetector)));

%T=Ttem(50:length(Ttem)-30));
%M=Mtem(50:length(Ttem)-30))
%Merr=Merrtem(50:length(Ttem)-30));

T=Ttem;
M=Mtem;
Merr=Merrtem;



Tc=63
b=.23
pfit=[b Tc]';
lower=[0 0];
upper=[1 70];

pfit=b;
lower=0;
upper=1;

options=optimset('Display','iter','MaxIter',50,'LargeScale','on','TolFun',.001);
%%p=lsqnonlin(@chisqrcalcinline,pfit,lower,upper,options,PointNum,Amps,ErrAmps,s);
%p=lsqnonlin(@chisqrcalcinline,pfit,lower,upper,options);



%Tc=abs(p(2));
%b=abs(p(1))

model=(1-T./Tc).^(b);

errorbar(T,M,Merr,'or'); hold on; plot(T,model); hold off;
xlabel('Temperature');ylabel('Intensity');



fid=fopen(['C:\ZnMn2O4\orderparameter.dat'],'w');
fprintf(fid,'%3.5g \t %3.5g \t %3.5g \r',...
    [T;M;Merr]);
fclose(fid)

function chisq=chisqrcalcinline(pfit)

global T
global M
global Merr
Tc=63;
b=abs(pfit(1));

%file='magnstrDyMn2O5incomremnk1k2_2';
%file='magnstrDyMn2O5incomremnloose2';
%8 is what we were using
%file='realspacestruct9'; %9 is pretty broad, allows for all planar structures
%file='realspacestruct10'; %10 we use everything, including 3 dimensions
%file='realspacestruct14'; %10 we use everything, including 3 dimensions


%model=feval(file,PointNum,s,pfit,scA1,scA2,frac);
%model=feval(file,PointNum,s,pfit);

model=(1-T./Tc).^b;

chisq=(M-model)./Merr;


return