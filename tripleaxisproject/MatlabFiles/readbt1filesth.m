function readbt1files()

clear global
instfile='c:\ZnMn2O4\5Kznm001.inst';
gsasfile='c:\ZnMn2O4\5kznm001.gsas';
gsasfile2='c:\ZnMn2O4\5kznm002.gsas';



%if exist([instfile])
%    fid=fopen([instfile],'r');%thedata=fscanf(fid,'%f %f %f',[3 inf])';fclose(fid);
%end



%while 1
           
%            line = fgetl(fid);if ~isstr(line), break, end
%            if findstr('ITYP',line)
%               [thetarange]=sscanf(line,'%*13c %*f %f %f %*f');
%               starttheta=thetarange(1)
%               endtheta=thetarange(2)
            
%            end
            
%end

%fclose(fid);


if exist([gsasfile])
    fid2=fopen([gsasfile],'r');%thedata=fscanf(fid,'%f %f %f',[3 inf])';fclose(fid);
    
end


datat=[];
while 1
           
            line = fgetl(fid2);if ~isstr(line), break, end
            if findstr('BANK',line)
               %token = strtok(line,'ITYP')
               [nrecords]=sscanf(line,'%*4s %*d %f %f %*5s %f %f' );
               ncols=nrecords(1)/nrecords(2);
               starttheta=nrecords(3)/100;
               steptheta=nrecords(4)/100;
               endtheta=starttheta+(nrecords(1)-1)*steptheta;
               while 1
                  line = fgetl(fid2);if ~isstr(line), break, end
                  [values]=sscanf(line,'%f');
                  datat=[datat values];
               end
               
               
            end
            
end

%data
fclose(fid2);

[n m]=size(datat);
datat=reshape(datat,1,n*m);
tarr=1:n*m;
data=datat(find(mod(tarr,2)~=0));
dataerr=datat(find(mod(tarr,2)==0));
theta=starttheta:steptheta:endtheta;
%length(theta)
%length(data)
%length(dataerr)
%plot(theta,data);

%secondfile



if exist([gsasfile2])
    fid2=fopen([gsasfile2],'r');%thedata=fscanf(fid,'%f %f %f',[3 inf])';fclose(fid);
    
end

nrecords=[];
datat2=[];
while 1
           
            line = fgetl(fid2);if ~isstr(line), break, end
            if findstr('BANK',line)
               %token = strtok(line,'ITYP')
               [nrecords]=sscanf(line,'%*4s %*d %f %f %*5s %f %f' );
               ncols=nrecords(1)/nrecords(2);
               starttheta2=nrecords(3)/100;
               steptheta2=nrecords(4)/100;
               endtheta2=starttheta2+(nrecords(1)-1)*steptheta2;
               theta2=starttheta2:steptheta2:endtheta2;
               while 1
                  line = fgetl(fid2);if ~isstr(line), break, end
                  [values]=sscanf(line,'%f');
                  datat2=[datat2 values];
               end
               
               
            end
            
end

%data
fclose(fid2);

[n m]=size(datat2);
datat2=reshape(datat2,1,n*m);
tarr2=1:n*m;
data2=datat2(find(mod(tarr2,2)~=0));
dataerr2=datat2(find(mod(tarr2,2)==0));

%plot(theta2,data2);



data3=data(1:length(data2))-data2;
dataerr3=sqrt(dataerr(1:length(data2)).^2+dataerr2.^2);
data3=data3;
global lambda;
lambda=1.5403;

%plot(theta3,data3,'s');
%axis([1 5 -500 4000]);

trange=find(theta2>26 & theta2<36);
%plot(theta2(trange),data3(trange));



global thfit;
global dfit;
global derr;

thfit=theta2(trange);
dfit=data3(trange)+105;
derr=dataerr3(trange);

global q
q=4*pi*sin(thfit/2*pi/180)/lambda;

%plot(q,dfit)

pfit=[-100 50 .03 28.2];
lower=[0 -Inf -Inf 27.5];
upper=[5 Inf Inf 29.5];

options=optimset('Display','iter','MaxIter',50,'LargeScale','on','TolFun',.001);
%%p=lsqnonlin(@chisqrcalcinline,pfit,lower,upper,options,PointNum,Amps,ErrAmps,s);
p=lsqnonlin(@chisqrcalcinline,pfit,lower,upper,options);

a=p(1); b=p(2); %c=p(3);
%model=a+b./(q-c).^1;


L=p(3);
thetab=p(4);
dtor=pi/180;

la=(2*L*sqrt(pi)/lambda)*(sin(thfit*dtor)-sin(thetab*dtor));
x=thfit;
F1=sqrt(L/sqrt(pi)/lambda);
for i=1:length(la)
    lai=i;
F2(i)=quad(@fintegral,0.,1.,1e-5);
end

%F2
%return

F=F1*F2;

model=a+(b*F./(sin(thfit*dtor/2)).^1.5)./sin(thfit*dtor);

p

plot(thfit,dfit,'s'); hold on;
plot(thfit,model);
hold off;


%chisq=((dfit-model)./derr);
%chisq=chisq*chisq


return

function chisq=chisqrcalcinline(pfit)

global thfit
global q
global dfit
global derr
global lambda

a=pfit(1);
b=pfit(2);
%c=pfit(3);
global la
global lai
L=pfit(3);
thetab=pfit(4);
dtor=pi/180;

%disp('la')
la=(2*L*sqrt(pi)/lambda)*(sin(thfit*dtor)-sin(thetab*dtor));
x=thfit;

F1=sqrt(L/sqrt(pi)/lambda);


for i=1:length(la)
    lai=i;
%F2(i)=quad(@fintegral,0.,1.,1e-5);
F2(i)=trapz(@fintegral,0.,1.);
end



F=F1*F2;

model=a+(b*F./(sin(thfit*dtor/2)).^1.5)./sin(thfit*dtor);


chisq=(dfit-model)./derr;


return


function F=fintegral(x);
global la
global lai
%disp('reached')
%F=zeros(size(la));
%x
%for i=1:length(la)
    %la(i)
    %(-x.^2.-la(i))
    %exp((-x.^2-la(i)).^2);
F=exp((-x.^2-la(lai)).^2);
return
end