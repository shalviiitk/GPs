function [s, dsde, lemda,svep, errf] = radial_return(s,lemdan,epn,dstrain,props,s0,harden)
% where s is stress
% lemdan is plastic multiplier
% epn is reference plastic strain
% s0 is yield strength
%Plane strain
errf=0;
iter=1;
tol=1e-10;
ex=props(1)*10^9;
nu=props(2);

dsde(1:3,1:3)=0;

dstrain3d(1:3,1:3)=0;

dstrain3d(1,1)=dstrain(1);
dstrain3d(2,2)=dstrain(2);
dstrain3d(1,2)=dstrain(3)/2;
dstrain3d(2,1)=dstrain(3)/2;

dstrainh=(dstrain3d(1,1)+dstrain3d(2,2)+dstrain3d(3,3))/3;
dstraindev=dstrain3d-dstrainh*eye(3);     %dev strain is equal to sum of elastic and plastic strain

xk=ex/(3*(1-2*nu));
xg=ex/(2*(1+nu));

s3d(1:3,1:3)=0;
s3d(1,1)=s(1);
s3d(2,2)=s(2);
s3d(3,3)=s(3);
s3d(1,2)=s(4);
s3d(2,1)=s(4);

ds=2*xg*dstraindev+3*xk*dstrainh*eye(3);
% ds=2*xg*dstraindev;
str=s3d+ds;                   %fixed trial stress or Elastic predictor

sh=(str(1,1)+str(2,2)+str(3,3))/3;
sdev0=str-sh*eye(3);

seqv0=0;

for i=1:3
    for j=1:3
        seqv0=seqv0+sdev0(i,j)*sdev0(i,j);   %??
    end
end
seqv0=sqrt(1.5*seqv0);

if abs(seqv0) > 0
    r0=1.5*sdev0/seqv0;
else
    r0=zeros(3,3);
end

ep(1,1)=epn(1);
ep(2,2)=epn(2);
ep(3,3)=epn(3);
ep(1,2)=epn(4);
ep(2,1)=epn(4);


dl=0;           % delta lemda

[sy, dsydq] = calc_yield_stress(lemdan,s0,harden);

f=seqv0-sy;
% f=0;
maxiter=50;
if f > 0
    
   while abs(f) > tol && iter <= maxiter
     
       ddl=(seqv0-3*xg*dl-sy)/(3*xg+dsydq);         % del lemda  
       dl=dl+ddl;
       s3d=str-2*xg*dl*r0;                          %Elastic predictor minus plastic corrector
       lemda=lemdan+dl;
       [sy, dsydq] = calc_yield_stress(lemda,s0,harden);
       f=seqv0-3*xg*dl-sy; 
       iter=iter+1;
       if abs(f) > tol
           errf=1;
           break
       end
   end
  
   upep=ep+dl*r0;
   svep(1)=upep(1,1);
   svep(2)=upep(2,2);
   svep(3)=upep(3,3);
   svep(4)=upep(1,2);
   dsde=calc_jacobian(xk,xg,r0,dl,dsydq,seqv0);
    
else
    
    lemda=lemdan;
    svep=epn; 
    s3d=str;
    dsde(1,1)=ex*(1-nu)/((1+nu)*(1-2*nu));
    dsde(1,2)=ex*nu/((1+nu)*(1-2*nu));
    dsde(2,1)=dsde(1,2);
    dsde(2,2)=dsde(1,1);
    dsde(3,3)=ex/(2*(1+nu)); 
    ds1=dsde*dstrain;
    
    err(1)=ds1(1)-ds(1,1);
    err(2)=ds1(2)-ds(2,2);
    err(3)=ds1(3)-ds(1,2);
%     err(4)=ds1(4)-ds(1,2);
end

s(1)=s3d(1,1);
s(2)=s3d(2,2);
s(3)=s3d(3,3);
s(4)=s3d(1,2);
end