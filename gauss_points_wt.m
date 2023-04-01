function [gauss_points,X1,Y1,weight,p,pt]=gauss_points_wt(nx,order,level,x1,y1,cx,cy)

% First level integration points
if level ==1
    pt=1;
    zai=1/3;
    eta=1/3;
  for i=1:nx/order
    X1(i)=(zai*(x1(i)-cx)+eta*(x1(i+1)-cx)+cx);
    Y1(i)=(zai*(y1(i)-cy)+eta*(y1(i+1)-cy)+cy);
  end
    w=0.5; 
    gauss_points(1:(nx/order)*pt,1:2)=[X1' Y1'];
    weight(1:pt)=w;
    p=[-1 -2 -3 -4 -5 -6];
end

if level==2
pt=3;
 zai=[1/6 2/3 1/6];
 eta=[1/6 1/6 2/3];
%  zai=[1/2 0 1/2];
%  eta=[1/2 1/2 0];
 for i=1:nx/order
     for j=1:length(zai)
     X1_(i,j)=(zai(j)*(x1(i)-cx)+eta(j)*(x1(i+1)-cx)+cx);
%      X1_2(i,j)=(zai(j)*(x2(i)-cx1)+eta(j)*(x1(i+1)-cx1)+cx1);
     Y1_(i,j)=(zai(j)*(y1(i)-cy)+eta(j)*(y1(i+1)-cy)+cy);
%      Y1_2(i,j)=(zai(j)*(y2(i)-cy1)+eta(j)*(y1(i+1)-cy1)+cy1);
     end
 end
X1=reshape(X1_',[1,nx/order*3]);
Y1=reshape(Y1_',[1,nx/order*3]);
w=1/6;
p=-1:nx-2;
    gauss_points(1:(nx/order)*pt,1:2)=[X1' Y1'];
    weight(1:(nx/order)*pt)=w';
end

if level==3
    pt=4;
zai=[1/3 1/5 3/5 1/5];                    %For degree 3 polynomials
eta=[1/3 1/5 1/5 3/5];
 for i=1:nx/order
     for j=1:length(zai)
     X1_(i,j)=(zai(j)*(x1(i)-cx)+eta(j)*(x1(i+1)-cx)+cx);
     Y1_(i,j)=(zai(j)*(y1(i)-cy)+eta(j)*(y1(i+1)-cy)+cy);
     end
 end
X1=reshape(X1_',[1,nx/order*pt]);
Y1=reshape(Y1_',[1,nx/order*pt]);
nn=1:pt:13;
for i=1:length(nn)
w(nn(i))=-27/96;
end
O=1:nx*pt;
rr=ismember(O,nn);
pp=find(~rr);
for i=1:length(pp)
w(pp(i))=25/96;
end
if (nx/order)==3
    n=0;
end
if (nx/order)==4
    n=1;
end
if (nx/order)==5
    n=2;
end
if (nx/order)==6
    n=3;
end
if (nx/order)==7
    n=4;
end
if (nx/order)==8
    n=5;
end
p=-1:2:pt+(2*n-1);
    gauss_points(1:(nx/order)*pt,1:2)=[X1' Y1'];
    weight(1:(nx/order)*pt)=w';
end

if level==4                                                                                   %For degree 4 polynomials
    pt=6;

zai=[0.091576213509771 0.091576213509771 0.816847572980459  0.445948490915965 0.445948490915965 0.108103018168070];
eta=[0.091576213509771 0.816847572980459 0.091576213509771  0.445948490915965 0.108103018168070 0.445948490915965];
 for i=1:nx/order
     for j=1:length(zai)
     X1_(i,j)=(zai(j)*(x1(i)-cx)+eta(j)*(x1(i+1)-cx)+cx);
     Y1_(i,j)=(zai(j)*(y1(i)-cy)+eta(j)*(y1(i+1)-cy)+cy);
     end
 end
X1=reshape(X1_',[1,nx/order*pt]);
Y1=reshape(Y1_',[1,nx/order*pt]);

mm=1:pt:pt*(nx/order)-(pt-1);

for i=1:nx/order
    jj=mm(i):mm(i)+pt-1;
    w(jj(1:3))=0.109951743655322*0.5;
    w(jj(4:pt))=0.223381589678011*0.5;
end
p=-1:4:4*nx+4;
    gauss_points(1:(nx/order)*pt,1:2)=[X1' Y1'];
    weight(1:(nx/order)*pt)=w';
end

if level==5                                                                                   %For degree 4 polynomials
    pt=7;

zai=[0.333333333333333 0.101286507323456 0.101286507323456  0.797426985353087 0.470142064105115 0.470142064105115 0.059715871789770];
eta=[0.333333333333333 0.101286507323456 0.797426985353087  0.101286507323456 0.470142064105115 0.059715871789770 0.470142064105115];
 for i=1:nx/order
     for j=1:length(zai)
     X1_(i,j)=(zai(j)*(x1(i)-cx)+eta(j)*(x1(i+1)-cx)+cx);
     Y1_(i,j)=(zai(j)*(y1(i)-cy)+eta(j)*(y1(i+1)-cy)+cy);
     end
 end
X1=reshape(X1_',[1,(nx/order)*pt]);
Y1=reshape(Y1_',[1,(nx/order)*pt]);

mm=1:pt:pt*(nx/order)-(pt-1);

for i=1:nx/order
    jj=mm(i):mm(i)+pt-1;
    w(jj(1))=0.225*0.5;
    w(jj(2:4))=0.125939180544827*0.5;
    w(jj(5:pt))=0.132394152788506*0.5;
end
    p=-1:5:7*(nx/order)+1;
end

if level==6                                                                                   %For degree 4 polynomials
    pt=12;
  zai=[0.063089014491502 0.063089014491502 0.873821971016996 0.249286745170910 0.249286745170910 0.501426509658179 0.310352451033785 0.053145049844816 0.310352451033785 0.053145049844816 0.636502499121399 0.636502499121399 ];
  eta=[0.063089014491502 0.873821971016996 0.063089014491502 0.249286745170910 0.501426509658179 0.249286745170910 0.053145049844816 0.310352451033785 0.636502499121399 0.636502499121399 0.053145049844816 0.310352451033785 ];
     for i=1:nx/order
     for j=1:length(zai)
     X1_(i,j)=(zai(j)*(x1(i)-cx)+eta(j)*(x1(i+1)-cx)+cx);
     Y1_(i,j)=(zai(j)*(y1(i)-cy)+eta(j)*(y1(i+1)-cy)+cy);
     end
 end
X1=reshape(X1_',[1,nx/order*pt]);
Y1=reshape(Y1_',[1,nx/order*pt]);
mm=1:pt:pt*(nx/order)-(pt-1);
for i=1:nx/order
    jj=mm(i):mm(i)+pt-1;
    w(jj(1:3))=0.050844906370207*0.5;
    w(jj(4:6))=0.116786275726379*0.5;
    w(jj(7:pt))=0.082851075618374*0.5;
end
    p=-1:10:12*nx+1;
        gauss_points(1:(nx/order)*pt,1:2)=[X1' Y1'];
    weight(1:(nx/order)*pt)=w';
end
if level==7                                                                                   %For degree 4 polynomials
    pt=16;
zai=[0.333333333333333 0.081414823414554  0.459292588292723 0.459292588292723 0.658861384496480 0.170569307751760 0.170569307751760 0.898905543365938 0.050547228317031 0.050547228317031 0.008394777409958  0.008394777409958 0.263112829634638 0.263112829634638 0.728492392955404 0.728492392955404];
eta=[0.333333333333333 0.459292588292723 0.459292588292723 0.081414823414554 0.170569307751760 0.170569307751760 0.658861384496480 0.050547228317031 0.050547228317031 0.898905543365938 0.263112829634638 0.728492392955404 0.008394777409958 0.728492392955404 0.008394777409958 0.263112829634638];     
for i=1:nx/order
     for j=1:length(zai)
     X1_(i,j)=(zai(j)*(x1(i)-cx)+eta(j)*(x1(i+1)-cx)+cx);
     Y1_(i,j)=(zai(j)*(y1(i)-cy)+eta(j)*(y1(i+1)-cy)+cy);
     end
 end
X1=reshape(X1_',[1,nx/order*pt]);
Y1=reshape(Y1_',[1,nx/order*pt]);
mm=1:pt:pt*(nx/order)-(pt-1);
for i=1:nx/order
    jj=mm(i):mm(i)+pt-1;
    w(jj(1))=0.144315607677787*0.5;
    w(jj(2:4))=0.095091634267285*0.5;
    w(jj(5:7))=0.103217370534718*0.5;
    w(jj(8:10))=0.032458497623198*0.5;
    w(jj(11:pt))=0.027230314174435*0.5;
end
    p=-1:(pt-2):pt*nx; 
    gauss_points(1:(nx/order)*pt,1:2)=[X1' Y1'];
    weight(1:(nx/order)*pt)=w';
end
   
end