clear
clc
%1.改气体名字，co2\h2o\o2;改absco的index
%2.改find(wave=)与对应的后两句wave范围
%3.改气体model
h2o_surface=10^4;
ke=1.380649*10^-23;
co2file='F:\Group-rt\OCO-2data\table\5.1\co2_v51.hdf';
h2ofile='F:\Group-rt\OCO-2data\table\5.1\h2o_v51.hdf';
absco=h5read(co2file,'/Gas_02_Absorption');
VMR=h5read(co2file,'/Broadener_01_VMR');
pressure=h5read(co2file,'/Pressure');
temperature=h5read(co2file,'/Temperature');
wave=h5read(co2file,'/Wavenumber');
p=load('1\labels_P.dat');
p=p*100;
t=load('1\labels_T.dat');
x=load('1\labels_x.dat');
x=x/1000000;
dx=load('1\labels_dx.dat');
h2o=h2o_surface*10.^(-0.2*dx);
%h2o=load('1\labels_x2.dat');
h2o=h2o/1000000;

[k1,k2]=find(wave==6180);%range-absco-wave
[k3,k4]=find(wave==6280);%range-absco-wave
absco=absco(k1:k3,:,:,:);%wave,x,t,p
wave=wave(k1:k3,1);

%————————rand version—————————
%x_grid=[200,600];
% x_grid=VMR;
% for i=1:sample
%     [p,t,x,t_grid1,t_grid2]=randptx(pressure,temperature,x_grid,pmodel);
%     outabsco(i,:)=inter4d(x,t,p,absco,VMR,t_grid1,t_grid2,pressure);
% %    kn=10^(-6)*p*x*10^(-6)/t/ke;%linear-based
% %    outabsco(i,:)=outabsco(i,:).*kn;
%     outprofile(i,:)=[p t x];
% end
% [x1,x2]=size(outprofile);
% [y1,y2]=size(outabsco);
% 
% C1={'%4e,';'%4e\n'};
% 
% for i=1:sample
%     fprintf(fid1,[C1{[ones(1,y2-1) 2]}],outabsco(i,:));
% end
% fclose(fid1);
% 
% C2={'%2e,';'%2e\n'};
% 
% for i=1:sample
%     fprintf(fid2,[C2{[ones(1,x2-1) 2]}],outprofile(i,:));
% end
% fclose(fid2);



% dlmwrite('co2absco_test.dat',outabsco,'precision','%.6e')
% dlmwrite('co2profile_test.dat',outprofile,'precision','%.3f')
%————————cal version—————————
for i=1:28
    vh2o=h2o(i)/(1-h2o(i));
%x=0.03;t=250;p(i)=pressure(i);
[k1,k2]=find((pressure-p(i))>0);
k1=k1(1);
t_grid1=temperature(:,k1-1);
t_grid2=temperature(:,k1);
outabsco(i,:)=inter4d(vh2o,t(i),p(i),absco,VMR,t_grid1,t_grid2,pressure);
kn=10^(-6)*p(i)*x(i)/t(i)/ke;%linear-based
outabsco(i,:)=outabsco(i,:).*kn;
end

outabsco_co2=outabsco;




co2file='F:\Group-rt\OCO-2data\table\5.1\h2o_v51.hdf';
absco=h5read(co2file,'/Gas_01_Absorption');
VMR=h5read(co2file,'/Broadener_01_VMR');
pressure=h5read(co2file,'/Pressure');
temperature=h5read(co2file,'/Temperature');
wave=h5read(co2file,'/Wavenumber');


fid=fopen('1\labels_absco_all.dat','w');
fid1=fopen('1\labels_absco_co2.dat','w');
fid2=fopen('1\labels_absco_h2o.dat','w');
[k1,k2]=find(wave==6180);%range-absco-wave
[k3,k4]=find(wave==6280);%range-absco-wave
absco=absco(k1:k3,:,:,:);%wave,x,t,p
wave=wave(k1:k3,1);


for i=1:28
    vh2o=h2o(i)/(1-h2o(i));
%x=0.03;t=250;p(i)=pressure(i);
[k1,k2]=find((pressure-p(i))>0);
k1=k1(1);
t_grid1=temperature(:,k1-1);
t_grid2=temperature(:,k1);
outabsco(i,:)=inter4d(vh2o,t(i),p(i),absco,VMR,t_grid1,t_grid2,pressure);
kn=10^(-6)*p(i)*h2o(i)/t(i)/ke;%linear-based
outabsco(i,:)=outabsco(i,:).*kn;
end
outabsco_h2o=outabsco;
outabsco=outabsco_h2o+outabsco_co2;
[y1,y2]=size(outabsco);
C1={'%4e,';'%4e\n'};
for i=1:28
    fprintf(fid,[C1{[ones(1,y2-1) 2]}],outabsco(i,:));
    fprintf(fid1,[C1{[ones(1,y2-1) 2]}],outabsco_co2(i,:));
    fprintf(fid2,[C1{[ones(1,y2-1) 2]}],outabsco_h2o(i,:));
end

% figure(1)
% set(gca,'FontName','Times New Roman' ,'FontSize',18);
% for i=1:15
%     plot(wave,outabsco(i,:),'linewidth',0.8)
%     hold on
%     legend_str{i}=['P=' num2str(p(i))];
% end
% legend(legend_str)
% xlabel('wavenumber[cm^{-1}]')
% ylabel('Absorption cross-section[cm^{2}/molec]')
% kn=10^(-6)*p*x*10^(-6)/t/ke;%linear-based
% outabsco(1,:)=outabsco(1,:).*kn;
% outprofile(1,:)=[p t x];
% out(:,1)=wave;out(:,2)=outabsco';
% fname = sprintf('p%.1ft%.2fx%.1f.dat', p,t,x);
% save(fname,'out','-ascii')
%———————————————————————
