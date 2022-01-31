function []=read_timeseries_Eyy_N_plot(dirnm,nV)
close all;
%nV=4; 
creationdate='10-Nov';
timecorrection=10;
cd(dirnm);
file=[dir('tsEyy_*e-01.mm');dir('tsEyy_*e+00.mm');dir('tsEyy_*e+01.mm')];
nm=extractfield(file, 'name');
Date=extractfield(file,'date')
%indx=find(contains(Date,'14-Nov'));
%storeTime=[0,str2double(extractBetween(nm(indx),'_','.mm'))];
storeTime=[0,str2double(extractBetween(nm,'_','.mm'))];
color=distinguishable_colors(nV);
intime=0;
flnm='Eyy_eigs_convergence';
figure(20)
hold on;
for i=1:length(storeTime)-1
	nm(i);
	[mat,rw,cl]=mmread(char(nm(i)));
	eigvals=zeros(nV,cl);
	dt=(storeTime(i+1)-storeTime(i))/cl;
	t=linspace(intime+dt,storeTime(i+1),cl);
	[t(1),t(end)]
	for j=1:cl
		cmat=reshape(mat(:,j),nV,nV);	
		eigvals(:,j)=eig(cmat);	
	end
	%eigvals(:,50:60)
	%pause()	
	for k=1:nV	
		%size(t)
		%size(eigvals(k,:))
		%scatter(t,log10(eigvals(k,:)),[],color(k,:))
		semilogy(t,(eigvals(k,:)),'color',color(k,:),'LineWidth',1.5)
		%plot(t,log10(eigvals(k,:)),'color',color(k,:),'LineWidth',1.5)

	end
	intime=t(end);
	pause(0.2);
end
savefig([flnm,'.fig']);
print('-painters','-dpdf',flnm);
cd ..
