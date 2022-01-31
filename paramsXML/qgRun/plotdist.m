
clear all;
nV=4;
timecorrection=10;
%cd(dirnm);
file=[dir('yT_*e-01.mm');dir('yT_*e+00.mm');dir('yT_*e+01.mm')];
%file=[dir('yT_1.7*e+01.mm')];
%file=[dir('yT_*e+00.mm')];
nm=extractfield(file, 'name');
Date=extractfield(file,'date')
%indx=find(contains(Date,'14-Nov'));
%storeTime=[0,str2double(extractBetween(nm(indx),'_','.mm'))];
storeTime=[0,str2double(extractBetween(nm,'_','.mm'))];
color=distinguishable_colors(nV);
intime=0;
%flnm='Eyy_eigs_convergence';
%figure(20)
hold on;
k=0;
for i=1:length(nm)
	nm(i);
	y=mmread(char(nm(i)));
  m=size(y,2);
  k=k+1;
  for j=1:1
    figure(j)
    subplot(5,5,k)
    %restruct=2;
    %subplot(mm/restruct,restruct,j,'Parent',p);
    histfit(y(:,j),30)
  end
  if k==25
    pause
    close all
    k=0;
  end  
end


% function plotdist(Y)
% n=size(Y,1)
% m=size(Y,2)
% i=10
% f=figure(i+3);
% f.Units='centimeters';
% f.OuterPosition=[0,0,22,29];
% f.PaperPositionMode='auto';
% f.PaperType='A4';
% p = uipanel('Parent',f,'BorderType','none'); 
% %p.Title = [method, sprintf(' ; nx=ny=  %d',nx)]; 
% p.TitlePosition = 'centertop'; 
% p.FontSize = 12;
% p.FontWeight = 'bold';
% for j=1:m
%   restruct=2;
%   if (m==1) 
%     restruct=1; 
%   end
%   mm=m;
%   if (rem(m,2)==1) 
%     mm=m+1; 
%   end
%   subplot(mm/restruct,restruct,j,'Parent',p);
%   histfit(Y(:,j),30)
%   %title(sprintf('Variance %))
% end
% end
