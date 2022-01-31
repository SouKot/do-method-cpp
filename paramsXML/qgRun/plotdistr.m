function plotdistr(Y,expect,Var)
n=size(Y,2)
m=size(Y,1)
i=10
f=figure(i+3);
f.Units='centimeters';
f.OuterPosition=[0,0,22,29];
f.PaperPositionMode='auto';
f.PaperType='A4';
p = uipanel('Parent',f,'BorderType','none'); 
%p.Title = [method, sprintf(' ; nx=ny=  %d',nx)]; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';
for j=1:m
  restruct=2;
  if (m==1) restruct=1; end;
  mm=m;
  if (rem(m,2)==1) mm=m+1; end;
  subplot(mm/restruct,restruct,j,'Parent',p);
  histfit(Y(j,:),30)
  title(sprintf('Variance %1.5g', Var(j)))
end
end
