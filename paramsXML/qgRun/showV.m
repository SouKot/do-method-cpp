function []=showV(V,var,method,nz,nx,ny,blksize,CASE,mp)
% []=showV(V,var,method,nz,nx,ny,blksize,CASE,mp);
% V  : nx*nz*ny*blksize x m matrix;
% var: m variances for printing above the plot
% method : string to be added on top of the plot
% nz : number of layers
% nx : number of points in x-direction
% ny : number of points in y-direction
% blksize : number of unknowns per cell/gridpoint
% CASE : indicates which crosscut plane should be visualized
% mp: number of vectors in V to be visualized.
%subplot = @(m,n,p) subtightplot (m, n, p, [0.1 -0.05], [0.06 0.06], [0.04 0.04]);
tt1=["modeVorticity_"+method,"modeStream_"+method];
tt2=["meanVorticity_"+method,"meanStream_"+method];
exp1=[-2,-4];
n=size(V,1);
m=size(V,2);
n=nx*nz*ny*blksize;
for i=1:blksize
  f=figure(i+3*CASE);
  %f.Units='centimeters';
  %f.OuterPosition=[0,0,22,29];
  %f.PaperPositionMode='auto';
  %f.PaperType='A5';
  %p = uipanel('Parent',f,'BorderType','none'); 
  %p.Title = [method, sprintf(' ; nx=ny=  %d',nx)]; 
  %p.TitlePosition = 'centertop'; 
  p.FontSize = 12;
  p.FontWeight = 'bold';
  restruct=2; 
  gc=custom_subplot_gap(ceil(mp/restruct),restruct,0.01, 0.05)
  for j=1:mp
    fld=reshape(V(i:blksize:n,j),nx,ny,nz);
    switch CASE 
    case 1
      for k=1:nz    
        %subplot(ceil(mp/restruct),restruct,k+(j-1)*nz);
        %subaxis(ceil(mp/restruct),restruct,k+(j-1)*nz,'sh',-0.1,'sv',0.0,'mb',0.0,'mt',0.0,'ml',0.0,'mr',0.0,'pl',0.0,'pr',0.0,'pt',0.0,'pb',0.0)
	%gc=smplot(ceil(mp/restruct),restruct,k+(j-1)*nz,'left',0,'right',0)
        axes(gc(j));
        contourf(transpose(fliplr(fld(:,:,k))),30);
        t=title(['variance = ',num2str(var(j),'%.2f')]);
        t.FontSize=11;
        set(gc(j),'Yticklabel',[]) 
        set(gc(j),'Xticklabel',[])
        %invl=ceil(length(h.LevelList)/5);
        %clabel(C,h,h.LevelList(1:invl:end));
        axis equal
        savefig([tt1(i)+'.fig']);
        print('-painters','-dpdf',[tt1(i)]);
	%holdaxis
           %cb = colorbar();
        %cb.Ruler.Exponent = exp1(i);
      end
    case 2
      restruct=1; 
      for k=1:nz    
        subplot(ceil(mp/restruct),restruct,k+(j-1)*nz);
        %colormap(redblue(45))
        %contourf(transpose(fliplr(fld(:,:,k))),80,'edgecolor','none')
        contourf(transpose(fliplr(fld(:,:,k))),30)
        set(gca,'Yticklabel',[]) 
        set(gca,'Xticklabel',[])
        axis equal
        savefig([tt2(i)+'.fig']);
        print('-painters','-dpdf',[tt2(i)]);
 %colorbar
      end
    case 3
      nxl=1;
      for k=1:nxl
        subplot(mp,nxl,k+(j-1)*nxl);
        fldtmp(:,:)=fld(k,:,:);
        %colormap(redblue(80))
	      %contourf(fldtmp',80,'edgecolor','none')
        contourf(fldtmp',30)
      end
    case 4
      for k=1:ny
        subplot(mp,ny,k+(j-1)*ny);
        fldtmp(:,:)=fld(:,k,:);
        %colormap(redblue(80))
        %contourf(fldtmp,80,'edgecolor','none')
        contourf(fldtmp,30)
      end
    end
    %if CASE==1
      %savefig([tt1(i)+'.fig']);
      %print('-painters','-dpdf',tt1(i));
    %end
  end
end
