function []=analyzeNplot(mnflnm,vflnm,coeflnm,nx,ny,Re,stochsize,vcol,flnm)

addpath /home/sourabh/fredwubs/MatlabScripts
%vcol=1; %1:4
%  fid=fopen('timeslicesDt1e-3ndtsub10m4Frc1.txt')
%  fid=fopen('timeslicesDt1e-3ndtsub10m4.txt')
%  fid=fopen('timeslicesSymForcing_dt1e-3_ndtsub10_m4.txt')
%timeslicesDt1e-3ndtsub10m4Frc10.txt')
  cnt=0;
ts=[];

  close all
  cnt=cnt+1;
  %[V,Y,x,t,magnitudes]=rdtimeslice(fid,stochsize);
  %[~,~,x]=textread([dir,'MeanSolnFinal.txt'],'%f %f %f','headerlines',1); 
  [V]=mmread(vflnm);
  [x]=mmread(mnflnm);
  [Y]=mmread([coeflnm]);
  % dir= 're_20_m_1_sfc_0_T_0.5_stochiter_10e4/'
  % [V]=mmread([dir,'v_4.79999989e-01.mm']);
  % [x]=mmread([dir,'mean4.79999989e-01.mm']);
  % [Y]=mmread([dir,'yT_4.79999989e-01.mm']);
  
  Y=Y';
  size(x);
    [n,m]=size(V);
  size(Y);
  %Y=Y';
  %pause
  excess=0; 
  %Vt=comp_moments(Y,4,excess);
  %excess=1
  [Vt,expect,Var]=comp_moments(Y,4,excess);
  Y=Vt'*Y;
  V=V*Vt;
  m=size(vcol,2);
%   DS=diag(Var);
%    
%   D=spdiags(kron([-1 2 -1],ones(nx,1)),[-1,0,1],nx,nx);
%   %Since psi has zeroes all around the first and last column/row should be made zero 
%   D([1,nx],:)=0;D(:,[1,nx])=0;
%   D=kron(speye(nx),D)+kron(D,speye(nx)); D=kron(D,sparse([0 0;0 1]));
%   [V,R]=qrM(V,D);
%   [P,Vd]=eig(R*DS*R');
%   [Var,iVar]=sort(diag(Vd),'descend');
%   fprintf('Variances for Laplacian-orthogonal basis\n')
%   Var'
%   %(R*DS*R')P=P*Vd -> V*R*DS*R'*V'= V*P*Vd*P'V' 
%   V=V*P(:,iVar);
%   m=4
%   Y=P(:,iVar)'*(R*Y);
% Var
% size(V)
% [x,V(:,vcol)];
% [0;Var(vcol)];
%if vcol~=0
  showV([V(:,vcol)],[Var(vcol)],flnm,1,nx,ny,2,1,m)
%else
  showV([x],[0],flnm,1,nx,ny,2,2,1)
%end
  %pause
  %indYp=find(Y(1,:) > 1.0); length(indYp);
  %indYn=find(Y(1,:) < -1.0); length(indYn);
  %Yp=Y(:,indYp); expectp=sum(Yp,2)/length(indYp);  
  %Yn=Y(:,indYn); expectn=sum(Yn,2)/length(indYn); 
 % showV([x+expectp(1)*V(:,1),x+expectn(1)*V(:,1)],[0;0],'Approx. Stable solutions',1,nx,nx,2,1,2)
  %showV(x,1,'Average sol',1,nx,nx,2,1,1)
  %pause
  %Y=P(:,iVar)'*(R*Y);
  
  %excess=0; 
  %Vt=comp_moments(Y,4,excess);
  %excess=1;
  %[Vt,expect,moments,Var]=comp_moments(Y,4,excess);
  %Y=Vt'*Y;
  %V=V*Vt;
  %Vt % should be identity
  
  %Var;
  %CoSkewn=moments.m3; 
  %CoKurt=1; 
  %CoExKurt=moments.m4;
  mV=size(V,2);
  %[Etdis(1:mV,cnt),EtmeanV(1:mV,cnt)]=energies(V,x,Var,CoSkewn,Re); 
  %ts=[ts,t];
  plotdistr(Y,expect,Var);

  mnsf=reshape(x(2:2:end),nx,ny);
  x=transpose(fliplr(mnsf));
  size(x); 
  [Fy,Fx]=gradient(x,1/(nx-1));
  disp([min(min(Fy))]);
  disp([max(max(Fy))]);
  disp(["minimum of v = ",num2str(min(min(Fx)))]);
  disp(["maximum of v = ",num2str(max(max(Fx)))]);




