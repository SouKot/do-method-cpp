[x]=mmread('MeanSolnFinal.mm');
[V]=mmread('V_BASE.mm');
[Y]=mmread('YTrans_COEFF.mm');
xx=linspace(0,1,128);
[ha,hfig]=tight_subplot_cm(2,2,[1.9,1.7],[1.3,0.4],[1.5,.2],18,24, 'yes');
lw=1;
load Exact_data.mat;
exact_sol=val_0_8;
axes(ha(1));
plot(xx', exact_sol.mean,'linewidth',lw);
ylab=ylabel('$\mathbf{E\mathit u}$','interpreter','latex','FontSize',14);
xlab=xlabel('$\mathbf{\mathit x}$','interpreter','latex','FontSize',14)
title('mean');
axes(ha(2));
plot(xx', exact_sol.variance,'linewidth',lw);
ylab=ylabel('$\mathbf{E\mathit u^2}$','interpreter','latex','FontSize',14);
xlab=xlabel('$\mathbf{\mathit x}$','interpreter','latex','FontSize',14);
title('variance');
axes(ha(3));
plot(xx', exact_sol.m3,'linewidth',lw);
ylab=ylabel('$\mathbf{E\mathit u^3}$','interpreter','latex','FontSize',14);
xlab=xlabel('$\mathbf{\mathit x}$','interpreter','latex','FontSize',14);
title('3rd moment');
axes(ha(4));
plot(xx', exact_sol.kurtosis_excess,'linewidth',lw)
ylab2=ylabel('$\mathbf{E\mathit u^4-3}(\mathbf{E} \mathbf{\mathit u^2)^2}$','interpreter','latex','FontSize',14)
%ylab2.Position(1)=ylab2.Position(1)+0.02
xlab=xlabel('$\mathbf{\mathit x}$','interpreter','latex','FontSize',14);
yticks((-2:1:1)*10^-5)
% Apply the labels
%yticklabels(newLabels);
title('excess 4th moment')
hold on



vy=V*Y';
%openfig exact_moments.fig
for i=1:4
  axes(ha(i))
  %subplot(2,2,i)
  if i==1
    hold on;
    plot(xx',x,'linewidth',lw)
    legend('Exact','DO');
   % legend('DO');
  else
    hold on;

    if i==4
      mom = sum(vy.^i,2)/size(Y,1)-3*mom2.^2;
    else
      mom=sum(vy.^i,2)/size(Y,1);
      if i==2
        mom2=mom;
      end
    end
    plot(xx',mom,'linewidth',lw);
  end
end


