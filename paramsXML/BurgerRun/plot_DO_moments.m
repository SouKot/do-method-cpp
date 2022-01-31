[x]=mmread('MeanSolnFinal.mm'); 
[V]=mmread('V_BASE.mm');
[Y]=mmread('YTrans_COEFF.mm');
xx=linspace(0,1,128);
vy=V*Y';
openfig exact_moments.fig
for i=1:4
  i
  subplot(2,2,i)
  if i==1
    hold on;
    plot(xx',x)
    legend('Exact','DO');
  else
    hold on;
    plot(xx',[ sum(vy.^i,2)/size(Y,1) ]);
  end
end
  
    
