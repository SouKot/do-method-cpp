function [V,expect,Var]=comp_moments(Y, Nmom, excess);
[m,n]=size(Y);
%compute Expectation
fprintf('expectation\n')
expect=sum(Y,2)/n;
%correct
Y=Y-kron(expect,ones(1,n));
% test sum(Y,2)

%Compute variance
[Q,S,V]=svd(Y',0) ;%U S V'
% if 0
%   P=S*V';
% else
%   Y=V'*Y;
%   P=S;
% end
%test: norm(Q*P-Y')
%fprintf('standard deviation\n')
StandDev=diag(S)/sqrt(n);
Var=StandDev.^2;
%fprintf('covariance matrix\n')
%P'*P/n
%rescale to make the standard deviation of Q identity. 
%Q=sqrt(n)*Q;
%P=P/sqrt(n);
%nextmoment=2;
%im=Q;
%stringm=',m'
% for moment=2:nextmoment-1
%   im=kron(im,ones(1,m)).*kron(ones(1,size(im,2)),Q);
%   stringm=[stringm,',m']
% end
% for moment=nextmoment:Nmom
%   im=kron(im,ones(1,m)).*kron(ones(1,size(im,2)),Q);
%   stringm=[stringm,',m']
%   T=sum(im)/n;
%   fprintf('Q moment %d\n',moment)
%   T=eval(['reshape(T',stringm,')']);
%   if excess
%     if moment==2
%       for i=1:m
% 	T(i,i)=T(i,i)-1;
%       end
%     end
%     if moment==4
%       for i=1:m
% 	for j=1:m
% 	  T(i,i,j,j)=T(i,i,j,j)-1;
% 	  T(i,j,i,j)=T(i,j,i,j)-1;
% 	  T(i,j,j,i)=T(i,j,j,i)-1;
% 	end
%       end
%     end
%   end
%   for l=1:moment
%     T=reshape(T,m^(moment-1),m);
%     T=T*P;
%     T=T';
%   end
%   fprintf('Moment %d\n',moment)
%   eval(['reshape(T',stringm,')'])
%   eval(['moments.m',int2str(moment),'=T;']);
end
