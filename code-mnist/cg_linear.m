function grad_pre = cg_linear(diff,epsilon,x,y,J)
%n=length(A);
n=length(y);
% if(nargin<5);x0=zeros(n,1);end
% if(nargin<4);epsilon=10^(-3);end
x0=zeros(n,1);
r0=diff.Fy(x,y)-diff.fyyv(x,y,x0);
r=r0;
d=r0;
m=zeros(n,1);
p=zeros(n,1);
grad=x0;
% min(n-1,200) = J
for k=0:J
    m(k+1)=r'*r/(d'*diff.fyyv(x,y,d));
    grad=grad+m(k+1)*d;
    r=diff.Fy(x,y)-diff.fyyv(x,y,grad);
    if(norm(r,inf)<=epsilon||k+1==n)
        break;
    end
    p(k+1)=norm(r)^2/norm(r0)^2;
    d=r+p(k+1)*d;
    r0=r;
end
grad_pre=grad;
end





