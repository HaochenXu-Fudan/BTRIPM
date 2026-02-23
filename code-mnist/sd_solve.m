function grad=sd_solve(diff,epsilon,x,y,J)
n=length(y);
x0=zeros(n,1);
iter=0;
x=x0;
r=diff.fxy(x,y,x0)-diff.fyyv(x,y,x0);
rr=dot(r,r);
while sqrt(rr)>=epsilon&&iter<J
    iter=iter+1;
    Ar=diff.fyy(x,t,r);
    alpha=rr/dot(r,Ar);
    x=x+alpha*r;
    r=r-alpha*Ar;
    rr=dot(r,r);
end
end

function q = cg_solver(diff,epsilon,x,y,J)
n=length(lbd);
q0=zeros(n,1);
i=0;
r=diff.fxy(x,y,q0)-diff.fyyv(x,y,q0);
delta_new=r'*r;
delta0=delta_new;


end


function[x,obj]=cg_solver(A,b,x0,epsilon,max_iter)
x=x0; %初始点
i=0;
r=b-A*x; % 计算残差 
d=r; %方向
delta_new=r'*r; %残差的内积
delta0=delta_new;
obj=0.5*x'*A*x-b'*x; %目标函数
while(i<max_iter)&&(delta_new/delta0>epsilon^2)
    q=A*d;
    alpha=delta_new/(d'*q); %计算系数 alpha
    x=x+alpha*d; %迭代一步（向下走一步）
    obj=[obj,0.5*x'*A*x-b'*x];
    r=b-A*x; % 更新计算残差
    delta_old=delta_new;
    delta_new=r'*r; %更新残差的内积
    beta=delta_new/delta_old; %更新下一步方向的系数
    d=r+beta*d; %更新下一步的方向
    i=i+1;
end
end



