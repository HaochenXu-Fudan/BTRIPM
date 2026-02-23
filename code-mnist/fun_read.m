function [F,f,diff,others] = fun_read(funname)
% p:features, c:categories, 
switch funname
    case{'20news'}      
%         Xtr Ytr Xval Yval sparse, n, p, c
        load('Xtr.mat');load('Xval.mat');load('Ytr.mat');load('Yval.mat');
        [nval,p] = size(Xval);[ntr,~] = size(Xtr);c = 20;
        diff = struct();
        F = @(x,y) Fval(x,y,Xval,Yval,nval,p,c);
        f = @(x,y) fval(x,y,Xtr,Ytr,ntr,p,c);
        diff.Fx = @(x,y) zeros(p*c,1);
        diff.Fy = @(x,y) diffFy(x,y,Xval,Yval,nval,p,c);
        diff.fx = @(x,y) exp(x).*(y.^2)/(2*p*c);
        diff.fy = @(x,y) difffy(x,y,Xtr,Ytr,ntr,p,c);
        diff.Fxx = @(x,y,v) zeros(p*c,1);
        diff.Fxy = @(x,y,v) zeros(p*c,1);
        diff.Fyx = @(x,y,v) zeros((p)*c,1);
        diff.Fyy = @(x,y,v) diffFyy(x,y,Xval,Yval,nval,p,c,v);
        diff.fxx = @(x,y,v) (exp(x).*(y.^2)).*v/(2*p*c);
        diff.fxy = @(x,y,v) (exp(x).*y).*v/(p*c);
        diff.fyx = @(x,y,v) (exp(x).*y).*v/(p*c);
%         diff.fyy = @(x,y) difffyy_CE(x,y,input_tr,m2,n2,gtr,gtr_nonzero)*const;
        diff.fyyv = @(x,y,v) difffyy_v(x,y,Xtr,Ytr,ntr,p,c,v);
        others = struct();
%     case{'20newssgd'}
% %         Xtr Ytr Xval Yval sparse, n, p, c
%         load('Xtr.mat');load('Xval.mat');load('Ytr.mat');load('Yval.mat');
%         [nval,p] = size(Xval);[ntr,~] = size(Xtr);c = 20;
%         
%         diff = struct();
%         F = @(x,y) Fval(x,y,Xval,Yval,nval,p,c);
%         f = @(x,y) fval(x,y,Xtr,Ytr,ntr,p,c);
%         diff.Fx = @(x,y) zeros(p*c,1);
%         diff.Fy = @(x,y) diffFy(x,y,Xval,Yval,nval,p,c);
%         diff.fx = @(x,y) exp(x).*(y.^2)/(2*p*c);
%         diff.fy = @(x,y) difffy(x,y,Xtr,Ytr,ntr,p,c);
%         diff.Fxx = @(x,y,v) zeros(p*c,1);
%         diff.Fxy = @(x,y,v) zeros(p*c,1);
%         diff.Fyx = @(x,y,v) zeros((p)*c,1);
%         diff.Fyy = @(x,y,v) diffFyy(x,y,Xval,Yval,nval,p,c,v);
%         diff.fxx = @(x,y,v) (exp(x).*(y.^2)).*v/(2*p*c);
%         diff.fxy = @(x,y,v) (exp(x).*y).*v/(p*c);
%         diff.fyx = @(x,y,v) (exp(x).*y).*v/(p*c);
% %         diff.fyy = @(x,y) difffyy_CE(x,y,input_tr,m2,n2,gtr,gtr_nonzero)*const;
%         diff.fyyv = @(x,y,v) difffyy_v(x,y,Xtr,Ytr,ntr,p,c,v);
%         others = struct();
end
end
function z = Fval(x,y,Xval,Yval,n,p,c)
W = reshape(y,p,c);
P = getP(W,Xval);
z = -ones(1,n)*((Yval.*log(P))*ones(c,1))/n;
end
function z = fval(x,y,Xtr,Ytr,n,p,c)
W = reshape(y,p,c);
P = getP(W,Xtr);
z = -ones(1,n)*((Ytr.*log(P))*ones(c,1))/n;
z = z +exp(x)'*(y.^2)/(2*c*p);
end
function z = diffFy(x,y,Xval,Yval,n,p,c)
W = reshape(y,p,c);
P = getP(W,Xval);
temp = Xval'*(P-Yval)/n;
z = temp(:);
end
function z = difffy(x,y,Xtr,Ytr,n,p,c)
W = reshape(y,p,c);
P = getP(W,Xtr);
temp = Xtr'*(P-Ytr)/n;
z = temp(:)+exp(x).*y/(c*p);
end
function z = diffFyy(x,y,Xval,Yval,n,p,c,v)
W = reshape(y,p,c);
P = getP(W,Xval);
v = reshape(v,p,c);
U1 = P.*(Xval*v);
U2 = Xval'*U1;
V1 = sum(U1,2);
V2 = (V1*ones(1,c)).*P;
V3 = (Xval'*V2);
z = U2(:)-V3(:);
z = z/n;
end
function z = difffyy_v(x,y,Xtr,Ytr,n,p,c,v)
W = reshape(y,p,c);
P = getP(W,Xtr);
v = reshape(v,p,c);
U1 = P.*(Xtr*v);
U2 = Xtr'*U1;
V1 = sum(U1,2);
V2 = (V1*ones(1,c)).*P;
V3 = (Xtr'*V2);
z = U2(:)-V3(:);
z = z/n;
end
% W:(p)*c, X:n*(p)
function y = getP(W,X)
c = 20;
XW = X*W;
XW = XW-max(XW,2);
expXW = exp(XW)+1e-12;
sumexpXW = sum(expXW,2);
y = expXW./(sumexpXW*ones(1,c))+1e-13;
end