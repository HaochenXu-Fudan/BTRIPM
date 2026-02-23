function [F,f,diff,others] = fun_read(funname)
switch funname
    case{'toy0'}
        others = struct();
        F = @(x,y) x^2+y^2;
        f = @(x,y) sin(x+y);
        diff = struct();
        diff.Fx = @(x,y) 2*x;
        diff.Fy = @(x,y) 2*y;
        diff.Fxx = @(x,y,v) 2*v;
        diff.Fyy = @(x,y,v) 2*v;
        diff.Fxy = @(x,y,v) 0;
        diff.Fyx = @(x,y,v) 0;
        diff.fx = @(x,y) cos(x+y);
        diff.fy = @(x,y) cos(x+y);
        diff.fxx = @(x,y,v) -sin(x+y)*v;
        diff.fyy = @(x,y) -sin(x+y);
        diff.fyyv = @(x,y,v) -sin(x+y)*v;
        diff.fxy = @(x,y,v) -sin(x+y)*v;
        diff.fyx = @(x,y,v) -sin(x+y)*v;
    case{'toy2'}
        others = struct();
        F = @(x,y) (x-2)^2+(y-2)^2;
        f = @(x,y) sin(x+y);
        diff = struct();
        diff.Fx = @(x,y) 2*(x-2);
        diff.Fy = @(x,y) 2*(y-2);
        diff.Fxx = @(x,y,v) 2*v;
        diff.Fyy = @(x,y,v) 2*v;
        diff.Fxy = @(x,y,v) 0;
        diff.Fyx = @(x,y,v) 0;
        diff.fx = @(x,y) cos(x+y);
        diff.fy = @(x,y) cos(x+y);
        diff.fxx = @(x,y,v) -sin(x+y)*v;
        diff.fyy = @(x,y) -sin(x+y);
        diff.fyyv = @(x,y,v) -sin(x+y)*v;
        diff.fxy = @(x,y,v) -sin(x+y)*v;
        diff.fyx = @(x,y,v) -sin(x+y)*v;
    case{'wifi'}
        others = struct();
%         数据导入和预处理
        data = importdata('data/wifi.txt');
        index = (data(:,8)<3);
        data(index,8) = 1;
        data(~index,8) = -1;
        [m,~] = size(data);
        data(:,1:7) = data(:,1:7)./(ones(m,1)*max(data(:,1:7) ));
%         划分训练集和测试集
%         CVO = cvpartition(m,'Kfold',2);
%         trainIndices = CVO.training(1);
%         testIndices = CVO.test(1);
        CVO = cvpartition(m,'holdout',0.5);
        trainIndices = CVO.training;
        testIndices = CVO.test;
%         数据第八列为输出, 下区分输入和输出。
        output_tr = data(trainIndices,8);
        input_tr = data(trainIndices,:);
        input_tr(:,8) = [];
        [m1,n] = size(input_tr);
        output_te = data(testIndices,8);
        input_te = data(testIndices,:);
        input_te(:,8) = [];
        [m2,~] = size(input_te);
%         上层变量写为x,下层变量写为y
        const = 1e-3;
        F = @(x,y) F_logistic(x,y,input_te,output_te,m2)*const;
        f = @(x,y) f_logistic(x,y,input_tr,output_tr,m1)*const;
        diff.Fx = @(x,y) 0;
        diff.Fy = @(x,y) diffFy_logistic(x,y,input_te,output_te,m2,n)*const;
        diff.Fxx = @(x,y,v) 0;
        diff.Fyy = @(x,y,v) diffFyy_logistic(x,y,v,input_te,output_te,m2)*const;
        diff.Fxy = @(x,y,v) 0;
        diff.Fyx = @(x,y,v) zeros(7,1);
        diff.fx = @(x,y) exp(x)*norm(y)^2*const;
        diff.fy = @(x,y) difffy_logistic(x,y,input_tr,output_tr,m1,n)*const;
        diff.fxx = @(x,y,v) exp(x)*norm(y)^2*const*v;
        diff.fyy = @(x,y) difffyy_logistic(x,y,input_tr,output_tr,m1,n)*const;
        diff.fyyv = @(x,y,v) difffyyv_logistic(x,y,v,input_tr,output_tr,m1,n)*const;
        diff.fxy = @(x,y,v) 2*exp(x)*(y'*v)*const;
        diff.fyx = @(x,y,v) 2*exp(x)*(y*v)*const;
    case{'crowdsourced'}
        others = struct();
%         data_1 = xlsread('data/crowdsourced_tr.xlsx');
%         data_2 = xlsread('data/crowdsourced_te.xlsx');
        data_1 = readmatrix('data/crowdsourced_tr.xlsx');
        data_2 = readmatrix('data/crowdsourced_te.xlsx');
        data = [data_1;data_2];
        [m,~] = size(data);
        R = randperm(m);
        data_tr = data(R(1:5000),:);
        data_te = data(R(5001:10000),:);
        
        input_tr = data_tr(:,2:29);
        output_tr = data_tr(:,1);
        input_te = data_te(:,2:29);
        output_te = data_te(:,1);
        
        [m1,n] = size(input_tr);
        [m2,~] = size(input_te);
        regular = max([input_tr;input_te]);
        input_tr = input_tr./(ones(m1,1)*max(regular));
        input_te = input_te./(ones(m2,1)*max(regular)); 
        
%         const = 1e-5;
        const = 1e-4;
         F = @(x,y) F_logistic(x,y,input_te,output_te,m2)*const;
        f = @(x,y) f_logistic(x,y,input_tr,output_tr,m1)*const;
        diff.Fx = @(x,y) 0;
        diff.Fy = @(x,y) diffFy_logistic(x,y,input_te,output_te,m2,n)*const;
        diff.Fxx = @(x,y,v) 0;
        diff.Fyy = @(x,y,v) diffFyy_logistic(x,y,v,input_te,output_te,m2)*const;
        diff.Fxy = @(x,y,v) 0;
        diff.Fyx = @(x,y,v) zeros(28,1);
        diff.fx = @(x,y) exp(x)*norm(y)^2*const;
        diff.fy = @(x,y) difffy_logistic(x,y,input_tr,output_tr,m1,n)*const;
        diff.fxx = @(x,y,v) exp(x)*norm(y)^2*const*v;
        diff.fyy = @(x,y) difffyy_logistic(x,y,input_tr,output_tr,m1,n)*const;
        diff.fyyv = @(x,y,v) difffyyv_logistic(x,y,v,input_tr,output_tr,m1,n)*const;
        diff.fxy = @(x,y,v) 2*exp(x)*(y'*v)*const;
        diff.fyx = @(x,y,v) 2*exp(x)*(y*v)*const;
    case{'mnist5000'}
        others = struct();
        mnist_image = loadMNISTImages('data/mnist_train_image.idx3-ubyte');
        mnist_label = loadMNISTLabels('data/mnist_train_label.idx1-ubyte');
        R = randperm(60000);
        input_tr = mnist_image(:,R(1:5000));
        spinput_tr = sparse(input_tr);
        [n2,m2] = size(input_tr);
%         m2: data n2: features
        output_tr = mnist_label(R(1:5000),1);
%         contaminate 10% of the labels
        Q = randperm(5000);
        q = 2500; % number of contaminated samples
%       larger the q, lower trust on lower problem, which leads to increase
%       of the lower loss.
        c = Q(1:q);
        others.c = c;
        for i=1:q
            a = rand;
            output_tr(c(i)) = mod(output_tr(c(i))+ceil(9*a),10);
        end
%         Noted! here test set is validation set in Liu.et al 2021
        input_te = mnist_image(:,R(5001:10000));
        output_te = mnist_label(R(5001:10000),1);
        % for reducing computation
        spinput_te = sparse(input_te);
        [n1,m1] = size(input_te);
        
%         classify the labels
        Btr = classify(output_tr);
        Btr =(Btr==1);
        
        Bte = classify(output_te);
%         bte = sum(Bte);
        Bte =(Bte==1);
        gtr = [repmat(input_tr,10,1);ones(10,m2)];
        gte = [repmat(input_te,10,1);ones(10,m1)];
%         spgte = sparse(gte);spgtr = sparse(gtr);
        sumyI_te = getsumyI(input_te,Bte,m1,n1);
        gtr_nonzero = (gtr~=0);
        
        const = 5e-4;
        F = @(x,y) F_CEloss(x,y,input_te,Bte,m1,n1)*const;
        f = @(x,y) f_CEloss(x,y,input_tr,Btr,m2,n2)*const;
%         others.getp_te = @(W) getp(W,input_te,m1,n1);
%         others.getp_tr = @(W) getp(W,input_tr,m2,n2);
        diff.Fx = @(x,y) zeros(m2,1);
        diff.Fy = @(x,y) diffFy_CE(x,y,input_te,m1,n1,sumyI_te)*const;
        diff.fy = @(x,y) difffy_CE(x,y,input_tr,Btr,m2,n2,gtr)*const;
        diff.Fxx = @(x,y,v) zeros(m2,1);
        diff.Fxy = @(x,y,v) zeros(m2,1);
        diff.Fyx = @(x,y,v) zeros(10*n1+10,1);
        diff.Fyy = @(x,y,v) diffFyy_CE(x,y,v,input_te,m1,n1,spinput_te)*const;
        diff.fx = @(x,y) difffx_CE(x,y,input_tr,Btr,m2,n2)*const;
        diff.fxx = @(x,y,v) diff.fx(x,y).*ddsigmoid(x).*v*const;
        diff.fxy = @(x,y,v) difffxy_CE(x,y,v,input_tr,Btr,m2,n2,spinput_tr)*const;
        diff.fyx = @(x,y,v) difffyx_CE(x,y,v,input_tr,Btr,m2,n2,spinput_tr)*const;
        diff.fyy = @(x,y) difffyy_CE(x,y,input_tr,m2,n2,gtr,gtr_nonzero)*const;
        diff.fyyv = @(x,y,v) difffyyv_CE(x,y,v,input_tr,m2,n2,spinput_tr)*const;
        case{'fashion_mnist'}
        others = struct();
        fmnist_image = loadMNISTImages('data/fashion_train-images-idx3-ubyte');
        fmnist_label = loadMNISTLabels('data/fashion_train-labels-idx1-ubyte');
        R = randperm(60000);
        input_tr = fmnist_image(:,R(1:5000));
        [n2,m2] = size(input_tr);
%         m2: data n2: features
        output_tr = fmnist_label(R(1:5000),1);
        % for reducing computation
        spinput_tr = sparse(input_tr);
%         contaminate 10% of the labels
        Q = randperm(5000);
        q = 2500; % number of contaminated samples
%       larger the q, lower trust on lower problem, which leads to increase
%       of the lower loss.
        c = Q(1:q);
        others.c = c;
        for i=1:q
            a = rand;
            output_tr(c(i)) = mod(output_tr(c(i))+ceil(9*a),10);
        end
%         this should be named as validation set instead.
        input_te = fmnist_image(:,R(5001:10000));
        output_te = fmnist_label(R(5001:10000),1);
        [n1,m1] = size(input_te);
        spinput_te = sparse(input_te);
        
%         classify the labels
        Btr = classify(output_tr);
        Btr =(Btr==1);
        
        Bte = classify(output_te);
%         bte = sum(Bte);
        Bte =(Bte==1);
        gtr = [repmat(input_tr,10,1);ones(10,m2)];
        gte = [repmat(input_te,10,1);ones(10,m1)];
        sumyI_te = getsumyI(input_te,Bte,m2,n2);
        gtr_nonzero = (gtr~=0);
        
        const = 4e-4;
        F = @(x,y) F_CEloss(x,y,input_te,Bte,m1,n1)*const;
        f = @(x,y) f_CEloss(x,y,input_tr,Btr,m2,n2)*const;
%         others.getp_te = @(W) getp(W,input_te,m1,n1);
%         others.getp_tr = @(W) getp(W,input_tr,m2,n2);
        diff.Fx = @(x,y) zeros(m2,1);
        diff.Fy = @(x,y) diffFy_CE(x,y,input_te,m1,n1,sumyI_te)*const;
        diff.fy = @(x,y) difffy_CE(x,y,input_tr,Btr,m2,n2,gtr)*const;
        diff.Fxx = @(x,y,v) zeros(m2,1);
        diff.Fxy = @(x,y,v) zeros(m2,1);
        diff.Fyx = @(x,y,v) zeros(10*n1+10,1);
        diff.Fyy = @(x,y,v) diffFyy_CE(x,y,v,input_te,m1,n1,spinput_te)*const;
        diff.fx = @(x,y) difffx_CE(x,y,input_tr,Btr,m2,n2)*const;
        diff.fxx = @(x,y,v) diff.fx(x,y).*ddsigmoid(x).*v*const;
        diff.fxy = @(x,y,v) difffxy_CE(x,y,v,input_tr,Btr,m2,n2,spinput_tr)*const;
        diff.fyx = @(x,y,v) difffyx_CE(x,y,v,input_tr,Btr,m2,n2,spinput_tr)*const;
        diff.fyy = @(x,y) difffyy_CE(x,y,input_tr,m2,n2,gtr,gtr_nonzero)*const;
        diff.fyyv = @(x,y,v) difffyyv_CE(x,y,v,input_tr,m2,n2,spinput_tr)*const;
end
end
%% auxiliary functions
%% logistic BLO model
function y = F_logistic(lambda,X,input_te,output_te,m)
% y = 0;
temp = exp(-output_te.*(input_te*X));
temp = g1_logistic(temp);
y = sum(temp);
% for i=1:m
%     temp = exp(-output_te(i)*input_te(i,:)*X);
%     temp2 = 1/(1+temp);
%     y = y+temp2;
% end
end
function y = f_logistic(lambda,X,input_tr,output_tr,m)
y = exp(lambda)*norm(X)^2;
temp = exp(-output_tr.*(input_tr*X));
temp = g1_logistic(temp);
y = sum(temp)+y;
% for i=1:m
%     temp = exp(-output_tr(i)*input_tr(i,:)*X);
%         temp2 = 1/(1+temp);
%     y = y+temp2;
% end
end
function y = diffFy_logistic(lambda,X,input_te,output_te,m,n)
temp = exp(-output_te.*(input_te*X));
temp = g2_logistic(temp).*output_te;
temp = ones(n,1)*temp';
temp = temp.*input_te';
y = sum(temp,2);
% y = 0;
% for i=1:m
%     temp = exp(-output_te(i)*input_te(i,:)*X);
%     temp2 = temp/(1+temp)^2;
%     y = y+temp2*output_te(i)*input_te(i,:)';
% end
end
function y = diffFyy_logistic(lambda,X,v,input_te,output_te,m)
y = 0;
temp = exp(-output_te.*(input_te*X));
temp = g3_logistic(temp).*(output_te.^2);
for i =1:m
    y = y + temp(i)*input_te(i,:)'*(input_te(i,:)*v);
end
% y = 0;
% for i=1:m
%     temp = exp(-output_te(i)*input_te(i,:)*X);
%     temp2 = (temp^2-temp)/(1+temp)^3;
%     y = y+temp2*output_te(i)^2*input_te(i,:)'*input_te(i,:);
% end
end
function y = difffy_logistic(lambda,X,input_tr,output_tr,m,n)
y = 2*exp(lambda)*X;
temp = exp(-output_tr.*(input_tr*X));
temp = g2_logistic(temp).*output_tr;
temp = ones(n,1)*temp';
temp = (temp.*input_tr');
y = sum(temp,2)+y;
% for i=1:m
%     temp = exp(-output_tr(i)*input_tr(i,:)*X);
%     temp2 = temp/(1+temp)^2;
%     y = y+temp2*output_tr(i)*input_tr(i,:)';
% end
end
function y = difffyy_logistic(lambda,X,input_tr,output_tr,m,n)
y = 2*exp(lambda)*eye(n,n);
temp = exp(-output_tr.*(input_tr*X));
temp = g3_logistic(temp).*(output_tr.^2);
for i =1:m
    y = y + temp(i)*input_tr(i,:)'*input_tr(i,:);
end
% for i=1:m
%     temp = exp(-output_tr(i)*input_tr(i,:)*X);
%     temp2 = (temp^2-temp)/(1+temp)^3;
%     y = y+temp2*output_tr(i)^2*input_tr(i,:)'*input_tr(i,:);
% end
end
function y = difffyyv_logistic(lambda,X,v,input_tr,output_tr,m,n)
y = 2*exp(lambda)*v;
temp = exp(-output_tr.*(input_tr*X));
temp = g3_logistic(temp).*(output_tr.^2);
for i =1:m
    y = y + temp(i)*input_tr(i,:)'*(input_tr(i,:)*v);
end
end
function y = g1_logistic(x)
y = log(1+x);
end
function y = g2_logistic(x)
y = -x./(1+x);
end
function y = g3_logistic(x)
y = x./(1+x).^2;
end
%% mnist
function y = classify(label)
y = zeros(length(label),10);
for i = 1:10
    y(:,i) = (label==(i-1));
end
end
function y = getp(W,X,m,n)
% here W is a vector with dimension of 784*10+10, X is the image
b = W(end-9:end);
W_ori = reshape(W(1:end-10),n,10)';
temp = W_ori*X+b*ones(1,m);
temp = exp(temp);
y = temp./(ones(10,1)*sum(temp,1));
% we would modify it if illness appears.
end
function y = getsumyI(input,B,m,n)
y = zeros(10*(n+1),1);
for i=1:10
    y(1+n*(i-1):n*i,1) = sum(input(:,B(:,i)),2);
    y(10*n+i,1) = sum(B(:,i));
end
end
% function y = weightedsumyI(weight,input,B,m,n)
% y = zeros(10*(n+1),1);
% for i=1:10
% %     y(1+n*(i-1):n*i,1) = sum((ones(n,1)*weight(B(:,i))).*input(:,B(:,i)),2);
%     y(1+n*(i-1):n*i,1) =input(:,B(:,i))*weight(B(:,i))';
%     y(10*n+i,1) = sum(weight(B(:,i)));
% end
% end

function y = F_CEloss(lambda,W,input_te,Bte,m1,n1)
P = getp(W,input_te,m1,n1);
y = 0;
for i=1:10
    y = y + sum(-log(P(i,Bte(:,i))));
end
end
function y = f_CEloss(lambda,W,input_tr,Btr,m2,n2)
b = sigmoid(lambda)';
P = getp(W,input_tr,m2,n2);
y = 0;
for i=1:10
    y = y + sum(-b(Btr(:,i)).*log(P(i,Btr(:,i))));
end
end
function y = diffFy_CE(lambda,W,input_te,m1,n1,sumyI_te)
P = getp(W,input_te,m1,n1);
y = zeros(n1*10+10,1);
for i =1:10
    y(1+n1*(i-1):i*n1) = input_te*P(i,:)';
end
y(10*n1+1:end) = sum(P,2);
y = y - sumyI_te;
end
function y = difffy_CE(lambda,W,input_tr,Btr,m2,n2,g)
P = getp(W,input_tr,m2,n2);
y = zeros(n2*10+10,1);
weight = sigmoid(lambda)';
wP = P.*(ones(10,1)*weight);%%
for i = 1:10
    b = Btr(:,i);
    wP(i,b) = wP(i,b) - weight(b);
    y(1+n2*(i-1):i*n2) = input_tr*wP(i,:)';
    y(10*n2+i) = sum(wP(i,:));
end
end
function y =difffx_CE(lambda,W,input_tr,Btr,m2,n2)
y = dsigmoid(lambda);
P = getp(W,input_tr,m2,n2);
for i =1:10
    b = Btr(:,i);
    y(b) = -y(b).*log(P(i,b)');
end
end
function y = diffFyy_CE(lambda,W,v,input_te,m1,n1,spinput_te)
y = zeros(10*n1+10,1);
P = getp(W,input_te,m1,n1);
% modified on 2022/1/17
% w2 = v'*(g.*[kron(P,ones(n1,1));P]);%
% gp = computegp(spinput_te,P,m1,n1);
% w2 = v'*gp;
[w2,z] = computegp(spinput_te,P,m1,n1,v,'yy');
% for i=1:10
%     w1 = P(i,:).*(v(1+n1*(i-1):n1*i)'*input_te+v(10*n1+i)-w2);
% %     y(1+n1*(i-1):n1*i) = sum(input_te*(ones(n1,1)*w1),2);
%     y(1+n1*(i-1):n1*i) = input_te*w1';
%     % this is added on 12/16
%     y(i+10*n1) = sum(w1);
% end
W1 = z+P.*((v(10*n1+1:end)*ones(1,m1))-(ones(10,1)*w2));
s = spinput_te*W1';
    % for i=1:10
    %      y(1+n1*(i-1):n1*i) = s(:,i);
    % end
y(1:10*n1)=reshape(s,10*n1,1);
y(1+10*n1:end) = sum(W1,2);
end
function y = difffxy_CE(lambda,W,v,input_tr,Btr,m2,n2,spinput_tr)
temp =  dsigmoid(lambda);
P = getp(W,input_tr,m2,n2);
% modified on 2022/1/17
% w = v'*(g.*[kron(P,ones(n2,1));P]);
% gp = computegp(spinput_tr,P,m2,n2);
% w = v'*gp;
[w,~] = computegp(spinput_tr,P,m2,n2,v,'yy');
y = w';
for i =1:10
    b =Btr(:,i);
    y(b) = y(b)-((v(1+(i-1)*n2:n2*i)'*spinput_tr(:,b))'+v(n2*10+i)*ones(sum(b),1));
end
y = y.*temp;
end
function y = difffyx_CE(lambda,W,v,input_tr,Btr,m2,n2,spinput_tr)
temp =  dsigmoid(lambda);
v = v.*temp;
P = getp(W,input_tr,m2,n2);
% modified on 2022/1/17
% w = (g.*[kron(P,ones(n2,1));P])*v;
% gp = computegp(spinput_tr,P,m2,n2,v,'yx');
% w = gp*v;
[w,~] = computegp(spinput_tr,P,m2,n2,v,'yx');
y = w;
for i =1:10
    b = Btr(:,i);
    y(n2*(i-1)+1:n2*i) = y(n2*(i-1)+1:n2*i) -input_tr(:,b)*v(b);
%     modified on 12/17
    y(10*n2+i) = y(10*n2+i)-sum(v(b));
end
end
function y = difffyy_CE(lambda,W,input_tr,m2,n2,g,g_nonzero)
weight = sigmoid(lambda);
P = getp(W,input_tr,m2,n2);
y = zeros(10*n2+10,10*n2+10);
wP = (ones(10,1)*weight').*P;
wG = g.*[kron(wP,ones(n2,1));wP];
G = g.*[kron(P,ones(n2,1));P];
D = zeros(10*n2+10,n2);
for i =1:m2
    b = g_nonzero(:,i);
    gib = G(b,i);
    wgib = wG(b,i);
    y(b,b) = y(b,b) - wgib*gib';
    bb = b(1:n2);
    D(b,bb) = wgib*input_tr(bb,i)'+D(b,bb);
%     hi = input_tr(:,i);
%     Pi = P(:,i);
%     D = generateD(hi,Pi);
%     y = y + weight(i)*D;
end
for i=1:10
    y(1+n2*(i-1):n2*i,1+n2*(i-1):n2*i)=y(1+n2*(i-1):n2*i,1+n2*(i-1):n2*i)+D(1+n2*(i-1):n2*i,:);
    y(10*n2+i,1+n2*(i-1):n2*i) = y(10*n2+i,1+n2*(i-1):n2*i) + D(10*n2+i,:);
end
y(1+n2*(i-1):n2*i,10*n2+i) = y(10*n2+i,1+n2*(i-1):n2*i);
y(10*n2+1:end,10*n2+1:end) = y(10*n2+1:end,10*n2+1:end)+diag(sum(wP,2));
% y = y+eye(10*n2,10*n2)*2e-3;% regulization
end

function y = difffyyv_CE(lambda,W,v,input_tr,m2,n2,spinput_tr)
weight = sigmoid(lambda);
y = zeros(10*n2+10,1);
P = getp(W,input_tr,m2,n2);
% modified on 2022/1/17
% w2 = v'*(g.*[kron(P,ones(n2,1));P]);%
% gp = computegp(spinput_tr,P,m2,n2);
% w2 = v'*gp;
[w2,z] = computegp(spinput_tr,P,m2,n2,v,'yy');
% for i=1:10
% %     w1 = P(i,:).*(v(1+n2*(i-1):n2*i)'*input_tr+v(10*n2+i)-w2);
%     w1 = z(i,:)+P(i,:)*v(10*n2+i)-w2.*P(i,:);
% % %     y(1+n2*(i-1):n2*i) = sum(input_tr*(ones(n2,1)*(w1.*weight)),2);
%     w1 = w1'.*weight;
%     y(1+n2*(i-1):n2*i) = input_tr*(w1);
% % %     y(1+n2*(i-1):n2*i) = input_tr*(w1'.*weight);
%     % modified on 12/17
%     y(10*n2+i) = sum(w1);
% end
W1 = z+P.*((v(10*n2+1:end)*ones(1,m2))-(ones(10,1)*w2));
W1 = W1.*(ones(10,1)*weight'); W1 = W1';
s = spinput_tr*W1;
for i=1:10
     y(1+n2*(i-1):n2*i) = s(:,i);
end
y(1+10*n2:end) = sum(W1)';
end

function y = sigmoid(x)
y = 1./(1+exp(-x));
end
function y = dsigmoid(x)
temp = exp(-x);
y = temp./(1+temp).^2;
end
function y = ddsigmoid(x)
% compute the sigmoid(x)''/sigmoid(x)'
temp = exp(-x);
y = (temp-1)./(1+temp);
end
% added on 2022/1/17 by ch
function [y,z] = computegp(spinput,P,m,n,v,mod)%modified on 2022/1/21
% z = zeros(11,m);
% for i=1:10
% %     y(1+n*(i-1):n*i,:) = full(spinput.*(ones(n,1)*P(i,:)));
%     z(i,:) = (v(1+n*(i-1):n*i)'*spinput).*P(i,:);
% end
% % y(10*n+1:10*n+10,:) = P;
% z(11,:) = v(1+n*10:end)'*P;
% y = sum(z);
if strcmp(mod,'yy')
%     v1 = [v(1:n)';v(1+n:2*n)';v(1+2*n:3*n)';v(1+3*n:4*n)';v(1+4*n:5*n)';v(1+5*n:6*n)';v(1+6*n:7*n)';v(1+7*n:8*n)';v(1+8*n:9*n)';v(1+9*n:10*n)'];
    v1 = reshape(v(1:10*n),n,10)';
    v2 = v(1+10*n:end)';
    z = (v1*spinput).*P;
    y = sum(z)+v2*P;
elseif strcmp(mod,'yx')
    z = spinput*(P'.*(v*ones(1,10)));
    y = reshape(z,10*n,1);
    y = [y;P*v];
end
end
