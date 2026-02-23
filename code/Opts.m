function opts = Opts(funname,alg_name)
opts = struct();
if strcmp(funname,'toy0')||strcmp(funname,'toy2')
    switch alg_name
        case{'bvfim'}
            opts.epsilon1 = 1e-2;
            opts.epsilon2 = 1e-1;
            opts.u = @(k) 10/1.02^k;
            % if x0 = 0, y0 =0, we can set both u of bvfim and pro are 10/1.04^k
            opts.tau0 = 1;
            opts.eta = 1.01;
            opts.flag = 0;
            opts.L = 1;
            opts.K = 500;
            opts.Tz = 50;
            opts.Ty = 25;
            opts.s1 = 0.1;
            opts.s2 = 0.01;
            opts.alpha = 0.01;
        case{'rhg'}
            opts.epsilon1 = 1e-6; 
            opts.epsilon2 = 1e-6; 
            opts.epsilon3 = 1e-6; 
            opts.s = 0.1;
            opts.s1 = 0.1;
            opts.s2 = 0.1;
            opts.K = 500;
            opts.T = 100;
            opts.Ty = 100;
        case{'cg'}
            opts.s1 = 1e-2;
            opts.s2 = 1e-2; 
            opts.K = 500;
            opts.T = 100;
            opts.epsilon1 = 1e-2;
            opts.epsilon2 = 1e-2;
        case{'pro'}
            opts.epsilon1 = 2e-3;
            opts.epsilon2 = 5e-2;
            opts.u = @(k) 60/1.04^k;
            opts.tau0 = 1;
            opts.eta = 1.02;
            opts.flag = 0;
            opts.L = 1;
            opts.K = 300;
            opts.Tz = 50;
            opts.Ty = 25;
            opts.s1 = 0.1;
            opts.s2 = 0.01;
            opts.alpha = 0.01;
        case{'Tregion'}
            opts.epsilon1 = 1e-8;
            opts.epsilon2 = 1e-4;
            opts.u = @(k) 4/1.4^k;
            opts.tau0 = 1;
            opts.eta = 1.02;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 100;
            opts.Tz = 100;
            opts.s1 = 0.1;
            opts.maxit = 5;
    end
elseif strcmp(funname,'wifi')
    switch alg_name
        case{'bvfim'}
            opts.epsilon1 = 1e-3;
            opts.epsilon2 = 1e-2;
            opts.u = @(k) 4/1.01^k;
            opts.tau0 = 1;
            opts.eta = 1.01;
            opts.flag = 0;
            opts.L = 1;
            opts.K = 500;
            opts.Tz = 50;
            opts.Ty = 25;
            opts.s1 = 0.1;
            opts.s2 = 0.01;
            opts.alpha = 0.01;
        case{'rhg'}
            opts.epsilon1 = 1e-10; 
            opts.epsilon2 = 1e-10; 
            opts.epsilon3 = 1e-10; 
            opts.s = 0.1;
            opts.s1 = 0.1;
            opts.s2 = 0.1;
            opts.K = 500;
            opts.T = 100;
            opts.Ty = 100;
        case{'trhg'}
            opts.epsilon1 = 1e-10; 
            opts.epsilon2 = 1e-10; 
            opts.epsilon3 = 1e-10; 
            opts.s = 0.1;
            opts.s1 = 0.1;
            opts.s2 = 0.1;
            opts.K = 500;
            opts.T = 50;
            opts.Ty = 50;
        case{'cg'}
            opts.s1 = 0.01; 
            opts.s2 = 0.01; 
            opts.K = 500;
            opts.T = 100;
            opts.epsilon1 = 1e-2;
            opts.epsilon2 = 1e-2;
        case{'Tregion'}
            opts.epsilon1 = 1e-6;
            opts.epsilon2 = 1e-3;
            opts.u = @(k) 4/1.5^k;
            opts.tau0 = 1;
            opts.eta = 1.02;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 20;
            opts.Tz = 50;
            opts.s1 = 0.1;
            opts.maxit = 5;
    end
elseif strcmp(funname,'crowdsourced')
    switch alg_name
        case{'bvfim'}
            opts.epsilon1 = 1e-3;
            opts.epsilon2 = 1e-2;
            opts.u = @(k) 4/1.02^k;
            opts.tau0 = 1;
            opts.eta = 1.01;
            opts.flag = 0;
            opts.L = 1;
            opts.K = 500;
            opts.Tz = 50;
            opts.Ty = 25;
            opts.s1 = 0.5;
            opts.s2 = 0.01;
            opts.alpha = 0.01;
         case{'rhg'}
            opts.epsilon1 = 1e-10; 
            opts.epsilon2 = 1e-10; 
            opts.epsilon3 = 1e-10; 
            opts.s = 0.5;
            opts.s1 = 0.01;
            opts.s2 = 0.01;
            opts.K = 500;
            opts.T = 100;
            opts.Ty = 100;
         case{'trhg'}
            opts.epsilon1 = 1e-10; 
            opts.epsilon2 = 1e-10; 
            opts.epsilon3 = 1e-10; 
            opts.s = 0.5;
            opts.s1 = 0.01;
            opts.s2 = 0.01;
            opts.K = 500;
            opts.T = 50;
            opts.Ty = 50;
        case{'cg'}
            opts.s1 = 0.01; 
            opts.s2 = 0.01; 
            opts.K = 500;
            opts.T = 100;
            opts.epsilon1 = 1e-4;
            opts.epsilon2 = 1e-4;
        case{'Tregion'}
            opts.epsilon1 = 1e-6;
            opts.epsilon2 = 2e-1;
            opts.u = @(k) 4/1.2^k;
            opts.tau0 = 0.8;
            opts.eta = 1.3;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 25;
            opts.Tz = 50;
            opts.s1 = 0.2;
            opts.maxit = 5;
%             opts.hessian = 'nonaccurate';%    
    end
elseif strcmp(funname,'mnist5000')
    switch alg_name
        case{'bvfim'}
            opts.epsilon1 = 2e-1;
            opts.epsilon2 = 5e-2;
            opts.u = @(k) 6/1.02^k;
            opts.zj = 0;
            opts.tau0 = 1;
            opts.eta = 1.03;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 500;
            opts.Tz = 50;
            opts.Ty = 25;
            opts.s1 = 0.4;
            opts.s2 = 0.01;
            opts.alpha = 0.01;
        case{'bvfimzj'}
            opts.epsilon1 = 1e-3;
            opts.epsilon2 = 1e-2;
            opts.u = @(k) 4/1.02^k;
            opts.tau0 = 1;
            opts.eta = 1.01;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 500;
            opts.Tz = 50;
            opts.s1 = 0.4;
            opts.s2 = 0.01;
            opts.alpha = 0.01;
            opts.zj=1;
        case{'Tregion'}
            opts.epsilon1 = 1e-6;
            opts.epsilon2 = 3e-1;
            opts.u = @(k) 10/1.2^k;
            opts.zj = 0;
            opts.tau0 = 1;
            opts.eta = 1.3;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 50;
            opts.Tz = 50;
            opts.s1 = 0.4;
             opts.maxit = 6;
            opts.Delta = 1e-2;
            opts.hessian = 'nonaccurate';%
        case{'Tregionzj'}
            opts.epsilon1 = 1e-6;
%             opts.epsilon2 = 3e-3;
            opts.epsilon2 = 3e-3;
            opts.u = @(k) 8/1.2^k;
            opts.zj=1;
            opts.tau0 = 1;
%             opts.eta = 1.3;
            opts.eta = 1.1;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 50;
            opts.Tz = 50;
            opts.s1 = 0.4;
%             opts.maxit = 8;
             opts.maxit = 6;
%             opts.Delta = 0.25;
            opts.Delta = 1e-2;
%             opts.Delta = 2;
            opts.hessian = 'nonaccurate';%
%             opts.hessian = 'guassian';
    end
elseif strcmp(funname,'fashion_mnist')
    switch alg_name
        case{'bvfim'}
            opts.epsilon1 = 2e-1;
            opts.epsilon2 = 5e-2;
            opts.u = @(k) 5/1.02^k;
            opts.zj=0;
            opts.tau0 = 1;
            opts.eta = 1.03;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 500;
            opts.Tz = 50;
            opts.Ty = 25;
            opts.s1 = 0.1;
            opts.s2 = 0.01;
            opts.alpha = 0.01;
        case{'bvfimzj'}
            opts.epsilon1 = 2e-1;
            opts.epsilon2 = 1e-2;
            opts.u = @(k) 5/1.02^k;
            opts.zj=1;
            opts.tau0 = 1;
            opts.eta = 1.01;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 500;
            opts.Tz = 50;
            opts.Ty = 25;
            opts.s1 = 0.1;
            opts.s2 = 0.01;
            opts.alpha = 0.01;
        case{'Tregionzj'}
            opts.epsilon1 = 1e-6;
            opts.epsilon2 = 3e-3;
            opts.u = @(k) 6/1.2^k;
            opts.zj=1;
            opts.tau0 = 1;
            opts.eta = 1.3;
%             opts.eta = 1.1;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 50;
            opts.Tz = 50;
            opts.s1 = 0.1;
            opts.maxit = 6;
%             opts.Delta = 0.25;
            opts.Delta = 1e-2;
            opts.hessian = 'nonaccurate';
%             opts.hessian = 'guassian';
        case{'Tregion'}
            opts.epsilon1 = 1e-6;
            opts.epsilon2 = 3e-1;
            opts.u = @(k) 6/1.2^k;
            opts.zj=0;
            opts.tau0 = 1;
            opts.eta = 1.3;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 50;
            opts.Tz = 50;
            opts.s1 = 0.1;
            opts.maxit = 6;
%             opts.Delta = 0.25;
            opts.Delta = 1e-2;
            opts.hessian = 'nonaccurate';
    end

end
end