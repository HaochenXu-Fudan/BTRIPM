function y = Opts(funname,algname)
opts = struct();
switch funname
    case{'20news'}
        switch algname
            case{'cg'}
            opts.s1 = 0.01;  
            opts.s2 = 0.01; 
            opts.K = 1;   
            opts.T = 50;  
            opts.J = 200;  
            opts.epsilon1 = 1e-4;   
            opts.epsilon2 = 1e-3;  
            case{'rhg'}
            opts.s = 0.05;  
            opts.s1 = 0.06;   
            opts.s2 = 0.01;  
            opts.K = 500;   
            opts.T = 50;  
            opts.Ty = 100;   
            opts.epsilon1 = 1e-6;  
            opts.epsilon2 = 1e-6; 
            opts.epsilon3 = 1e-6;
            case{'trhg'}
            opts.s = 0.05;  
            opts.s1 = 0.06;   
            opts.s2 = 0.01;  
            opts.K = 250;   
            opts.T = 50;  
            opts.Ty = 100;   
            opts.epsilon1 = 1e-6;  
            opts.epsilon2 = 1e-6; 
            opts.epsilon3 = 1e-6;
            case{'ITD'}
            opts.s1 = 0.01;   
            opts.J = 200;   
            opts.epsilon1 = 1e-4;   
            opts.alpha=0.01; 
            opts.t=10;
            opts.K=50;
            case{'bvfim'}
            opts.epsilon1 = 2e-1;
            opts.epsilon2 = 5e-2;
            opts.u = @(k) 10/1.01^k;
            opts.zj = 0;
            opts.tau0 = 1;
            opts.eta = 1.01;
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
            opts.s1 = 0.1;
            opts.s2 = 0.01;
            opts.alpha = 0.01;
            opts.zj=1;
        case{'Tregion'}
            opts.epsilon1 = 1e-6;
            opts.epsilon2 = 2e-1;
            opts.u = @(k) 4/1.2^k;
            opts.zj = 0;
            opts.tau0 = 1;
            opts.eta = 1.1;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 50;
            opts.Tz = 50;
            opts.s1 = 0.4;
            opts.maxit = 10;
            opts.Delta = 2e-1;
            opts.hessian = 'nonaccurate';%
        case{'Tregionzj'}
            opts.epsilon1 = 1e-6;
            opts.epsilon2 = 3e-3;
            opts.u = @(k) 8/1.2^k;
            opts.zj=1;
            opts.tau0 = 1;
            opts.eta = 1.1;
            opts.flag = 1;
            opts.L = 1;
            opts.K = 50;
            opts.Tz = 50;
            opts.s1 = 0.1;
            opts.maxit = 6;
            opts.Delta = 1e-2;
            opts.hessian = 'nonaccurate';%
        end
end
y = opts;
end