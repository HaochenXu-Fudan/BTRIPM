function error = errorcompute(funname,input)
error = struct();
switch funname
    case{'toy0'}
        error.F = abs(input.F-pi^2/8);
        error.f = abs(input.f+1);
        error.x = abs(input.x+pi/4);
        error.y = abs(input.y+pi/4);
    case{'toy2'}
        error.F = abs(input.F-2*(3*pi/4-2)^2);
        error.f = abs(input.f+1);
        error.x = abs(input.x-3*pi/4);
        error.y = abs(input.y-3*pi/4);
end
end