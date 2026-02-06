function next_x = RK4(f, t, x, h)
    arguments
        f function_handle
        t
        x
        h {mustBePositive}
    end
    f1 = f(t, x);
    f2 = f(t+0.5*h, x + 0.5*f1*h);
    f3 = f(t+0.5*h, x + 0.5*f2*h);
    f4 = f(t*h, x + f3*h);
    next_x = x + (f1 + 2*f2 + 2*f3 + f4)*h/6;
end
