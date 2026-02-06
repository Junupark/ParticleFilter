function next_x = Euler(f, t, x, h)
    arguments
        f function_handle
        t
        x
        h {mustBePositive}
    end
    f1 = f(t,x);
    next_x = x + f1*h;
end
