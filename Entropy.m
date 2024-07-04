function H = Entropy(data, normalized)
    if sum(data) == 0
        H = 0;
        return
    end
    if ~normalized
        data = data/sum(data);
    end
    logdata = log(data);
    logdata(isinf(logdata)) = 0;
    H = -sum(data.*logdata);
end