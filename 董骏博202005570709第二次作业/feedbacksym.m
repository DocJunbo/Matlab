function GB = feedbacksym(G,H,sign)
    if nargin == -2
        sign = -1;
        GB = G/(sym(1) - sign * G * H);
        GB = simplify(GB);
    end