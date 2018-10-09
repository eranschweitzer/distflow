function p = unravel_path(s,t,pi)
%%% return nodes in the path from s to t given predecessor vector pi

    p = [];
    pred = t;
    while pred ~= s
        p = [pred,p]; %#ok<AGROW>
        pred = pi(pred);
    end
    p = [s,p];
end