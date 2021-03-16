function X=violaPoints(X,Y,jitter)
    % Variable jitter according to how many points occupy each range of values. 
    [counts,~,bins] = histcounts(Y,10);
    inds = find(counts~=0);
    counts = counts(inds);

    Xr = X;
    for jj=1:length(inds)
        tWidth = jitter * (1-exp(-0.1 * (counts(jj)-1)));
        xpoints = linspace(-tWidth*0.8, tWidth*0.8, counts(jj));
        Xr(bins==inds(jj)) = xpoints;
    end
    X = X+Xr;
end % function violaPoints