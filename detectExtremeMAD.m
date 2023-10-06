function [t_, x_, up] = detectExtremeMAD(t, x, par)
    m = zeros(size(x));     % m is sliding-window median value
    for i = 1 : length(t)
        temp = x(t>=t(i)-par.tWin/2 & t<t(i)+par.tWin/2);
        m(i) = median(temp);
    end
    x_ = x-m;               % median-removed time series
    % x_ = (x_ - mean(x_, 'omitnan')) / std(x_, 'omitnan');
%     x_ = x_ - median(x_, 'omitnan');
%     mad_ = median(abs(x_), 'omitnan');
    x_ = x_ - median(x_(~isnan(x_)));
    mad_ = median(abs(x_(~isnan(x_))));
    x_ = x_/mad_;%*erfinv(1/2)*sqrt(2);
    up = x_ > par.thr;      % logical values of threshold crossing
    t_ = {}; n = 0; tLast = -Inf;
    for i = 1 : length(up)  % scan for events
        if up(i) && t(i)-tLast > par.ref
            n = n + 1;
            t_{n} = t(i);
            tLast = t(i);
        end
    end
    t_ = cell2mat(t_);