function metrics = compute_position_metrics(H,R)
%COMPUTE_POSITION_METRICS Compute rank, PDOP, FIM, and position bound.

N = size(H,1);

assert(size(H,2) == 4, ...
    'H must have four columns.');

assert(isequal(size(R),[N N]), ...
    'H and R dimensions are inconsistent.');

metrics.rankH = rank(H);

if metrics.rankH < 4
    metrics.PDOP       = Inf;
    metrics.J          = nan(4);
    metrics.CRB        = nan(4);
    metrics.alphaCRB_m = Inf;
    return;
end

%% PDOP
Q = (H.'*H) \ eye(4);

metrics.PDOP = ...
    sqrt(trace(Q(1:3,1:3)));

%% Fisher information matrix
J = H.' * (R \ H);

%% Cramer-Rao covariance matrix
C = J \ eye(4);

metrics.J   = J;
metrics.CRB = C;

metrics.alphaCRB_m = ...
    sqrt(trace(C(1:3,1:3)));

end