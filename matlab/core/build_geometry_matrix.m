function H = build_geometry_matrix(candidatePool,rRxECEF_m)
%BUILD_GEOMETRY_MATRIX Build the pseudorange geometry matrix.
%
% Each row is
%
%   h_i^T = [u_i^T  1]
%
% where u_i is the line-of-sight unit vector from the
% transmitter toward the receiver.

rRxECEF_m = rRxECEF_m(:);

N = height(candidatePool);
H = zeros(N,4);

for i = 1:N

    rTxECEF_m = 1000 * [ ...
        candidatePool.ECEF_X_km(i);
        candidatePool.ECEF_Y_km(i);
        candidatePool.ECEF_Z_km(i)];

    delta = rRxECEF_m - rTxECEF_m;

    u = delta / norm(delta);

    H(i,:) = [u.' 1];
end

end