function xp = netCS_Noise2(A,y,lambda, excite_idx, inhib_idx)
% A is a measurement matrix
% y is a results vector (measure of response to stim at patched cell)
% lambda is a regularization term that trades off sparsity of solution
% for closeness of fit. Larger lambda produces sparses solutions. default
% to 1/4
% xp is the solution to the minimization problem which represents
% connection strength
% min and max values for connection strength optional

%%%%%%%%%%%%%%%%%%%%%%%%
% Problem definition
% minimize    lambda*||x||(L-1) + ||A*x - b||(L-2)
%%%%%%%%%%%%%%%%%%%%%%%%

% solve the LP
% tic
% m trials, n cells
%     subject to
%         x >= -20;
%         x <= 10;
[m, n] = size(A);
cvx_begin quiet
    variable x(n);
    minimize(0.5*norm(A*x - y,2) + lambda*norm(x,1));
    subject to
        x(excite_idx) >= 0;
        x(inhib_idx) <= 0;
cvx_end
xp = x;
end