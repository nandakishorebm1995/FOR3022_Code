function K = Tlaplacematrix2d(dt,N)
%% 2D-LAPLACE MATRIX (K) AND REDUCED-ORDER BASIS COMPUTATION


ondiag = -2;
offdiag = 1;


D = sparse(1:N,1:N,ondiag*ones(1,N),N,N);
E = sparse(2:N,1:N-1,offdiag*ones(1,N-1),N,N);
K = (E+D+E')./(dt^2);

end