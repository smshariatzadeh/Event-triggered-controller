
function Z_approx = approx_minRPIset(  Ak, W, n_order, alpha)
    % approximate minimal robust positively invarient set (Z)
    % which takes the form of Z = alpha*(W + Ak*W + Ak^2*W + ... Ak^n_ordr*W).
    % where + denotes Minkowski addition.
    %
    %input:
    %   Ak: (closed loop) system matrix 
    %   W : system noise (as Polyhedron )
    %   n : degree of approximate (n >=1)
    %       choose bigger n for best approximation 
    %   alpha :a magnifier constant  (alpha>=1)
    %
    %output:
    %  Zapprox: approximation of minimal robust positively invarient set as
    %           polyhedron
    %
    % We could obtain approximate minimal RPI (Z) by computing an infinite geometric series,
    % which is not practicall to get.
    % we approximate this by trancating the polynomial.

    if (n_order <1)
       error ('degree of approximate must be (n >=1)')
    end
    if (alpha<1)
       error ('magnifier constant must be gearter than one(alpha>=1)') 
    end
    Z_approx = W;

    for n = 1:n_order
        Z_approx = Z_approx + Ak^n*W;
    end
    Z_approx = Z_approx*alpha;
    % which takes the form of Z = alpha*(W + Ak*W + Ak^2*W + ... Ak^n_ordr*W).
    % where + denotes Minkowski addition.
end
