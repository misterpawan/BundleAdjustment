% Function for Gauss Seidel as a smoother - used for pre and post smoothing

% A  : [in]  Fine grid problem matrix - J^TJ (Hessian Matrix)
% b  : [in]  Right hand side
% x :  [out]  the solution after gauss siedel iteration
% iter_cnt : [out] iteration counter

function [x, iter_cnt] = GS_Iter(A, b, iter_eps, max_iter, x0)
	n = size(A, 1);
	iter_cnt = 0;
	err = 1;
	x = x0;
	
	while ((err >= iter_eps) && (iter_cnt < max_iter))
		x0 = x;
		
		for i = 1 : n
			A_ii = A(i, i);
			x(i) = b(i) - A(i, :) * x + A_ii * x(i);
			x(i) = x(i) / A_ii;
		end
		
		iter_cnt = iter_cnt + 1;
		err = max(abs(x0 - x));
	end
end

