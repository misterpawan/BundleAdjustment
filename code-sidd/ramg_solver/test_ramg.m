% Testing, Recursive V-cycle Algebraic Multigrid Method on BA problem

function test_ramg()

	load JTJ.dat            % loading approximate hessian matrix csr format

   	A = spconvert(JTJ);     % converting the data in matrix form
    	[m,n] = size(A);        % finding the size of JTJ matrix
    	b = rand(m,1);          % taking right hand side b to mx1 vector

	[x, vc_cnt, rel_res] = ramg(A, b, @GS_Iter, 1, 1, 1e-10, 20);

end

