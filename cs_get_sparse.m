function s = cs_get_sparse(PSI,C,y)
%CS_GET_SPARSE Retrieve sparse vector `s` from compressed measurements
cvx_begin quiet;
    variable s(size(PSI,2));
    minimize( norm(s,1) );
    subject to
        C*PSI*s == y;
cvx_end;

fprintf('[CVX] Status = %s. Took %.2f seconds\n',cvx_status,cvx_cputime);
end

