function nlp = build_opti(nseg,t,A,e,w)
%BUILD_OPTI Build the symbolic computational graph for optimization
    tm = [t.^0; t.^1; t.^2; t.^3; t.^4]';

    opti = casadi.Opti();

    % -- 4th order spline coefficients
    ax = opti.variable(5,nseg);
    ay = opti.variable(5,nseg);
    % -- Time point discretization for each segment
    tp = opti.variable(nseg*(numel(t)-1));
    % -- Slack variables for C1 continuity
    sl = opti.variable(nseg-1);
    
    wp_x = opti.parameter(nseg+1,1);
    wp_y = opti.parameter(nseg+1,1);
    x_lim = opti.parameter();
    y_lim = opti.parameter();
    t_lim = opti.parameter();
    v_max = opti.parameter();
    a1 = opti.parameter();
    a2 = opti.parameter();
    
    % -- Assure C0 continuity of spline
    for i = 1:nseg
        opti.subject_to( ax(1,i) == wp_x(i) );
        opti.subject_to( ay(1,i) == wp_y(i) );
        opti.subject_to( sum(ax(:,i)) == wp_x(i+1) );
        opti.subject_to( sum(ay(:,i)) == wp_y(i+1) );
    end

    % -- Soft C1 continuity
    for i = 1:nseg-1
        opti.subject_to( ax(2,i+1)-(0:4)*ax(:,i)+sl(i) == 0 );
        opti.subject_to( ay(2,i+1)-(0:4)*ay(:,i)+sl(i) == 0 );
    end
    opti.subject_to( { sl(:) >= 0, sl(:) <= 0.05 } );
    
    % -- Constrain trajectory to flow domain
    x_vals = tm*ax; y_vals = tm*ay;
    opti.subject_to( { x_vals(:) >= 0, x_vals(:) <= x_lim } );
    opti.subject_to( { y_vals(:) >= 0, y_vals(:) <= y_lim } );

    % -- Constrain velocity
    % `x_vals` and `y_vals` are matrices, each column representing the
    % corresponding coordinates. First element in (i+1)th column and last
    % element in ith column are equal. We get rid of the last row and 
    % flatten the resulting truncated matrix for imposing velocity 
    % constriant
    x_vals_flat = x_vals(1:end-1,:);
    x_vals_flat = x_vals_flat(:);
    y_vals_flat = y_vals(1:end-1,:);
    y_vals_flat = y_vals_flat(:);
    for i = 2:length(tp)
        v = sum_square([(x_vals_flat(i)-x_vals_flat(i-1))/(tp(i)-tp(i-1)); 
                        (y_vals_flat(i)-y_vals_flat(i-1))/(tp(i)-tp(i-1))]);
        opti.subject_to( { v >= 0, v <= v_max^2 } );
    end

    % -- Constrain timepoints to be positive and monotonically increasing
    opti.subject_to( tp(:) >= 0 );
    for i = 1:length(tp)-1
        opti.subject_to( tp(i) < tp(i+1) );
    end
    opti.subject_to( tp(end) <= t_lim );

    % -- Compute actuation energy cost
    E = 0;
    for i = 2:length(tp)
        v = [(x_vals_flat(i)-x_vals_flat(i-1))/(tp(i)-tp(i-1));
             (y_vals_flat(i)-y_vals_flat(i-1))/(tp(i)-tp(i-1))];
        [~,~,~,ux_,uy_] = double_gyre(x_vals_flat(i),y_vals_flat(i),tp(i),A,e,w);
        E = E + sum_square(v - [ux_; uy_]);
    end

    obj = a1*tp(end) + a2*E;
    opti.minimize( obj );

    opti.set_initial(tp, linspace(0,100,numel(tp))');
    opti.set_initial(ax, ones(5,nseg));
    opti.set_initial(ay, ones(5,nseg));

    options = struct;
    options.ipopt.tol = 0.05;
    options.ipopt.print_level = 4;

    opti.solver('ipopt',options);
    
    nlp = struct;
    
    nlp.opti = opti;

    nlp.vars.ax = ax;
    nlp.vars.ay = ay;
    nlp.vars.tp = tp;
    
    nlp.params.wp_x = wp_x;
    nlp.params.wp_y = wp_y;
    nlp.params.x_lim = x_lim;
    nlp.params.y_lim = y_lim;
    nlp.params.t_lim = t_lim;
    nlp.params.v_max = v_max;
    nlp.params.a1 = a1;
    nlp.params.a2 = a2;

    nlp.traj_x = x_vals_flat;
    nlp.traj_y = y_vals_flat;

    nlp.obj = obj;
end
