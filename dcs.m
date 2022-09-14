clc; close all; clear all;

x_lim = 2; nx = 50;
y_lim = 1; ny = 25;
t_lim = 40;

dx = x_lim/nx;
dy = y_lim/ny;
dt = 0.01;

xspan = 0:dx:x_lim-dx;
yspan = 0:dy:y_lim-dy;
tspan = 0:dt:t_lim;

nt = length(tspan);

%% Flow Data

Af = 0.1;     % controls velocity magnitude
e = 0.25;     % dictates magnitude of oscillation in x-direction
w = 2*pi/10;  % angular oscillation frequency

f_flatten = @(mat3) reshape(mat3,numel(mat3(:,:,1)),[]);

[Xgrid,Ygrid,~,ux,uy] = double_gyre(xspan,yspan,tspan,Af,e,w);

X = vertcat(f_flatten(ux),f_flatten(uy));
u_bar = mean(X,2);  % temporal mean velocity
X = X - u_bar;      % centered flow velocity
fprintf('-- Built data matrix\n');

%% Proper Orthogonal Decomposition (POD)

fprintf('-- Running POD...\n');
% `A`   - Time coefficients
% `PHI` - Spatial modes (will be used as transform basis as well)
% `SIG` - Eigenvalues of spatial correlation matrix
tic;
[PHI,SIG,~] = svd(X);
A = X'*PHI;
SIG = diag(SIG).^2;
toc;

nmodes = 5;  % first 5 modes are sufficient for us

fprintf('-- Using r = %d POD modes\n',nmodes);

%% Optimal Waypoint Selection

nwaypoints = 5;
ntrials = 20;
wp_error = zeros(1,ntrials);
iwaypoints = zeros(nwaypoints,ntrials);

for i = 1:ntrials
    fprintf('-- Waypoint selection (%d/%d)\n',i,ntrials);
    C = randn(2*nwaypoints,size(PHI,1));  % random measurement matrix
    iloc = randi(nx*ny, [nwaypoints 1]);  % sensing locations
    y = X([iloc; iloc+(nx*ny)],1);        % extend indices for ux and uy
    s = cs_get_sparse(PHI,C,y);

    % -- Pick locations corresponding to `nwaypoints` largest entries in `s`
    [~,iloc] = sort(s,'descend');
    iloc = iloc(1:nwaypoints);

    % -- Reconstruct flow at selected sensor locations and compute error
    Xr = 0.*X;
    for k = 1:nmodes
        Xr(iloc,:) = Xr(iloc,:) + PHI(iloc,k)*A(:,k)';
    end
    wp_error(i) = norm(Xr(iloc,:) - X(iloc,:),2);
    iwaypoints(:,i) = iloc;
end

% -- Select waypoints with minimal reconstruction error
[~,imin] = min(wp_error);
y = iwaypoints(:,imin);
y(y>nx*ny) = y(y>nx*ny) - nx*ny;

%% Trajectory Optimization

% -- Parameters
x_begin = 0.10; y_begin = 0.05;
x_end = 1.90; y_end = 0.95;
ntrials = 1;
v_max = 0.7;
a1 = 1;  % penalty for duration
a2 = 1;  % penalty for energy cost

f_shuffle = @(vec) vec(randperm(length(vec)));

Xgrid_flat = f_flatten(Xgrid);
Ygrid_flat = f_flatten(Ygrid);

nseg = nwaypoints+1;

t = 0:0.05:1;  % spline parameter

prev_obj_val = inf;
traj_ax = zeros(5,nseg);
traj_ay = zeros(5,nseg);

nlp = build_opti(nseg,t,Af,e,w);

for i = 1:ntrials
    fprintf('-- Trajectory Optimization (%d/%d)\n',i,ntrials);
    y = f_shuffle(y);

    wp_x = [x_begin; Xgrid_flat(y); x_end];
    wp_y = [y_begin; Ygrid_flat(y); y_end];

    wp_x = [x_begin; 1.5; 1.55; 0.75; 0.3; 1.2; x_end;];
    wp_y = [y_begin; 0.15; 0.75; 0.4; 0.8; 0.2; y_end;];

    nlp.opti.set_value(nlp.params.wp_x, wp_x);
    nlp.opti.set_value(nlp.params.wp_y, wp_y);
    nlp.opti.set_value(nlp.params.x_lim, x_lim);
    nlp.opti.set_value(nlp.params.y_lim, y_lim);
    nlp.opti.set_value(nlp.params.v_max, v_max);
    nlp.opti.set_value(nlp.params.a1, a1);
    nlp.opti.set_value(nlp.params.a2, a2);

    try
        sol = nlp.opti.solve();

        if sol.value(nlp.obj) < prev_obj_val
            prev_obj_val = sol.value(nlp.obj);
            traj_ax = sol.value(nlp.vars.ax);
            traj_ay = sol.value(nlp.vars.ay);
            traj_x = sol.value(nlp.traj_x);
            traj_y = sol.value(nlp.traj_y);
            traj_tp = sol.value(nlp.vars.tp);
        end
    catch ME
        continue;
    end
end

%% Flow Reconstruction

x_seg = traj_x(traj_tp < 0.75*t_lim);
y_seg = traj_y(traj_tp < 0.75*t_lim);

iloc = round(ny*x_seg/dx) + round(y_seg/dy) + 1;
iloc = unique(iloc,'stable');
iloc = [iloc; iloc+nx*ny];
iloc = iloc(iloc <= 2*nx*ny);

yr = 0.*X(:,1);
for k = 1:nmodes
    yr(iloc) = yr(iloc) + PHI(iloc,k)*A(end,k)';
end
yr = yr(iloc);

C = randn(numel(yr),size(PHI,1));
s = cs_get_sparse(PHI,C,yr);
ur = PHI*s + u_bar;
urx_flat = ur(1:nx*ny); ury_flat = ur(nx*ny+1:end);
urx = reshape(urx_flat,[ny nx]); ury = reshape(ury_flat,[ny nx]);
ux_true = reshape(X(1:nx*ny,end) + u_bar(1:nx*ny),[ny nx]);
emat = urx - ux_true;

% -- Normalize
f_norm = @(mat) -1+2.*(mat-min(mat,[],'all'))./...
                (max(mat,[],'all')-min(mat,[],'all'));

ux_true = f_norm(ux_true);
urx = f_norm(urx);

%% Plots

close all;

% -- POD mode energy distribution plot
tke = 100*SIG/sum(SIG);
figure;
plot(tke(1:10),'o--k','MarkerFaceColor','k','MarkerSize',4), grid on;
xlabel('Mode','FontWeight','bold');
ylabel('Energy (%)','FontWeight','bold');
axis([1 10 0 100]);
xline(nmodes,'-r');
text(nmodes+0.1,60,sprintf('r = %d',nmodes));
cum_energy = sum(tke(1:nmodes));
cum_energy = fix(100*cum_energy)/100;
text(nmodes+0.1,55,sprintf('Cum. Energy = ~%.2f%%',cum_energy));
title('Contribution of First 10 POD modes to Energy','fontsize',12);
set(gcf,'Position',[200 200 800 600]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,'figures/pod_mode_energy_dist.png');

% -- Waypoint selection - reconstruction error plot
[wp_error,isorted] = sort(wp_error,'descend');
figure;
plot(1:size(wp_error,2), wp_error, 's--k','MarkerFaceColor','k'), grid on;
xlabel('Trials','FontWeight','bold');
ylabel('Error','FontWeight','bold');
xlim([1 size(wp_error,2)]);
title('Waypoint Trials v/s Reconstruction Error','FontSize',12);
set(gcf,'Position',[200 200 800 600]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,'figures/waypoint_selection.png');

% -- Trajectory plot
t = 0:0.01:1;  % spline parameter
tm = [t.^0; t.^1; t.^2; t.^3; t.^4]';
x_vals = tm*traj_ax; y_vals = tm*traj_ay;
figure; colormap winter;
quiverC2D(Xgrid,Ygrid,ux(:,:,1),uy(:,:,1),'Colorbar',true);
hold on; box on;
scatter(wp_x(1),wp_y(1),'sk','MarkerFaceColor','k');
scatter(wp_x(end),wp_y(end),'^k','MarkerFaceColor','k');
scatter(wp_x(2:end-1),wp_y(2:end-1),10,'ok','MarkerFaceColor','k');
text(wp_x(2:end-1),wp_y(2:end-1)-0.03,num2str((1:nwaypoints)'));
plot(x_vals(:),y_vals(:),'k');
xlabel('X','FontWeight','bold');
ylabel('Y','FontWeight','bold');
axis([0 x_lim 0 y_lim]);
title('Optimized Trajectory','FontSize',12);
pbaspect([x_lim y_lim 1]);
set(gcf,'Position',[200 200 960 540]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,'figures/trajectory_C1.png');

% -- Reconstruction
figure; colormap jet;
plt = pcolor(Xgrid,Ygrid,ux_true);
plt.EdgeColor = 'none';
plt.FaceColor = 'interp';
xticks(0:0.5:x_lim); yticks(0:0.5:y_lim);
title(sprintf('Reference Flow Map (u-component) at t = %.2fs',...
      tspan(end)),'FontSize',12);
pbaspect([x_lim y_lim 1]); colorbar;
set(gcf,'Position',[200 200 960 540]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,'figures/flow_u_ref.png');

figure; colormap jet;
plt = pcolor(Xgrid,Ygrid,urx);
plt.EdgeColor = 'none';
plt.FaceColor = 'interp';
xticks(0:0.5:x_lim); yticks(0:0.5:y_lim);
title(sprintf('Reconstructed Flow Map (u-component) at t = %.2fs',...
      tspan(end)),'FontSize',12);
pbaspect([x_lim y_lim 1]); colorbar;
set(gcf,'Position',[200 200 960 540]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,'figures/flow_u_recons.png');

figure; colormap jet;
plt = pcolor(Xgrid,Ygrid,emat);
plt.EdgeColor = 'none';
plt.FaceColor = 'interp';
xticks(0:0.5:x_lim); yticks(0:0.5:y_lim);
title(sprintf(...
      'Reconstruction Error Distribution (u-component) at t = %.2fs',...
      tspan(end)),'FontSize',12);
pbaspect([x_lim y_lim 1]); colorbar;
set(gcf,'Position',[200 200 960 540]);
set(gcf,'PaperPositionMode','auto');
saveas(gcf,'figures/flow_u_error.png');