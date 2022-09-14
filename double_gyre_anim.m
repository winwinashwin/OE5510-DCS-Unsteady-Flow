clc; close all; clear all;

x_lim = 2; nx = 50;
y_lim = 1; ny = 25;
t_lim = 10;

dx = x_lim/nx;
dy = y_lim/ny;
dt = 0.01;

xspan = 0:dx:x_lim;
yspan = 0:dy:y_lim;
tspan = 0:dt:t_lim;

nt = length(tspan);

%% Flow Data

A = 0.1;      % controls velocity magnitude
e = 0.25;     % dictates magnitude of oscillation in x-direction
w = 2*pi/10;  % angular oscillation frequency

[Xgrid,Ygrid,stream_fn,ux,uy] = double_gyre(xspan,yspan,tspan,A,e,w);

%% Plot

close all;

% Double-gyre flow animation
figure;
colormap jet;
set(gcf,'Position',[200 200 960 540]);
set(gcf,'PaperPositionMode','auto');
set(gcf,'Color','w');

filename = 'figures/double_gyre_flow_anim.gif';

for i = 1:10:nt
    contourf(Xgrid,Ygrid,stream_fn(:,:,i),50,'EdgeColor','none');
    colorbar;
    title(sprintf('Double-gyre flow - Stream function (t=%.2fs)',tspan(i)));
    axis([xspan(1) xspan(end) yspan(1) yspan(end)]);
    pbaspect([x_lim y_lim 1]);
    drawnow;
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end
