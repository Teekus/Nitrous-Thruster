% Solver for ALL MY EQUATIONS
close all
clear all

tic

%%%%%%%%%%%%%%%%%%%% Video Code %%%%%%%%%%%%%
window_parameters =  [ 50          74        1184         876];

video_record = false;
% SETUP VIDEO STUFF
videoFrameRate = 25;    %frame rate of the video
Directory_path = pwd;   %current path
video_title = '0.32_exit';

%the following is code that automatically creates a videos folder in the
%working folder (if it doesn't exist) and automatically gives your video a
% unique title that includes a date and time (but the time is hard to read)

if video_record
    if ~exist('videos','dir')   %check if the 'videos' folder in the current path does not exist
        mkdir(Directory_path,'videos')  %make the videos folder
        disp(['Created new folder ''videos'' in ',Directory_path,'... That''s where your simulation videos will be.'])
    end
    disp(['Creating new .avi video in ',Directory_path,'\videos'])
    video_string = [Directory_path,'\videos\',video_title,'_',datestr(now,1),'_',datestr(now,30)];  %give video unique name with date and time (so you don't overwrite old videos)
    writerObj = VideoWriter(video_string,'MPEG-4');  %create the video object, .avi
    set(writerObj, 'FrameRate',videoFrameRate)  %set the framerate
    open(writerObj);    %open the video
end




% Gas Properties (PASS THESE IN)
rho_N2O    = double(1.98);
Ru         = double(1.987);
u          = double(0.05);
k          = double(40/1000);
var_u      = int32(1);
u_min      = double(u);
u_max      = double(0.27);

disp('Function is trying to reach steady_state for lowest sigma at a new length')
% MATERIAL PROPERTIES (PASS THESE IN)
dens_M     = double(0.12);
mat        = 'Copper';
sigma      = double(433);
h          = double(10);
kb         = double(401*dens_M*0.33);
rho_b      = double(dens_M*8940);
cp_b       = double(390);


% Initilization

% For velocity, must us the smaller alpha
Average_CP = 38.6891;
Average_MW = 0.0440;

% Determining which is the limiting factor in alpha

L          = double(0.0635);

alpha_b    = kb./(rho_b*cp_b);
alpha_g    = k./(rho_N2O*(Average_CP/Average_MW));
alpha_min = min(alpha_b,alpha_g);
alpha_max = max(alpha_b,alpha_g);

dx      = double(alpha_min/(2.5*u_max)); % inserted a randam factor in here
x_nodes = double(floor(((L-1e-10)/dx) + 1));


if x_nodes < 500
    x_nodes = double(500);
    dx     = double((L-1e-10)/(x_nodes - 1));
end

dt         = double((0.5*dx^2)/(alpha_max));
final_time = 200;

t_steps    = double(floor((final_time/dt) + 1));
t_mesh = linspace(0, final_time, t_steps);




% x_mesh_old = x_mesh;
x_mesh = double(linspace(1e-10, L, x_nodes));



%%% THIS IS FOR CONTINUING SOLUTIONS AT HIGHER RESOLUTIONS IF NEEDED
% 
% % Interpolation for a higher resolution.
% for mesh_in = 1:x_nodes
%     gas_temp(mesh_in) = interp1(x_mesh_old, gas, x_mesh(mesh_in));
%     block_temp(mesh_in) = interp1(x_mesh_old, block, x_mesh(mesh_in));
%     dens_temp(mesh_in) = interp1(x_mesh_old, dens, x_mesh(mesh_in));
%     gasp_temp(mesh_in) = interp1(x_mesh_old, gas_p, x_mesh(mesh_in));
%     blockp_temp(mesh_in) = interp1(x_mesh_old, block_p, x_mesh(mesh_in));
%     densp_temp(mesh_in) = interp1(x_mesh_old, dens_p, x_mesh(mesh_in));
% 
% end
% 
% gas = double(gas_temp);
% block = double(block_temp);
% dens = double(dens_temp);
% 
% %
% gas_p = double(gasp_temp);
% block_p = double(blockp_temp);
% dens_p = double(densp_temp);



%%%% HEAT LOSS, SET TO 0 FOR NOW
% r      = 0.01;
%
% % heat_loss
% k_q = 3;
% L_q = 0.001;


%%%% INITIALIZING ALL VARIABLES TO BE PASSED IN
plot_step = int32(2000);
size_pl      = floor(t_steps/plot_step);

%
gas_pl    =  double(300*ones(1,size_pl*x_nodes));
block_pl  =  double(1200*ones(1,size_pl*x_nodes));
dens_pl   =  double(0*ones(1,size_pl*x_nodes));

gas    = double(300*(ones(1,x_nodes)));
% gas    = double(linspace(300,1200, x_nodes)); 
block  = double(1200*(ones(1,x_nodes)));
dens   = double(0*(ones(1,x_nodes)));

gas_p   = gas;
block_p = block;
dens_p  = dens;

% gas stuff
del_2_term = ones(1,x_nodes);
hg_b       = ones(1,x_nodes);
convec     = ones(1,x_nodes);
decomp     = ones(1,x_nodes);
beta       = ones(1,x_nodes);

cp         = double(ones(1, x_nodes)); 

% block terms

del_2_termb = ones(1,x_nodes);
hb_g        = ones(1,x_nodes);

% Density Terms

del_rho     = ones(1, x_nodes);
rho_change  = ones(1, x_nodes);


fac = double(0.13);
limit = double(2);
fac2 = double(0.75);

count  = int32(zeros(1,2)); 

num_vec = int32(500);
u_vec   = double(linspace(u_min, u_max, num_vec)); 
time_to_change = int32(1000);
main_u_vec = double(ones(1,size_pl));


disp('Array structures and variables have been defined')
disp('Entering main loop....')

% mex MEX_MATLAB_VELOCITYCHANGE.c

disp(' ');

disp('Successfully mexed matlab file and now moving into main loop') 


disp(' ');
disp(' ');
 
tic
MEX_MATLAB_VELOCITYCHANGE(rho_N2O, u, k, sigma, h, kb, ...
    rho_b, cp_b, cp, L, dx, dt, ...
    x_nodes, t_steps, x_mesh, plot_step, count, fac, limit, fac2, ...
    gas_pl, block_pl, dens_pl, gas, block,...
    dens, var_u ,u_min,u_max,time_to_change,u_vec,main_u_vec, num_vec);

toc
save('CU_BaseProperties_Lowest_SIGMA_MAX_DENSITY_0635length_027')
disp('Function is trying to reach steady_state for lowest sigma at a new length')
count = count(1) + (count(1) - length(main_u_vec));
count = length(main_u_vec); 
gas_pl = reshape(gas_pl,[ x_nodes,count]); 
gas_pl = gas_pl';
gas_pl(end,:) =  gas; 
% 
block_pl = reshape(block_pl,[ x_nodes,count]); 
block_pl = block_pl'; 
block_pl(end, :) = block;
% 
dens_pl = reshape(dens_pl,[ x_nodes,count]);
dens_pl = dens_pl';
dens_pl(end,:) =  dens;

Main_plot(x_mesh, gas_pl, block_pl, dens_pl,50,main_u_vec,kb,k,sigma,mat,video_record)
% plot(gas)
% 
% mex MEX_MATLAB.c
% 
% cp_single = double(cp);
% MEX_MATLAB(double(rho_N2O), double(u), double(k), double(sigma), double(h), double(kb), ...
%     double(rho_b), double(cp_b), double(L), double(dx), double(dt), ...
%     int32(x_nodes), int32(t_steps), double(x_mesh), int32(plot_step), double(fac), double(limit), double(fac2), ...
%     double(gas_pl), double(block_pl), double(dens_pl), double(gas), double(block),...
%     double(dens), cp_single, int32(0),int32(0),int32(0),int32(0));