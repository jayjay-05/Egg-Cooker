%pde solver
clear all;

% define parameters

T_initial = 30;
T_water = 100;

% egg parameters
kval = .541;                  
c = 3165;       %celcius
rho = 1071.2;

%water parameters
%k = 0.6;                  
%c = 4186;       %celcius
%rho = 1000;

alpha = kval/(c*rho);   %thermal diffusivity

%rchicken = 50mm;
%rostrich=150mm;
%rqual= 27 mm
 
R =50/1000;         % radius of egg [m]
t_max = 4000;      %s

%make the mesh
N=10;   % cut space into N sections
M=1750; % cut time  into M sections
dr=R/N; 
dt=t_max/M; % grid spacing


R_range = 0:dr:R;
time_range = 0:dt:t_max;



PDE = zeros(length(time_range), length(R_range)+1);

%Initialize BC
PDE(1, :) = T_initial;
PDE(:,length(PDE(1,:)))= T_water;

F = (alpha*dt)/(R*(dr^2));

length_time = length(PDE(:,1));
length_r = length(PDE(1,:));


for k=1:length_time -1
    for i =2 :length_r-1

       current_val = PDE(k,i);
       next_r_val = PDE(k, i+1);
       prev_r_val = PDE(k,i-1);
       next_t_val = (1-2*F)*current_val + (F*(next_r_val +prev_r_val));
       PDE(k+1,i)= next_t_val;

    end
    PDE(k+1, 1) = PDE(k+1,3);            %ghost point

end


Temp_chart = PDE;
Temp_chart(:,1)=[];


%Find Time taken to cook
time_count =0;
for ind_t =1: length_time
    if(min(Temp_chart(ind_t,:))>=70)
        time_count = time_range(ind_t);
        break
    end

end

cook_time = time_count+10;
cook_time_minutes = cook_time/60;

R_range = R_range./2;

%{
figure(1)
mesh(R_range,time_range,Temp_chart); 
  shading interp
colormap('jet')
xlabel('x'); ylabel('t'); zlabel('T(x,t)'); colorbar
title("Temperature distribution over time")


figure(4)
plot(Temp_chart(1:100:end, :)')
xlabel('x'); 
ylabel('T(x,t)');
grid on
title("Temperature Distribution every 100 timesteps")
%}



%PLOTINGG COOKING EGG
[x,y,z] = sphere(100);

figure('Color', 'w');

for ind_t = 1:size(Temp_chart,1)

    hax = axes('Position',[0 0 1 1]);

    for ii = length(R_range):-1:1
        h = surf(x*R_range(ii),y*R_range(ii),z*0+Temp_chart(ind_t,ii));  % sphere centered at origin
        set(h, 'EdgeColor', 'None');
        hold on;
      
    end

    view(0,-90);
    axis equal;
    set(hax, 'Visible', 'Off', 'CLim', [min(Temp_chart(:)) max(Temp_chart(:))]);
    colorbar
    pause(0.5);

end




