function avg_curr = PART23_Func(delV, Wb)

L = 200e-9;
W = 100e-9;
Lb = 40e-9;
condIn = 10^-2;
[cMap,Vmap,Ex,Ey,Jx,Jy,Current, Box_top_rec, Box_bottom_rec] = A2_PART2_Func(L, W, Lb, Wb, condIn, delV);

figure('Name', 'V(x,y)')
surf(Vmap','linestyle', 'none');
view(0,90)
title('Voltage - V(x,y) Graph') 
xlabel('Length [nm]');
ylabel('Width [nm]');
rectangle('Position', Box_top_rec)
rectangle('Position', Box_bottom_rec)

figure('Name', 'E Field')
quiver(-Ex,-Ey)
title('Electric Field - E(x,y) Graph')
xlabel('Length [nm]');
ylabel('Width [nm]');

%MONTE CARLO CODE

%Features
scattering = 1;
boxes = 1;
random_velocity = 1;
electric_field = 1;

%Constants
m_o = 9.11e-31; %Rest masss of electron [kg]
m_ef = 0.26*m_o; %Effective mass of an electron
kB = 1.3806e-23; %Bolzmann's constant [m^2 kg s^-2 k^-1]
T_lat = 300; %Lattice Temperature [k]
q_e = 1.602e-19; %Charge of an electron [C]

v_th = sqrt((2*kB*T_lat)/(m_ef));

region_x = 200e-9; %Width of region [m]
region_y = 100e-9; %Height of region [m]

num_p = 1000; %Number of particles
num_dt = 1000; %Number of time steps
num_plot = 10; %Number of particles to plot

dt = ((1/200)*region_x)/v_th; %Time step [sec]
tow_mn = 0.2e-12; %Mean time between collisions [sec]
%tow_mn = dt*10;
dx(1,:) = zeros(1,num_p); %X Position change
dx(1,:) = zeros(1,num_p); %X Position change

%Position and Velocity Arrays
Px(1,:) = zeros(1,num_p); %X-coordinate
Py(1,:) = zeros(1,num_p); %Y-coordinate
Vx(1,:) = zeros(1,num_p); %X-velocity
Vy(1,:) = zeros(1,num_p); %Y-velocity

%Boxes
Lc = L/2;
Wc = W/2;

Box1 = [(Lc - Lb/2) (Lc + Lb/2) 0 Wb 1]; %x1 x2 y1 y2 BC (1 = spec, 0 = diff)
Box1_rec = [ Box1(1) Box1(3) (Box1(2)-Box1(1)) (Box1(4)-Box1(3))]; %xywh
Box2 = [(Lc - Lb/2) (Lc + Lb/2) (W-Wb) W 1]; %x1 x2 y1 y2 BC (1 = spec, 0 = diff)
Box2_rec = [ Box2(1) Box2(3) (Box2(2)-Box2(1)) (Box2(4)-Box2(3))]; %xywh

%Top and bottom spectral vs. diffusive
top_bc = 1; %0 = diff, 1 = spec
bottom_bc = 1; %0 = diff, 1 = spec
%Generate colours for number of plotted particles
colours = zeros(3,num_plot);
for d=1:num_plot
    colours(:,d) = [rand;rand;rand];
end
%colours = {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'};

V_rand_store(1,:) = zeros(1,num_p); 

for j=1:num_p
   %Assign each particle a random position to start not inside the boxes
   keep_trying = 1;
   
   while(keep_trying == 1)
      Px(j) = rand*region_x;
      Py(j) = rand*region_y;
      if((Px(j) > Box1(1) && Px(j) < Box1(2)) && (Py(j) < Box1(4))) %Inside box1
          keep_trying = 1;
      elseif((Px(j) > Box2(1) && Px(j) < Box2(2)) && (Py(j) > Box1(3))) %Inside box1 %Inside Box2
          keep_trying = 1;
      else
          keep_trying = 0;
      end
      
      if(boxes == 0) %Turn off keep trying if boxes off
          keep_trying = 0;
      end
   end
 
   %Random velocity from dist with random direction
   %Random Velocity
   v_rand = sqrt(kB*T_lat/m_ef)*randn+v_th; %Normal Dist
   V_rand_store(j) = v_rand;
   theta = rand*360; %random direction
   
   if(random_velocity == 1)
       Vx(j) = v_rand*cosd(theta);
       Vy(j) = v_rand*sind(theta);
   else
       Vx(j) = v_th*cosd(theta);
       Vy(j) = v_th*sind(theta);
   end
   
end
figure('Name', 'Velocity Histogram')
histogram(V_rand_store);
xlabel('Velocity (m/s)')
ylabel('Number of Particles with Velocity')
title('Particle Velocity Histogram')
%Next position of plotted particles
Px_next = zeros(1,num_p);
Py_next = zeros(1,num_p);

%Old position of plotted particles
Px_old = zeros(1,num_p);
Py_old = zeros(1,num_p);
T = zeros(1,num_p);
T_avg = zeros(1,num_dt);
time = zeros(1,num_dt);

%Time since last scattered
mft = zeros(1,num_p);
avg_time_array = zeros(1,num_p);

%Acceleration
ax(1,:) = zeros(1,num_p);
ay(1,:) = zeros(1,num_p);

%Concentration
n = 10^15 / 0.0001; %[cm^2] to [m^2]
Jx = zeros(1,num_dt);
time_curr = zeros(1,num_dt);
running_time = 0;

figure('Name', 'P2: Particle Traj')
for i=1:num_dt
    %Max Velocity
    %{
    mag_v = sqrt((Vx + ax.*mft).^2 + (Vy + ay.*mft).^2);
    max_v = max(mag_v);
    dt = ((1/200)*region_x)/(max_v); %Time step [sec]
    %}
    
    %Current
    Jx(i) = q_e.*n.*mean(Vx);
    time_curr(i) = running_time;
    running_time = running_time + dt;
    
    Pscat = 1 - exp(-(dt)/(tow_mn)); %Scattering Probability Function
    %Scattering
    rand_scat = rand(1,num_p);

    Iscat = find(Pscat > rand_scat);
    for z=1:length(Iscat)
        v_rand = sqrt(kB*T_lat/m_ef)*randn+v_th; %Normal Dist
        theta = rand*360; %random direction
        Vx(Iscat(z)) = v_rand*cosd(theta);
        Vy(Iscat(z)) = v_rand*sind(theta);

        if(avg_time_array(Iscat(z)) ~= 0) %Not first value
            avg_time_array(Iscat(z)) = mean([avg_time_array(Iscat(z)) mft(Iscat(z))]);
        else
            avg_time_array(Iscat(z)) = mft(Iscat(z));
        end
        mft(Iscat(z)) = 0;
    end
    
    %Get next position
    Vx_next = Vx + ax.*dt;
    Px_next = Px + (Vx_next).*dt;
    Vy_next = Vy + ax.*dt;
    Py_next = Py + (Vy_next).*dt;
    
%Boundary Conditions
    %Left and Right
    Px(Px_next > region_x) = 0;
    Px(Px_next < 0) = region_x;
    
    %Top
    Vy(Py >= region_y) = -1*(Vy(Py >= region_y));
    
    %Bottom
    Vy(Py <=0) = -1*(Vy(Py <=0));

    %Box 1
    Px_th_lr_box1 = (Px_next >= Box1(1) & Px < Box1(1));
    Px_th_rl_box1 = (Px_next <= Box1(2) & Px > Box1(2));
    
    Py_th_td_box1 = (Py_next <= Box1(4) & Py > Box1(4));
    
    Px_in_x_box1 = (Px_next > Box1(1) & Px_next < Box1(2)); %If x component in box
    Py_in_y_box1 = (Py_next > Box1(3) & Py_next < Box1(4)); %If y component in box
    from_out_x_box1 = (Px < Box1(1) | Px > Box1(2)); 
    from_in_x_box1 = (Px >= Box1(1) & Px <= Box1(2));
    from_bet_y_boxs = (Py >= Box1(4) & Py <= Box2(3));
    
    bounce_x_box1 = Px_in_x_box1 & Py_in_y_box1 & from_out_x_box1 & (Px_th_lr_box1 | Px_th_rl_box1);
    %bounce_y_box1 = Px_in_x_box1 & Py_in_y_box1 & from_in_x_box1 & Py_th_td_box1;

    %bounce_x_box1 = Px_in_x_box1 & Py_in_y_box1 & from_out_x_box1;
    bounce_y_box1 = Px_in_x_box1 & Py_in_y_box1 & from_in_x_box1 ;
    
    %Negate Velocity
    Vx(bounce_x_box1) = -1*(Vx(bounce_x_box1));
    Vy(bounce_y_box1) = -1*(Vy(bounce_y_box1));

    %Change position to boundary
    Px(bounce_x_box1 & Px < Box1(1)) = Box1(1); %from left
    Px(bounce_x_box1 & Px > Box1(2)) = Box1(2); %from right
    Py(bounce_y_box1) = Box1(4);
    
    %Box2
    Px_th_lr_box2 = (Px_next >= Box2(1) & Px < Box2(1));
    Px_th_rl_box2 = (Px_next <= Box2(2) & Px > Box2(2));
    
    Py_th_bu_box2 = (Py_next >= Box2(3) & Py < Box2(3));
    
    Px_in_x_box2 = (Px_next > Box2(1) & Px_next < Box2(2)); %If x component in box
    Py_in_y_box2 = (Py_next > Box2(3) & Py_next < Box2(4)); %If y component in box
    from_out_x_box2 = (Px < Box2(1) | Px > Box2(2)); 
    from_in_x_box2 = (Px >= Box2(1) & Px <= Box2(2));
    
    bounce_x_box2 = Px_in_x_box2 & Py_in_y_box2 & from_out_x_box2 & (Px_th_lr_box2 | Px_th_rl_box2);
    %bounce_y_box2 = Px_in_x_box2 & Py_in_y_box2 & from_in_x_box2 & Py_th_bu_box2;

    %bounce_x_box2 = Px_in_x_box2 & Py_in_y_box2 & from_out_x_box2;
    bounce_y_box2 = Px_in_x_box2 & Py_in_y_box2 & from_in_x_box2;
    %Negate Velocity
    Vx(bounce_x_box2) = -1*(Vx(bounce_x_box2));
    Vy(bounce_y_box2) = -1*(Vy(bounce_y_box2));

    %Change position to boundary
    Px(bounce_x_box2 & Px < Box2(1)) = Box2(1); %from left
    Px(bounce_x_box2 & Px > Box2(2)) = Box2(2); %from right

    Py(bounce_y_box2) = Box2(3);
    
    %Update Acceleration
    for g=1:num_p
        x = round(Py(g)/1e-9);
        y = round(Px(g)/1e-9);
        if(x <= 0)
            x = 1;
        elseif(x > 100)
            x = 100;
        end

        if(y <= 0)
            y = 1;
        elseif (y > 200)
            y = 200;
        end
        ax(g)= (-(Ex(x,y)/1e-9)*q_e)/m_ef;  
        ay(g) = (-(Ey(x,y)/1e-9)*q_e)/m_ef;

    end
    
    %Update all position
    Px_old = Px;
    Vx = Vx + ax.*dt;
    Px = Px + (Vx).*dt;
    Py_old = Py;
    Vy = Vy + ay.*dt;
    Py = Py + (Vy).*dt;
    
    mft = mft + dt; %Add time since last bounce or scatter
    
    %Calculate the semiconductor temperature
    v_avg = mean(sqrt(Vx.^2 + Vy.^2));
    T_avg(i) = (v_avg^2*m_ef)/(2*kB);
    
    %Calculate the temperature of each electron
    T = ((sqrt(Vx.^2 + Vy.^2)).^2.*m_ef)./(2*kB);
    
    time(i) = i*dt;
    %Plot subet of particles
    for k=1:num_plot
        hold on;
        plot([Px_old(k) Px(k)],[Py_old(k) Py(k)], 'Color', colours(:,k));
    end
    if(i == 1)
        if(boxes == 1)
            rectangle('Position', Box1_rec)
            rectangle('Position', Box2_rec)
        end
        xlabel('X Position')
        ylabel('Y Position')
        xlim([0 region_x])
        ylim([0 region_y])
    end
    title(['Particle Trajectories (Avg Temp = ',num2str(T_avg(i)), ' K)'])
    

    fprintf('Iteration = %i \n', i);
    pause(0.05);
    
end

Px_bined = zeros(1,num_p);
Py_bined = zeros(1,num_p);
Num_Px = zeros(1,40);
Num_Py = zeros(1,20);
for i = 1:num_p
    Px_bined(i) = round(Px(i)/5e-9);
    Py_bined(i) = round(Py(i)/5e-9);
    
    for j = 1:40
        if(Px_bined(i) == j)
            Num_Px(j) = Num_Px(j) + 1;
        end
    end
    for k = 1:20
       if(Py_bined(i) == k)
            Num_Py(k) = Num_Py(k) + 1;
       end 
    end
% end

%Electron Density Map
figure('Name', 'P3: Density')
[den, cent] = hist3([Px',Py'], 'Nbins', [40 20]);
[X,Y] = meshgrid(cent{1}, cent{2});
surf(X,Y,den')
colorbar
xlabel('Position X')
ylabel('Position Y')
title('Electron Density')

%Temperature Map
%Bin the electrons
bin_length = 200e-9/40;
bin_width = 100e-9/20;
temp_binned = zeros(40,20);
for i=1:num_p
    for j=1:length(cent{1}) %X bin
         if(j == 1)
            x_in_bin = (Px > 0) & (Px <= (cent{1}(j)+bin_length/2));
         else
            x_in_bin = (Px > (cent{1}(j-1)-bin_length/2)) & (Px <= (cent{1}(j)+bin_length/2));
         end
       for k = 1:length(cent{2}) %Y bin
          if(k == 1)
              y_in_bin = (Py > 0) & (Py <= (cent{2}(k)+bin_width/2));
          else
              y_in_bin = (Py > (cent{2}(k-1)-bin_width/2)) & (Py <= (cent{2}(k)+bin_width/2));
          end
          
          curr_sum = sum(T(x_in_bin & y_in_bin));
          if(curr_sum == 0)
              temp_binned(j,k) = 0;
          else
              temp_binned(j,k) = mean(T(x_in_bin & y_in_bin));
          end
          
       end
    end
end
figure('Name', 'P3: Temp')

surf(X,Y,temp_binned')
colorbar
xlabel('Position X')
ylabel('Position Y')
title('Temperature Map')

%Average Current of all time
avg_curr = mean(Jx);

end

