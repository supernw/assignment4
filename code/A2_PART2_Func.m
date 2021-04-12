function [cMap,Vmap,Ex,Ey,Jx,Jy,Current, Box_top_rec, Box_bottom_rec] = A2_PART2_Func(L_i, W_i, Lb_i, Wb_i, condIn, V_BC)

L = L_i*1e9;
W = W_i*1e9;

nx = L; 
ny = W;

Lb = Lb_i*1e9;
Wb = Wb_i*1e9;

Lc = L/2; %Centre of length
Wc = W/2; %Centre of width
Box_top = [(Lc-(Lb/2)) (Lc + (Lb/2)) (W - Wb) W]; %x1 x2 y1 y2
Box_top_rec = [ Box_top(1) Box_top(3) (Box_top(2)-Box_top(1)) (Box_top(4)-Box_top(3))]; %xywh

Box_bottom = [(Lc-(Lb/2)) (Lc + (Lb/2)) 0 (0 + Wb)]; %x1 x2 y1 y2
Box_bottom_rec = [ Box_bottom(1) Box_bottom(3) (Box_bottom(2)-Box_bottom(1)) (Box_bottom(4)-Box_bottom(3))]; %xywh

%Conductivity Map
cMap = zeros(nx,ny);
%condIn = 10^(-2); %Inside Box
condOut = 1; %Outside Box
for i = 1:nx
    for j = 1:ny
        if( (i >= Box_top(1) && i <= Box_top(2) ) && (j >= Box_top(3) && j <= Box_top(4))) %Inside Top Box
            cMap(i,j) = condIn;   
        elseif ( (i >= Box_bottom(1) && i <= Box_bottom(2) ) && (j >= Box_bottom(3) && j <= Box_bottom(4))) %Inside Bottom Box
            cMap(i,j) = condIn;
        else
            cMap(i,j) = condOut;
        end
    end
end

G = sparse(nx*ny);
F = zeros(1,nx*ny);

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny; %Node Mapping
      
        if i == 1 %Set Left Voltage to 1
           %G(n,:) = 0;
           G(n,n) = 1;
           F(n) = V_BC;
        elseif i == nx %Set Right Voltage to 0
            %G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 0;
        elseif j == 1 %Bottom - Insulated
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nyp = n + 1;
            
            %Average centre conductivity with surrounding
            rxm = (cMap(i,j) + cMap(i-1,j))/2;
            rxp = (cMap(i,j) + cMap(i+1,j))/2;
            ryp = (cMap(i,j) + cMap(i,j+1))/2;
            
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp; 
            
            G(n,n) = -(rxm + rxp + ryp);
        elseif j == ny
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = n - 1;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2;
            rxp = (cMap(i,j) + cMap(i+1,j))/2;
            rym = (cMap(i,j) + cMap(i,j-1))/2;
            
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            
            G(n,n) =-(rxm + rxp + rym); %Only 3 resistors
        else %Bulk
            nxp = j + i*ny;
            nxm = j + (i-2)*ny;
            nyp = n + 1;
            nym = n - 1;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2;
            rxp = (cMap(i,j) + cMap(i+1,j))/2;
            rym = (cMap(i,j) + cMap(i,j-1))/2;
            ryp = (cMap(i,j) + cMap(i,j+1))/2;
            
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            G(n,n) = -(rxm + rxp + rym + ryp);  
        end
    end
end

V = G\F';

%Reverse Node Mapping
Vmap = zeros(nx,ny);
for i = 1:nx
    for j =1:ny
        n = j+(i-1)*ny;
        Vmap(i,j) = V(n);
    end
end


%Ex and Ey
[Ex, Ey] = gradient(Vmap');

Jx = -Ex.*cMap';
Jy = -Ey.*cMap';

Jxt = Jx';
C0 = sum(Jxt(1,:)); %Current at x=0
Cnx = sum(Jxt(nx,:)); %Current at x = nx
Current = (C0+Cnx)/2;

end

