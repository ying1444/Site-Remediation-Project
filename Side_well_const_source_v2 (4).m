clear
clc

format long

X_plume_edge=1105/2; %ft %plume edge dist from center in X
Y_plume_edge=709/2; %ft %plume edge dist from center in Y
X_hot_edge=10; %ft %max concentration area dist from center in X
Y_hot_edge=5; %ft %max concentration area dist from center in Y
Aquifer_thickness=66; %ft

No_well =2; %number of wells
Well_flow=2*192.5; %flow rate per well from gpm to ft3/day

v_darcy=0.0098; %ft/day %aqifer flow rate
Q=No_well*Well_flow/(Y_plume_edge*2*Aquifer_thickness); %ft^3/day %well pumping rate normalized over capture area
porosity=0.3;
seepage_v=(v_darcy)/porosity;
seepage_Q=Q/porosity; %seepage velocity for advection taken as normalized seepage velocity caused by pumps



C_plume_edge=10*28.3; %concentration at plume edge from ug/L to ug/ft3
C_hot_edge=26340*28.3; %concentration at max concentration area from ug/L to ug/ft3


vol_tank=160*0.035/10;%ft3

density_TCE = 1.61*1000*1000*1000*28.3; %ug/ft3
MW_TCE = 131.4; %g/mol
diffusion_TCE=0.0303; %cm2/hr
K_TCE=0.0259*0.787; %ft/day
long_dispersivity= 0.32*((Y_plume_edge*2)^0.83);
trans_dispersivity = 0.3*long_dispersivity;
C_TCE_sat=1100*1000*28.3; %ug/ft

density_ben = 0.876*1000*1000*1000*28.3; %ug/ft3
MW_ben = 78.11; %g/mol
diffusion_ben=0.0666; %cm2/hr
K_ben=0.013*0.787; %ft/day
long_dispersivity= 0.32*((Y_plume_edge*2)^0.83); %ft
trans_dispersivity = 0.3*long_dispersivity;
C_ben_sat=1780*1000*28.3; %ug/ft3

Dx=long_dispersivity*seepage_v;
Dy=0.3*Dx; %ft/day %dispersion coefficient in Y


% solving for imaginary instantaneous contaminant input concentrations and time elapsed
%  to form initial plume
syms CoA_imaginary_initial t_imaginary_initial
eqns = [
   C_plume_edge==((CoA_imaginary_initial/(4*pi*t_imaginary_initial*sqrt(Dx*Dy)))*exp((-((0-seepage_v*t_imaginary_initial)^2)/(4*t_imaginary_initial*Dx))-(Y_plume_edge^2/(4*t_imaginary_initial*Dy)))),   
     C_hot_edge==((CoA_imaginary_initial/(4*pi*t_imaginary_initial*sqrt(Dx*Dy)))*exp((-((0-seepage_v*t_imaginary_initial)^2)/(4*t_imaginary_initial*Dx))-(Y_hot_edge^2/(4*t_imaginary_initial*Dy)))),
 ];
vars =[CoA_imaginary_initial t_imaginary_initial];
[sol_CoA, sol_t,] = solve(eqns,vars);
eval_sol_CoA=eval(sol_CoA); %imaginary initial point contaminant mass input
eval_sol_t=eval(sol_t); %imaginary time elapsed since contaminant mass input

((eval_sol_CoA/(4*pi*eval_sol_t*sqrt(Dx*Dy)))*exp((-((0-seepage_v*eval_sol_t)^2)/(4*eval_sol_t*Dx))-(Y_hot_edge^2/(4*eval_sol_t*Dy)))) %quick check to see if the max concentration of your initial plume is what you want it o be.
i=1;
j=1;

dx=20; 

dy=20;
dt=15; %should be somewhat less than dx/seepage_Q for stability of code

%actual initialization of plume using imaginary CoA and time

    for x=-553:dx:533
    
        for y=-355:dy:355
       
        C_initial(i,j)=((eval_sol_CoA/(4*pi*eval_sol_t*sqrt(Dx*Dy)))*exp((-((x-seepage_v*eval_sol_t)^2)/(4*eval_sol_t*Dx))-(y^2/(4*eval_sol_t*Dy))));
        i=i+1;
        end
    j=j+1;
    i=1;
    end
    
%can input breakpoint right after this to stop code and check out a contour
%plot of initial plume
contourf(C_initial)
colorbar;
 


        
        
i=2;
    j=2;
    t_increment=1;
    

    [m,n]=size(C_initial); %getting size of initial plume matrix (m rows by n columns)

    C=zeros(m+2,n+1); % initiating a larger zero matrix of largersize for boundary imaginary nodes
 
%following code is putting initial plume matrix into larger sero matrix
    % Get sizes
    [rowsBig, columnsBig] = size(C);
    [rowsSmall, columnsSmall] = size(C_initial);
    % Specify upper left row, column of where
    % we'd like to paste the small matrix.
    row1 = 2;
    column1 = 1;
    % Determine lower right location.
    row2 = row1 + rowsSmall - 1;
    column2 = column1 + columnsSmall - 1;
    % See if it will fit.
    if row2 <= rowsBig
    % It will fit, so paste it.
        C(row1:row2, column1:column2) = C_initial;
    else
     % It won't fit
        warningMessage = sprintf('That will not fit.\nThe lower right coordinate would be at row %d, column %d.',...
        row2, column2);
        uiwait(warndlg(warningMessage));
    end
    
well_1=round(m/2-112/dx); %well 1 positino
well_2=round(m/2+112/dx); %well 2 position




[m,n]=size(C); %getting size of new concentration matrix
C_ben(:,:,t_increment)=C*0.27; % first time step will take on the initial plume concentrations
C_TCE(:,:,t_increment)=C*0.73;
dCdt_ben=zeros(m,n,t_increment); %initializing dCdtt matrix with zeros because special midification required

% next block of code is for frame capturing and exporting of videos. ignore
% until results are good.
set(gca, 'nextplot', 'replacechildren');  
caxis manual; 
caxis([0 7e5]);
%v = VideoWriter('10 year no source.avi');
%v.FrameRate=30;
%open(v);
%F(t_increment) = getframe(gcf);



density_NAPL=density_TCE*0.72+density_ben*0.27;

V_trans=(4)^2*Aquifer_thickness;
Mass_total_NAPL= density_NAPL*vol_tank;
mass_ben(t_increment)=0.27*Mass_total_NAPL;
mass_TCE(t_increment)=0.73*Mass_total_NAPL;


xb=(mass_ben(t_increment)/MW_ben)/((mass_ben(t_increment)/MW_ben)+(mass_TCE(t_increment)/MW_TCE));
xt=(mass_TCE(t_increment)/MW_TCE)/((mass_ben(t_increment)/MW_ben)+(mass_TCE(t_increment)/MW_TCE));


AV=1/Aquifer_thickness;


    
Dx=long_dispersivity*seepage_Q;
Dy=0.3*Dx; %ft/day %dispersion coefficient in Y

    
    
while   any(any(C_TCE(:,1:n-1,t_increment)>0.7*28.3))>0 %while concentrations anywhere are above 10 ug/L keep running through loops
       
        
     
        C_ben(:,1,t_increment)=0; %left bound always 0
        C_ben(1,:,t_increment)=C_ben(3,:,t_increment);% top row imaginary nodes are mirror image of 3rd row
        C_ben(m,:,t_increment)=C_ben(m-2,:,t_increment); %bottom row imaginary nodes are mirror images of ttwo rows above
        C_ben(:,n,t_increment)=C_ben(:,n-2,t_increment)+(2*seepage_Q*C_ben(:,n-1,t_increment)*dx/Dx); %no advction or disperssion flux out passed wells.
        
        %C(well_1,n,t_increment) = ((-seepage_Q*(pi*0.25^2)*C(well_1,n-1,t_increment)+seepage_v*C(well_1,n-1,t_increment))*2*dx/Dx)+C(well_1,n-2,t_increment);
        %C(well_2,n,t_increment) = ((-seepage_Q*(pi*0.25^2)*C(well_2,n-1,t_increment)+seepage_v*C(well_2,n-1,t_increment))*2*dx/Dx)+C(well_2,n-2,t_increment);
       
    
    
       
    for j=2:1:n-1 
            
        if j==n-1;
                
                %forced dCdt in last column in each time step to be the
                %same as numbers as a colume before
            for i=2:1:m-1
                if C_ben(i,j,t_increment)>0
                
                    dCdt_ben(i,j,t_increment)=dCdt_ben(i,j-1,t_increment);
                
                    C_ben(i,j,t_increment+1)=C_ben(i,j,t_increment)+dt*dCdt_ben(i,j,t_increment);
                else
                    
                    dCdt_ben(i,j,t_increment)=0;
                
                    C_ben(i,j,t_increment+1)=0;
                end
            end
                
        else
            for i=2:1:m-1
                    %advection
                    if C_ben(i,j,t_increment)>0
       
                        dCdx_ben(i,j,t_increment+1) = (C_ben(i,j+1,t_increment)-C_ben(i,j-1,t_increment))/(2*dx);

                        %Dispersion
                        d2Cdx2_ben(i,j,t_increment) = (C_ben(i,j+1,t_increment)-(2*C_ben(i,j,t_increment))+C_ben(i,j-1,t_increment))/(dx^2);
                        d2Cdy2_ben(i,j,t_increment) = (C_ben(i+1,j,t_increment)-(2*C_ben(i,j,t_increment))+C_ben(i-1,j,t_increment))/(dy^2);

                        dCdt_ben(i,j,t_increment)=-seepage_Q*dCdx_ben(i,j,t_increment)+Dx*d2Cdx2_ben(i,j,t_increment)+Dy*d2Cdy2_ben(i,j,t_increment);
                        % solving for new Cs using dCdt
                        C_ben(i,j,t_increment+1)=C_ben(i,j,t_increment)+dt*dCdt_ben(i,j,t_increment);
                        
                    else
                        dCdx_ben(i,j,t_increment+1) = 0;

                        %Dispersion
                        d2Cdx2_ben(i,j,t_increment) = 0;
                        d2Cdy2_ben(i,j,t_increment) = 0;

                        dCdt_ben(i,j,t_increment)=0;
                        % solving for new Cs using dCdt
                        C_ben(i,j,t_increment+1)=C_ben(i,j,t_increment)+dt*dCdt_ben(i,j,t_increment);
                    end
                        
            end
        end
            
    end
    
    for j=round(n/2):1:round(n/2)+1
            if j==n-1;
                
                %forced dCdt in last column in each time step to be the
                %same as numbers as a colume before
                for i=2:1:m-1
                

                dCdt_ben(i,j,t_increment)=dCdt_ben(i,j-1,t_increment);
                
                C_ben(i,j,t_increment+1)=C_ben(i,j,t_increment)+dt*dCdt_ben(i,j,t_increment);
                end
                
            else
                for i=round(m/2):1:round(m/2)+1
                    
                    
                
                C_ben_eq=xb*C_ben_sat;
                
                    if C_ben(i,j,t_increment)>= C_ben_eq
                        dC_ben(t_increment)=0;
                   
                        C_ben(i,j,t_increment+1)=C_ben(i,j,t_increment)+dCdt_ben(i,j,t_increment)*dt;
                
                    else
                        dC_ben(t_increment)=K_ben*AV*(C_ben_eq-C_ben(i,j,t_increment));
                        C_ben(i,j,t_increment+1)=C_ben(i,j,t_increment)+dC_ben(t_increment)*dt+dCdt_ben(i,j,t_increment)*dt;
                    end
                
                end
            end

    end

    
    
%TCE
    
    
         C_TCE(:,1,t_increment)=0; %left bound always 0
        C_TCE(1,:,t_increment)=C_TCE(3,:,t_increment);% top row imaginary nodes are mirror image of 3rd row
        C_TCE(m,:,t_increment)=C_TCE(m-2,:,t_increment); %bottom row imaginary nodes are mirror images of ttwo rows above
        C_TCE(:,n,t_increment)=C_TCE(:,n-2,t_increment)+(2*seepage_Q*C_TCE(:,n-1,t_increment)*dx/Dx); %no advction or disperssion flux out passed wells.
        
        %C(well_1,n,t_increment) = ((-seepage_Q*(pi*0.25^2)*C(well_1,n-1,t_increment)+seepage_v*C(well_1,n-1,t_increment))*2*dx/Dx)+C(well_1,n-2,t_increment);
        %C(well_2,n,t_increment) = ((-seepage_Q*(pi*0.25^2)*C(well_2,n-1,t_increment)+seepage_v*C(well_2,n-1,t_increment))*2*dx/Dx)+C(well_2,n-2,t_increment);
       
    
    
       
    for j=2:1:n-1 
            
        if j==n-1;
                
                %forced dCdt in last column in each time step to be the
                %same as numbers as a colume before
            for i=2:1:m-1
  
                if C_TCE(i,j,t_increment)>0
                
                    dCdt_TCE(i,j,t_increment)=dCdt_TCE(i,j-1,t_increment);
                
                    C_TCE(i,j,t_increment+1)=C_TCE(i,j,t_increment)+dt*dCdt_TCE(i,j,t_increment);
                else
                    
                    dCdt_TCE(i,j,t_increment)=0;
                
                    C_TCE(i,j,t_increment+1)=0;
                end
            end
                
        else
            for i=2:1:m-1
                    %advection
                    if C_TCE(i,j,t_increment)>0
       
                        dCdx_TCE(i,j,t_increment+1) = (C_TCE(i,j+1,t_increment)-C_TCE(i,j-1,t_increment))/(2*dx);

                        %Dispersion
                        d2Cdx2_TCE(i,j,t_increment) = (C_TCE(i,j+1,t_increment)-(2*C_TCE(i,j,t_increment))+C_TCE(i,j-1,t_increment))/(dx^2);
                        d2Cdy2_TCE(i,j,t_increment) = (C_TCE(i+1,j,t_increment)-(2*C_TCE(i,j,t_increment))+C_TCE(i-1,j,t_increment))/(dy^2);

                        dCdt_TCE(i,j,t_increment)=-seepage_Q*dCdx_TCE(i,j,t_increment)+Dx*d2Cdx2_TCE(i,j,t_increment)+Dy*d2Cdy2_TCE(i,j,t_increment);
                        % solving for new Cs using dCdt
                        C_TCE(i,j,t_increment+1)=C_TCE(i,j,t_increment)+dt*dCdt_TCE(i,j,t_increment);
                        
                    else
                        dCdx_TCE(i,j,t_increment+1) = 0;

                        %Dispersion
                        d2Cdx2_TCE(i,j,t_increment) = 0;
                        d2Cdy2_TCE(i,j,t_increment) = 0;

                        dCdt_TCE(i,j,t_increment)=0;
                        % solving for new Cs using dCdt
                        C_TCE(i,j,t_increment+1)=C_TCE(i,j,t_increment)+dt*dCdt_TCE(i,j,t_increment);
                    end
            end
        end
            
    end
    
    for j=round(n/2):1:round(n/2)+1
            if j==n-1;
                
                %forced dCdt in last column in each time step to be the
                %same as numbers as a colume before
                for i=2:1:m-1
                

                dCdt_TCE(i,j,t_increment)=dCdt_TCE(i,j-1,t_increment);
                
                C_TCE(i,j,t_increment+1)=C_TCE(i,j,t_increment)+dt*dCdt_TCE(i,j,t_increment);
                end
                
            else
                for i=round(m/2):1:round(m/2)+1
                    
                
                    C_TCE_eq=xt*C_TCE_sat;
                
                    if C_TCE(i,j,t_increment)>= C_TCE_eq
                        dC_TCE(t_increment)=0;
                    
                   
                        C_TCE(i,j,t_increment+1)=C_TCE(i,j,t_increment)+dCdt_TCE(i,j,t_increment)*dt;
                        
                    else
                        dC_TCE(t_increment)=K_TCE*AV*(C_TCE_eq-C_TCE(i,j,t_increment));
                        C_TCE(i,j,t_increment+1)=C_TCE(i,j,t_increment)+dC_TCE(t_increment)*dt+dCdt_TCE(i,j,t_increment)*dt;
                    end
                
                end
            end

    end

            
  
    if mass_ben(t_increment)>0
   
  
        mass_ben(t_increment+1)=-dC_ben(t_increment)*dt*V_trans+mass_ben(t_increment);
    
        xb=(mass_ben(t_increment)/MW_ben)/((mass_ben(t_increment)/MW_ben)+(mass_TCE(t_increment)/MW_TCE));
    
    else
        
        mass_ben(t_increment+1)=0;
        xb=0; 
    end
    
     if mass_TCE(t_increment)>0
   
 
        mass_TCE(t_increment+1)=-dC_TCE(t_increment)*dt*V_trans+mass_TCE(t_increment);
        xt=(mass_TCE(t_increment)/MW_TCE)/((mass_ben(t_increment)/MW_ben)+(mass_TCE(t_increment)/MW_TCE));
    else
        mass_TCE(t_increment+1)=0;
    	xt=0;
        
    end
    
    
    figure(1)
    contourf(C_ben(:,1:n-1,t_increment))
    colorbar;
    drawnow;
    title('Benzene Plume')
    xlabel('distance X direction ') 
    ylabel('distance Y direction') 
     
    
    figure(4)
    contourf(C_TCE(:,1:n-1,t_increment))
    colorbar;
    drawnow;

    title('TCE Plume')
    xlabel('distance X direction ') 
    ylabel('distance Y direction') 
        
    %F(t_increment) = getframe(gcf);
    
        
    Average_C_ben(t_increment)=mean(mean(C_ben(:,1:n-6,t_increment))); %checking average concentration levels     
    Average_C_TCE(t_increment)=mean(mean(C_TCE(:,1:n-6,t_increment)));
    t_increment = t_increment+1;
    
    j=1;
    time(t_increment)=t_increment*dt;
    C_well_ben=(C_ben(:,n-1,t_increment));
    C_well_ave_ben(t_increment)=mean(C_well_ben);  
    
    C_well_TCE=(C_TCE(:,n-1,t_increment));
    C_well_ave_TCE(t_increment)=mean(C_well_TCE);  
end
    
    


%time=t_increment*dt; %days 
%more video writing stuff
%writeVideo(v,F)
%close(v);



figure(4)
plot(time,Average_C_ben)
title('Average benzene Concentration - 2 Side Wells (Numerical)')
xlabel('Time (days)') 
ylabel('Average Concentration (ug/ft3)')

figure(5)
plot(time,Average_C_TCE)
title('Average TCE Concentration - 2 Side Wells (Numerical)')
xlabel('Time (days)') 
ylabel('Average Concentration (ug/ft3)')



