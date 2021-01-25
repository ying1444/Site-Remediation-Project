clear


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

Dy=2; %ft/day %dispersion coefficient in Y
Dx=(Dy*(X_plume_edge/3)^2/(Y_plume_edge/3)^2); %ft/day %dispersion coefficient in X


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

dx=5; 

dy=5;
dt=2; %should be somewhat less than dx/seepage_Q for stability of code

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
C(:,:,1)=C; % first time step will take on the initial plume concentrations
dCdt=zeros(m,n,t_increment); %initializing dCdtt matrix with zeros because special midification required

% next block of code is for frame capturing and exporting of videos. ignore
% until results are good.
set(gca, 'nextplot', 'replacechildren');  
caxis manual; 
caxis([0 7e5]);
%v = VideoWriter('10 year no source.avi');
%v.FrameRate=30;
%open(v);
%F(t_increment) = getframe(gcf);


while   any(any(C(:,1:n-1,t_increment)>10*28.3))>0 %while concentrations anywhere are above 10 ug/L keep running through loops
        
     
        C(:,1,t_increment)=0; %left bound always 0
        C(1,:,t_increment)=C(3,:,t_increment);% top row imaginary nodes are mirror image of 3rd row
        C(m,:,t_increment)=C(m-2,:,t_increment); %bottom row imaginary nodes are mirror images of ttwo rows above
        C(:,n,t_increment)=C(:,n-2,t_increment)+(2*seepage_Q*C(:,n-1,t_increment)*dx/Dx); %no advction or disperssion flux out passed wells.
        
        %C(well_1,n,t_increment) = ((-seepage_Q*(pi*0.25^2)*C(well_1,n-1,t_increment)+seepage_v*C(well_1,n-1,t_increment))*2*dx/Dx)+C(well_1,n-2,t_increment);
        %C(well_2,n,t_increment) = ((-seepage_Q*(pi*0.25^2)*C(well_2,n-1,t_increment)+seepage_v*C(well_2,n-1,t_increment))*2*dx/Dx)+C(well_2,n-2,t_increment);
       
    
    
       
    for j=2:1:n-1 
            
        if j==n-1;
                
                %forced dCdt in last column in each time step to be the
                %same as numbers as a colume before
            for i=2:1:m-1
                

                dCdt(i,j,t_increment)=dCdt(i,j-1,t_increment);
                
                C(i,j,t_increment+1)=C(i,j,t_increment)+dt*dCdt(i,j,t_increment);
            end
                
        else
            for i=2:1:m-1
                    %advection
                dCdx(i,j,t_increment) = (C(i,j+1,t_increment)-C(i,j-1,t_increment))/(2*dx);

                %Dispersion
                d2Cdx2(i,j,t_increment) = (C(i,j+1,t_increment)-(2*C(i,j,t_increment))+C(i,j-1,t_increment))/(dx^2);
                d2Cdy2(i,j,t_increment) = (C(i+1,j,t_increment)-(2*C(i,j,t_increment))+C(i-1,j,t_increment))/(dy^2);

                dCdt(i,j,t_increment)=-seepage_Q*dCdx(i,j,t_increment)+Dx*d2Cdx2(i,j,t_increment)+Dy*d2Cdy2(i,j,t_increment);
                % solving for new Cs using dCdt
                C(i,j,t_increment+1)=C(i,j,t_increment)+dt*dCdt(i,j,t_increment);
            end
        end
            
    end
    figure(1)
    contourf(C(:,1:n-1,t_increment))
    colorbar;
    %drawnow;

    title('Removed Source 10 Year Plume')
    xlabel('distance X direction ') 
    ylabel('distance Y direction') 
        
    %F(t_increment) = getframe(gcf);
    
        
    Average_C(t_increment)=mean(mean(C(:,1:n-6,t_increment))); %checking average concentration levels     
        
    t_increment = t_increment+1;
    
    j=1;
    time(t_increment)=t_increment*dt;
    C_well=(C(:,n-1,t_increment));
    C_well_ave(t_increment)=mean(C_well);    
end
    
%time=t_increment*dt; %days 
%more video writing stuff
%writeVideo(v,F)
%close(v);

figure(2)
plot(Average_C);        
title('Average Concentraton - 2 Side Wells (Numerical)')
xlabel('Time (days)') 
ylabel('Average Concentration at in Capture Zone (ug/ft3)');


figure(3)
plot(time,C_well_ave)
title('Well Concentration - 2 Side Wells (Numerical)')
xlabel('Time (days)') 
ylabel('Average Concentration at well (ug/ft3)') 



