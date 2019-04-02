clc
clear all;
close all;

partdiam       = 50;    % mm
part_thickness = 3;     % mm

sphdiam = 50;        % mm
Vo      = -300;       % m/s

dt      = 0.00025;     % time step in seconds
Tt      = 2.5;        % simulation time in seconds

E   =  10.0E1;    % N/mm^2 (200 GPa)
nu  =  0.3;      % Poisson ratio
rho =  78E-07;   % kg/mm^3 (7800 kg/m^3)
g   =  9.81 ;    % m/s^2. Gravity: (9.81 m/s2)


kn = E*part_thickness/(sqrt(3)*(1-nu));
ks = E*part_thickness*(1-3*nu)/(sqrt(3)*(1-nu^2)); % Cheng, et al (2009). "New discrete element models for elastoplastic problems"

box_w  = 1000;    % mm
box_d  = 1500;    % mm
% % box_w  = 300;    % mm
% % box_d  = 300;    % mm


Elev    = box_d + sphdiam/2 + 50;        % mm

%% Plot particles' initial state
np_x = box_w/partdiam - 1;
np_z = box_d/partdiam - 1;

G      = E/(2*(1+nu));
Tr     = (pi*(partdiam/2)*sqrt(rho/G)*1000/(0.1631*nu+0.8766))*0.2;  % Rayleigh time
dt_req =  Tr*0.2  % 0.2 to 0.4 Tr
dt

v_coords=zeros((np_x+1)*(np_z+1),2);

v_x= partdiam/2:partdiam:box_w;
v_z= partdiam/2:partdiam:box_d;

[XX ZZ] = meshgrid(v_x,v_z);

figure(1000)
hold on


cnt=1;
for i=1:np_x+1

    for j=1:np_z+1

        v_coords(cnt,1)=  XX(j,i);
        v_coords(cnt,2)=  ZZ(j,i);
        cnt=cnt+1;

        circle2d(XX(j,i),ZZ(j,i),partdiam/2,'b', '-');

    end

end

v_coords(cnt,1)=  box_w/2; % hitting disc
v_coords(cnt,2)=  Elev;
circle2d(box_w/2,Elev,sphdiam/2,'r', '-');
    ylim([0 1.5*(Elev)])
    xlim([(-1.5*Elev/2 + box_w/2) (1.5*Elev/2 + box_w/2) ])
grid on
hold off
axis square

% marching time
NTimeSteps = floor(Tt/dt);
Npart = cnt;
ro = partdiam/2;

coords_o      = v_coords;
coords_act    = coords_o;
Ut_m1         = zeros(Npart,3); % displacements at t-1
Ut_p1         = zeros(Npart,3); % displacements at t+1
Ut            = zeros(Npart,3); % current displacements
Vt_mhalf      = zeros(Npart,3); % velocities at t=t-dt/2
Vt_plushalf   = zeros(Npart,3); % velocities at t=t+dt/2

NormFces_t = zeros(Npart,Npart);
TangFces_t = zeros(Npart,Npart);


mpart = rho*part_thickness*pi*ro^2;  % kg
Ipart = 1/2*mpart*ro^2;
M_part= mpart*[1 0 0;0 1 0;0 0 1/2*ro^2];

Vt_mhalf(Npart,2) = Vo;
% % base_part = find(coords_act(:,2)==ro);
% for i=1: NTimeSteps
tic
for i=1:NTimeSteps

    for j=1:Npart

        cc_part    = coords_act(j,:);
        uu_part    = Ut(j,:);
        vv_part    = Vt_mhalf(j,:);
        % %         uu_part_m1 = Ut_m1(j,:);

        if i==38 && j==Npart
            po=90;
        end

        % find contacts
        contacts  = ( coords_act(:,1) <= (cc_part(1)+partdiam) ).*( coords_act(:,1) >= (cc_part(1)-partdiam) ).*( coords_act(:,2) <= (cc_part(2)+partdiam) ).*( coords_act(:,2) >= (cc_part(2)-partdiam) );
        cont_part = find(contacts==1);
        % get contact forces

        [ CntFzes, Normal_up, Tang_up ] = contForces2(cc_part, uu_part, vv_part, coords_act(cont_part,:), Ut(cont_part,:), Vt_mhalf(cont_part,:), NormFces_t(:,j), TangFces_t(:,j), ks*1000, kn*1000, ro, dt, cont_part);

        NormFces_t(:,j) = Normal_up;
        TangFces_t(:,j) = Tang_up;

        % obtain upart at t+1
        if length(cont_part)==1
            vFzes = CntFzes;
        else
            vFzes = sum(CntFzes')' ;
        end

        vFzes(2) = vFzes(2) - mpart*g;                                    % add gravity. Forces measured in N. Negative downwards
        Vt_plushalf(j,:) =  Vt_mhalf(j,:) + dt * (inv(M_part) * vFzes)' ; % These are in m/s
        Ut_p1(j,:) =  Ut(j,:) + (dt * Vt_plushalf(j,:));                  % These are in meters

    end

    coords_act = coords_o + Ut_p1(:,1:2)*1000;

    for mm=1:Npart
        % check for base particles
%         if ( coords_act(mm,2) < ro && Vt_plushalf(mm,2)<0)
           if ( coords_act(mm,2) < ro )
            Vt_plushalf(mm,2)=0;
            coords_act(mm,2) = ro;
        end
    end

    % % % %     coords_act = coords_o + Ut_p1(:,1:2)*1000;
    % % % %
    % % % %     for mm=1:Npart
    % % % %         % check for base particles
    % % % %         if ( coords_act(mm,2) < ro && Vt_plushalf(mm,2)<0)
    % % % % %         if ( coords_act(mm,2) < ro )
    % % % %             Vt_plushalf(mm,2)=0;
    % % % %             coords_act(mm,2) = ro;
    % % % %         end
    % % % %     end


    %update coords and displ
    Ut    = Ut_p1;
    Vt_mhalf = Vt_plushalf;

    %      pause(0.00001)
    figure(1000)
    aa=coords_act;
    for kk=1:Npart
        if kk==Npart
            circle2d(aa(kk,1),aa(kk,2),partdiam/2,'m', '-');
        else
            circle2d(aa(kk,1),aa(kk,2),partdiam/2,'b', '-');
        end
        hold on
    end
%     xlim([0 box_w])
    ylim([0 1.5*(Elev)])
    xlim([(-1.5*Elev/2 + box_w/2) (1.5*Elev/2 + box_w/2) ])
    grid on
axis square
    hold off


    po=90;

end
toc
%}
