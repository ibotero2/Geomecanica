/*
clear all;
close all;
*/
#include <stdio.h>
#include <math.h>
#include<conio.h>
//#include <iostream>

//using namespace std;

int main()
{
    float partdiam = 100;    // mm
    float part_thickness = 1;     // mm
    float sphdiam = 100;        // mm
    float Vo      = -300;       // m/s
    float dt      = 0.0005;     // time step in seconds
    float Tt      = 2.5;        // simulation time in seconds
    float E   =  200000.0;    // N/mm^2 (200 GPa)
    float nu  =  0.3;      // Poisson ratio
    double rho =  0.0000078;   // kg/mm^3 (7800 kg/m^3)
    float g   =  9.81 ;    // m/s^2. Gravity: (9.81 m/s2)

    float kn = E*part_thickness/(sqrt(3)*(1-nu));
    float ks = E*part_thickness*(1-3*nu)/(sqrt(3)*(1.0-nu*nu)); // Cheng, et al (2009). "New discrete element models for elastoplastic problems"

    float box_w  = 600;    // mm
    float box_d  = 400;    // mm
    // box_w  = 300;    // mm
    // box_d  = 300;    % mm

    float Elev    = box_d + sphdiam/2 + 50;        // mm

// Plot particles' initial state

    int np_x = box_w/partdiam - 1;
    int np_z = box_d/partdiam - 1;

    float G      = E/(2*(1+nu));
    float pi     =atan(1)*4;
    float Tr     = (pi*(partdiam/2)*sqrt(rho/G)*1000/(0.1631*nu+0.8766))*0.2;  // Rayleigh time
    float dt_req =  Tr*0.2;  // 0.2 to 0.4 Tr
//dt
//   printf("%d", np_x );

    float v_x[np_x+2],v_z[np_z+2];

    int i=0;
    for (i=0; i<np_x+1;i++)
        v_x[i]= partdiam/2 + i*partdiam;
        i=i+1;
    int h=0;
    for (h=0; h<np_z+1;h++)
        v_z[h]= partdiam/2 +h*partdiam;
        i=h+1;

  float XX [np_z+1][np_x+1], ZZ[np_z+1][np_x+1], v_coords[(np_x+1)*(np_z+1)][2];

	for (h=0; h<np_z+1;h++){
        for (i=0; i<np_x+1;i++){
            XX[h][i]=partdiam/2 + i*partdiam;
            ZZ[h][i]=partdiam/2 + h*partdiam;
    	}
	}
	int j=0;
	int cnt=0;
	for( i=0; i<np_x+1; i++){
	    for (j=0;j<np_z+1;j++){
	        v_coords[cnt][0]=  XX[j][i];    //Se llena con la malla
	        v_coords[cnt][1]=  ZZ[j][i];   //Se llena con la malla
	        cnt=cnt+1;
		}
	}
	v_coords[cnt][0]=box_w/2;//hitting disc
	v_coords[cnt][1]=  Elev;


    //marching time
	float NTimeSteps= 5;
	//NTimeSteps = floor(Tt/dt);
	int Npart = cnt;
	float ro = partdiam/2;
	float coords_o [cnt][2];
	float coords_act [cnt][2];

	for (i=0;i<cnt;i++){
		for (j=0; j<2;j++){
			coords_o[i][j]=v_coords[i][j];
			coords_act[i][j]= coords_o [i][j];
		}
	}
	float mpart = rho*part_thickness*pi*(ro*ro);  // kg
	float Ipart = (0.5)*mpart*(ro*ro);
	float M_part [3][3]= {{1*mpart, 0, 0},{0, 1*mpart, 0},{ 0,0, mpart*0.5*ro*ro}};


	float	Ut_m1[Npart][3], Ut_p1[Npart][3], Ut[Npart][3], Vt_mhalf[Npart][3], Vt_plushalf[Npart][3];
    float NormFces_t[Npart][Npart], TangFces_t[Npart][Npart];
	//displacements at t-1, displacements at t+1, current displacements, velocities at t=t-dt/2, velocities at t=t+dt/2

	Vt_mhalf[Npart][2] = Vo;      //Hasta aca va bien no hemos empezado contactos y fuerzas
	float se_part[cnt];  //varia de tamaño ¡¡¡¡
	int ww=0;
	for (i=0;i<cnt;i++){
		if (coords_act[i][1]==ro);   //No se si esta bien
			se_part[ww]=i;
			ww=ww+1;
}
    //se_part = find(coords_act(:,2)==ro);
	float cc_part[2], uu_part [2], vv_part[2];

	for (i=0;i<NTimeSteps;i++){
	    for (j=0;j<Npart;j++){
	        cc_part [0]   = coords_act[j][0];
	        cc_part [1]   = coords_act[j][1];
	        uu_part [0]   = Ut[j][0];
	        uu_part [1]   = Ut[j][1];
	        vv_part [0]   = Vt_mhalf[j][0];
	        vv_part [1]   = Vt_mhalf[j][1];
	        //      uu_part_m1 = Ut_m1(j,:);

	        if (i==38 && j==Npart);
	            int po=90;
		float contacts [cnt]; 
		float cont_part [cnt];  //¡¡¡Varia de tamaño
	    // find contacts
	
        // get contact forces Entra funcion que no sabemos   Nicolas lo esta haciendo
        int aa=0;
        for (aa=0;aa<cnt;aa++) {
			//if (( (coords_act[aa][0]) <= (cc_part[0]+partdiam*2) )*(( coords_act[aa][1] <= (cc_part[1]+partdiam*2) ))*(( coords_act[aa][1] <= (cc_part[1]+partdiam*2) ))*( coords_act[aa][1] >= (cc_part[1]-partdiam*2) ))
            //    contacts [aa] =1;
            //else
            //    contacts [aa] =0;
            contacts[aa] = (( (coords_act[aa][0]) <= (cc_part[0]+partdiam*2) )*(( coords_act[aa][1] <= (cc_part[1]+partdiam*2) ))*(( coords_act[aa][1] <= (cc_part[1]+partdiam*2) ))*( coords_act[aa][1] >= (cc_part[1]-partdiam*2) ));
			int jj=0;
			if (contacts[aa] ==1);
				cont_part[jj]=aa;
				jj=jj+1;
	}
		

		int mm=3;
		float CntFzes[3][mm];
        //[ CntFzes, Normal_up, Tang_up ] = contForces2(cc_part, uu_part, vv_part, coords_act(cont_part,:), Ut(cont_part,:), Vt_mhalf(cont_part,:), NormFces_t(:,j), TangFces_t(:,j), ks*1000, kn*1000, ro, dt, cont_part);
//Funcion junajo
        //NormFces_t(:,j) = Normal_up;  //No se usan 
        //TangFces_t(:,j) = Tang_up;
		/*
        // obtain upart at t+1
		if (sizeof(cont_part)==1);
            float vFzes[3][mm];   //No conocemos el tamaño
			vFzes [0][mm]= CntFzes[0][mm];
			vFzes [1][mm]= CntFzes[1][mm];
			vFzes [2][mm]= CntFzes[2][mm];
        else
        	float	vFzes [3][mm];*/
        	
            //vFzes = sum(CntFzes')' ;
          //  mm_j es mm de funcion que no conocemos
       /*     for(i=0;mm_j<;i++){                     No se si el for va dentro de los otros o no ¿?
            	vFzes[0][i]= vFzes[0][i]+ CntFzes[0][i];
            	vFzes[1][i]= vFzes[1][i]+ CntFzes[1][i];
            	vFzes[2][i]= vFzes[2][i]+ CntFzes[2][i];
        	}*/

        //vFzes[1] = vFzes[1]- mpart*g;                                    //add gravity. Forces measured in N. Negative downwards

    /*    for (i=0;i<cnt;i++){
	    	for (j=0;j<3;j++){
	    		  //      Vt_plushalf(j,:) =  Vt_mhalf(j,:) + dt * (inv(M_part) * vFzes)' ; // These are in m/s  Como hacer la inversa en c
				Ut_p1[j][0] =  Ut[j][0] + (dt * Vt_plushalf[j][0]);                  // These are in meters
				Ut_p1[j][1] =  Ut[j][1] + (dt * Vt_plushalf[j][1]);                  // These are in meters
				Ut_p1[j][2] =  Ut[j][2] + (dt * Vt_plushalf[j][2]);                  // These are in meters
			}
		}/*

/*		for (i=0;i<cnt;i++){
	    	for (j=0;j<2;j++){
    			coords_act[i][0] = coords_o[i][0] + Ut_p1(:,1:2)*1000;    // que son esos puntos ??
    			coords_act[i][1] = coords_o[i][0] + Ut_p1(:,1:2)*1000;    // que son esos puntos ??
		}
		}*/
/*
    	for (mm=0;mm<Npart;mm++){

        // check for base particles
//  %        if ( coords_act(mm,2) < ro && Vt_plushalf(mm,2)<0)
        	if ( coords_act[mm][1] < ro );
           	 	Vt_plushalf[mm][1]=0;
            	coords_act[mm][1] = ro;
  		}
*/
{  //  % % % %     coords_act = coords_o + Ut_p1(:,1:2)*1000;    //Comentarios Doriam
   // % % % %
  //  % % % %     for mm=1:Npart
  //  % % % %         % check for base particles
  //  % % % %         if ( coords_act(mm,2) < ro && Vt_plushalf(mm,2)<0)
  //  % % % % %         if ( coords_act(mm,2) < ro )
  //  % % % %             Vt_plushalf(mm,2)=0;
  //  % % % %             coords_act(mm,2) = ro;
  //  % % % %         end
  //  % % % %     end
}

    //update coords and displ
    
    	for (i=0;i<cnt;i++){
	    	for (j=0;j<3;j++){
   	 			Ut [i][j]   = Ut_p1[i][j];
    			Vt_mhalf [i][j]= Vt_plushalf[i][j];			
		}
	}
{    //     pause(0.00001)  // Graficas ¿Donde las vamos a hacer?????????
    //figure(1000)
    //aa=coords_act;
    //for kk=1:Npart
    //    if kk==Npart
    //        circle2d(aa(kk,1),aa(kk,2),partdiam/2,'m', '-');
    //    else
    //        circle2d(aa(kk,1),aa(kk,2),partdiam/2,'b', '-');
    //    end
    //    hold on
    //end
//%     xlim([0 box_w])
//    ylim([0 1.5*(Elev)])
//    xlim([(-1.5*Elev/2 + box_w/2) (1.5*Elev/2 + box_w/2) ])
//    grid on
//axis square
 //   hold off
}

   // int po=90;
	return 0;
	}
}
    
}


