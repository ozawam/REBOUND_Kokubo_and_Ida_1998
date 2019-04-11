#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"

void additional_forces(struct reb_simulation* const r);
void heartbeat(struct reb_simulation* r);

double tmax = 300000.0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();

//Setup constants.
    r->dt                   = 1e-4;//0.001*2.*M_PI;                // initial timestep

// Setup integrator.
    r->integrator           = REB_INTEGRATOR_IAS15;
//   r->integrator    = REB_INTEGRATOR_MERCURIUS;
//    r->integrator           = REB_INTEGRATOR_LEAPFROG;

// Setup collision.
    r->collision            = REB_COLLISION_DIRECT;
 //   r->collision    = REB_COLLISION_TREE;
    r->collision_resolve    = reb_collision_resolve_merge;        // Choose merger collision routine.

// Setup callback function for velocity dependent forces.
    r->additional_forces     = additional_forces;    r->force_is_velocity_dependent = 1;

// Setup callback function for outputs.
    r->heartbeat            = heartbeat;
    r->hash_ctr             = 0;
    r->usleep        = 10000;        // Slow down integration (for visualization only)

 //   double boxsize = 3.0;
 //   reb_configure_box(r,boxsize,1,1,1);


   // Star
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.1;
    reb_add(r, star);

   // Planetesimal disk parameters
    double M_sun=1.989*pow(10.0,33.0);//g
    double total_disk_mass = 4000.0* 3.0*pow(10.0,23.0)/M_sun;
    int N_planetesimals = 1000;
    double planetesimal_mass = total_disk_mass/N_planetesimals;
    double delta_a = 0.085;
    double central_a = 1.0;
    double amin = central_a - (delta_a/2.0), amax = central_a + (delta_a/2.0);   //planet at inner edge of disk
    double powerlaw = -1.5;
    
    r->N_active = N_planetesimals + 10;

    // Generate Planetesimal Disk
     for (int i=0;i< N_planetesimals;i++){
        struct reb_particle pt = {0};
        double a    = reb_random_powerlaw(amin,amax,powerlaw);
        double e    = reb_random_rayleigh(0.0014);
        double inc  = reb_random_rayleigh(0.0007);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        double phi     = reb_random_uniform(0,2.*M_PI);

        double rho = 2.0;
        double one_au = 1.49597870*pow(10.0,13.0);//cm
  
        double radius_f = 1.0;//fold enlargement(radius enhancement)

        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, e, inc, Omega, apsis, phi);
        pt.m = planetesimal_mass;
        pt.r = pow((3.0*planetesimal_mass*M_sun)/(4.0*M_PI*rho) ,1.0/3.0)/ one_au * radius_f; //unit is AU.
        pt.lastcollision = 0;
        reb_add(r, pt);
    }

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    reb_integrate(r, INFINITY);
  //  reb_integrate(r, tmax);

}

void additional_forces(struct reb_simulation* const r){
    // velocity dependent drag force.
	struct reb_particle* const particles = r->particles;
	//Setup constants. Units AU,Msun,T G=1 
	double _rho_planet = 2.0;//g/cm^3
	double rho_planet = _rho_planet * 1.49597871e13 * 1.49597871e13 * 1.49597871e13 / 1.9884e33; // [Msun/AU^3]
	double mass_sun = particles[0].m;
	double lum_sun = 1.0;
	//posとvelはstar. を代用すること
	double PI = 4.0*atan(1.0);
	double alpha = 11.0 / 4.0; // exponent of rho
        double beta = 0.5; // exponent of temperature
	double Cd = 1.0; // coeficient aero drag
	double f_gas = 1.0; // scaling factor of gas density
	double tau_gas = 1e6 * 2.0 * PI; // 1e6[yr] * 2PI[T]/[yr]

	//Calculation of Aero Drag
	//coef
	double coef_rho_gas = exp(- r->t /tau_gas) * 2400.0 * f_gas * 1.49597871e13 * 1.49597871e13 / (0.047*2 * 1.9884e33) * pow(lum_sun, -0.125) * pow(mass_sun, 0.5);
	double coef_cs_vk = 1.0/29.78 * pow(lum_sun, 0.125) * pow(mass_sun, -0.5); // about 0.0018 for solar luminosity and solar mass
	// 1.0/29.78 is cs/vk at 1AU (cs=1km/sec, vkep=29.78km/sec)
	double coef_acc_gd = 0.5*Cd*PI;

	//double dragcoefficient = 1e-10;

	const int N = r->N;
	for (int i=1;i<N;i++){
        //variable
	double r_sq = particles[i].x*particles[i].x + particles[i].y*particles[i].y;
	double inv_r = 1.0 / sqrt(r_sq);
	double r = r_sq * inv_r;
	double r_sqrt = sqrt(r);
	double r_4rt = sqrt(r_sqrt);
	double ev[3] = {-particles[i].y*inv_r, particles[i].x*inv_r, 0.0}; // unit vector of kepler velocity
	double vkep[3] = { sqrt(mass_sun * inv_r) * ev[0],sqrt(mass_sun * inv_r) * ev[1],sqrt(mass_sun * inv_r) * ev[2],}; // kepler velocity //0 x, 1 y, 2 z
	double cs_vk = coef_cs_vk * r_4rt;
	double eta = 0.5*(alpha+beta)*cs_vk*cs_vk;
	double vgas[3] = {(1.0 - eta)*vkep[0],(1.0 - eta)*vkep[1],(1.0 - eta)*vkep[2]};
	double u[3] = {particles[i].vx - vgas[0],particles[i].vx - vgas[1],particles[i].vx - vgas[2]};
	double rplanet = cbrt(3.0*particles[i].m/(4.0*PI*rho_planet));
	double rho_gas = coef_rho_gas * inv_r * inv_r * inv_r * r_4rt;

        double sys_acc_gd[3] = { (-coef_acc_gd * rplanet * rplanet * rho_gas * sqrt(u[0]*u[0]) * u[0]) / particles[i].m, (-coef_acc_gd * rplanet * rplanet * rho_gas * sqrt(u[1]*u[1]) * u[1]) / particles[i].m, (-coef_acc_gd * rplanet * rplanet * rho_gas * sqrt(u[2]*u[2]) * u[2]) / particles[i].m } ;  
	particles[i].ax += sys_acc_gd[0];
        particles[i].ay += sys_acc_gd[1];
        particles[i].az += sys_acc_gd[2];
    }
}

void heartbeat(struct reb_simulation* r){
     //   int snap_n = (int)(r->t)/(100.0*2.0*M_PI);
        int snap_n = (r->hash_ctr);
        char SNAP[30];
    if (reb_output_check(r, 1.0*2.0*M_PI)){  
        reb_output_timing(r, tmax);
           }

    if (reb_output_check(r, 10.0*2.0*M_PI)){
       sprintf(SNAP,"output/snap01/snap%05d.dat",snap_n);
       reb_output_orbits(r,SNAP);
   //    reb_output_orbits(r,"output/snap01/snap.dat");
       r->hash_ctr ++;
    }
}

