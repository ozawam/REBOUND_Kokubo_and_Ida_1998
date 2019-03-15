#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    r->dt                   = 0.001*2.*M_PI;                // initial timestep
 //   r->integrator           = REB_INTEGRATOR_IAS15;
 //   r->integrator    = REB_INTEGRATOR_MERCURIUS;
    r->integrator           = REB_INTEGRATOR_LEAPFROG;
    r->collision            = REB_COLLISION_DIRECT;
 //   r->collision    = REB_COLLISION_TREE;
    r->collision_resolve    = reb_collision_resolve_merge;        // Choose merger collision routine.
    r->heartbeat            = heartbeat;
    r->hash_ctr             = 0;

 //   double boxsize = 3.0;
 //   reb_configure_box(r,boxsize,1,1,1);

   // Star
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.1;
    reb_add(r, star);
  /*
    // Add planets
    int N_planets = 7;
    for (int i=0;i<N_planets;i++){
        double a = 1.+(double)i/(double)(N_planets-1);        // semi major axis in AU
        double v = sqrt(1./a);                     // velocity (circular orbit)
        struct reb_particle planet = {0};
        planet.m = 1e-4; 
        planet.r = 4e-2;                     // radius in AU (it is unphysically large in this example)
        planet.lastcollision = 0;                // The first time particles can collide with each other
        planet.x = a; 
        planet.vy = v;
        reb_add(r, planet); 
    }

*/
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
  
        double radius_f = 6.0;//fold enlargement(radius enhancement)

        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, e, inc, Omega, apsis, phi);
        pt.m = planetesimal_mass;
        pt.r = pow((3.0*planetesimal_mass*M_sun)/(4.0*M_PI*rho) ,1.0/3.0)/ one_au * radius_f; //unit is AU.
        pt.lastcollision = 0;
        reb_add(r, pt);
    }

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

  //  reb_integrate(r, INFINITY);
    reb_integrate(r, 30000);

}

void heartbeat(struct reb_simulation* r){
     //   int snap_n = (int)(r->t)/(100.0*2.0*M_PI);
        int snap_n = (r->hash_ctr);
        char SNAP[30];
    if (reb_output_check(r, 1.0*2.0*M_PI)){  
        reb_output_timing(r, 30000);
           }

    if (reb_output_check(r, 10.0*2.0*M_PI)){
       sprintf(SNAP,"output/snap01/snap%05d.dat",snap_n);
       reb_output_orbits(r,SNAP);
   //    reb_output_orbits(r,"output/snap01/snap.dat");
       r->hash_ctr ++;
    }
}

