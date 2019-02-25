#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    r->dt                   = 0.01*2.*M_PI;                // initial timestep
    r->integrator           = REB_INTEGRATOR_IAS15;
    r->collision            = REB_COLLISION_DIRECT;
    r->collision_resolve    = reb_collision_resolve_merge;        // Choose merger collision routine.
    r->heartbeat            = heartbeat;

    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.1;
    reb_add(r, star);
    
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
    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

  //  reb_integrate(r, INFINITY);
    reb_integrate(r, 30000);

}

void heartbeat(struct reb_simulation* r){
        int snap_n = 0;
        char SNAP[30];
    if (reb_output_check(r, 10.*2.*M_PI)){  
        reb_output_timing(r, 0);
       sprintf(SNAP,"output/snap01/snap%05d.dat",snap_n);
       reb_output_orbits(r,SNAP);
   //    reb_output_orbits(r,"output/snap01/snap.dat");
           snap_n ++;
    }
}

