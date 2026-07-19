#include "parcels.h"

#include "math.h"

const int ngrid = 2;

typedef struct
{
  double *lon;
  double *lon_nextloop;
  double *lat;
  double *lat_nextloop;
  double *depth;
  double *depth_nextloop;
  double *time;
  double *time_nextloop;
  long *id;
  double *dt;
  int *obs_written;
  int *state;
  int *ngrids;
  int *xi;
  int *yi;
  int *zi;
  int *ti;
}  JITParticlep;

static inline StatusCode AdvectionRK4(JITParticlep *particles, int pnum, double time, CField *U1, CField *V3, CField *U2, CField *V4)
{
  float parcels_tmpvar7 = 0;
  float parcels_tmpvar6 = 0;
  float parcels_tmpvar5 = 0;
  float parcels_tmpvar4 = 0;
  float parcels_tmpvar3 = 0;
  float parcels_tmpvar2 = 0;
  float parcels_tmpvar1 = 0;
  float parcels_tmpvar0 = 0;
  type_coord v4 = 0;
  type_coord u4 = 0;
  type_coord lat3 = 0;
  type_coord lon3 = 0;
  type_coord v3 = 0;
  type_coord u3 = 0;
  type_coord lat2 = 0;
  type_coord lon2 = 0;
  type_coord v2 = 0;
  type_coord u2 = 0;
  type_coord lat1 = 0;
  type_coord lon1 = 0;
  type_coord v1 = 0;
  type_coord u1 = 0;
  int parcels_interp_state;
  type_coord particle_dlon = 0;
  particles->lon[pnum] = particles->lon_nextloop[pnum];
  type_coord particle_dlat = 0;
  particles->lat[pnum] = particles->lat_nextloop[pnum];
  type_coord particle_ddepth = 0;
  particles->depth[pnum] = particles->depth_nextloop[pnum];
  particles->time[pnum] = particles->time_nextloop[pnum];
  ;
  while (1==1)
  {
    particles->state[pnum] = temporal_interpolationUV(particles->lon[pnum], particles->lat[pnum], particles->depth[pnum], time, U1, V3, &particles->xi[pnum*ngrid], &particles->yi[pnum*ngrid], &particles->zi[pnum*ngrid], &particles->ti[pnum*ngrid], &parcels_tmpvar0, &parcels_tmpvar1, CGRID_VELOCITY, NEMO);
    {
    }
    if (particles->state[pnum] != ERROROUTOFBOUNDS )
    {
      CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
      break;
    }
    particles->state[pnum] = temporal_interpolationUV(particles->lon[pnum], particles->lat[pnum], particles->depth[pnum], time, U2, V4, &particles->xi[pnum*ngrid], &particles->yi[pnum*ngrid], &particles->zi[pnum*ngrid], &particles->ti[pnum*ngrid], &parcels_tmpvar0, &parcels_tmpvar1, CGRID_VELOCITY, NEMO);
    {
    }
    if (particles->state[pnum] != ERROROUTOFBOUNDS )
    {
      CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
      break;
    }
    CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
    break;
  }
  u1 = parcels_tmpvar0;
  v1 = parcels_tmpvar1;
  lon1 = (particles->lon[pnum] + ((u1 * 0.5) * particles->dt[pnum]));
  lat1 = (particles->lat[pnum] + ((v1 * 0.5) * particles->dt[pnum]));
  while (1==1)
  {
    particles->state[pnum] = temporal_interpolationUV(lon1, lat1, particles->depth[pnum], (time + (0.5 * particles->dt[pnum])), U1, V3, &particles->xi[pnum*ngrid], &particles->yi[pnum*ngrid], &particles->zi[pnum*ngrid], &particles->ti[pnum*ngrid], &parcels_tmpvar2, &parcels_tmpvar3, CGRID_VELOCITY, NEMO);
    {
    }
    if (particles->state[pnum] != ERROROUTOFBOUNDS )
    {
      CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
      break;
    }
    particles->state[pnum] = temporal_interpolationUV(lon1, lat1, particles->depth[pnum], (time + (0.5 * particles->dt[pnum])), U2, V4, &particles->xi[pnum*ngrid], &particles->yi[pnum*ngrid], &particles->zi[pnum*ngrid], &particles->ti[pnum*ngrid], &parcels_tmpvar2, &parcels_tmpvar3, CGRID_VELOCITY, NEMO);
    {
    }
    if (particles->state[pnum] != ERROROUTOFBOUNDS )
    {
      CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
      break;
    }
    CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
    break;
  }
  u2 = parcels_tmpvar2;
  v2 = parcels_tmpvar3;
  lon2 = (particles->lon[pnum] + ((u2 * 0.5) * particles->dt[pnum]));
  lat2 = (particles->lat[pnum] + ((v2 * 0.5) * particles->dt[pnum]));
  while (1==1)
  {
    particles->state[pnum] = temporal_interpolationUV(lon2, lat2, particles->depth[pnum], (time + (0.5 * particles->dt[pnum])), U1, V3, &particles->xi[pnum*ngrid], &particles->yi[pnum*ngrid], &particles->zi[pnum*ngrid], &particles->ti[pnum*ngrid], &parcels_tmpvar4, &parcels_tmpvar5, CGRID_VELOCITY, NEMO);
    {
    }
    if (particles->state[pnum] != ERROROUTOFBOUNDS )
    {
      CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
      break;
    }
    particles->state[pnum] = temporal_interpolationUV(lon2, lat2, particles->depth[pnum], (time + (0.5 * particles->dt[pnum])), U2, V4, &particles->xi[pnum*ngrid], &particles->yi[pnum*ngrid], &particles->zi[pnum*ngrid], &particles->ti[pnum*ngrid], &parcels_tmpvar4, &parcels_tmpvar5, CGRID_VELOCITY, NEMO);
    {
    }
    if (particles->state[pnum] != ERROROUTOFBOUNDS )
    {
      CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
      break;
    }
    CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
    break;
  }
  u3 = parcels_tmpvar4;
  v3 = parcels_tmpvar5;
  lon3 = (particles->lon[pnum] + (u3 * particles->dt[pnum]));
  lat3 = (particles->lat[pnum] + (v3 * particles->dt[pnum]));
  while (1==1)
  {
    particles->state[pnum] = temporal_interpolationUV(lon3, lat3, particles->depth[pnum], (time + particles->dt[pnum]), U1, V3, &particles->xi[pnum*ngrid], &particles->yi[pnum*ngrid], &particles->zi[pnum*ngrid], &particles->ti[pnum*ngrid], &parcels_tmpvar6, &parcels_tmpvar7, CGRID_VELOCITY, NEMO);
    {
    }
    if (particles->state[pnum] != ERROROUTOFBOUNDS )
    {
      CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
      break;
    }
    particles->state[pnum] = temporal_interpolationUV(lon3, lat3, particles->depth[pnum], (time + particles->dt[pnum]), U2, V4, &particles->xi[pnum*ngrid], &particles->yi[pnum*ngrid], &particles->zi[pnum*ngrid], &particles->ti[pnum*ngrid], &parcels_tmpvar6, &parcels_tmpvar7, CGRID_VELOCITY, NEMO);
    {
    }
    if (particles->state[pnum] != ERROROUTOFBOUNDS )
    {
      CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
      break;
    }
    CHECKSTATUS_KERNELLOOP(particles->state[pnum]);
    break;
  }
  u4 = parcels_tmpvar6;
  v4 = parcels_tmpvar7;
  particle_dlon += (((((u1 + (2 * u2)) + (2 * u3)) + u4) / 6.0) * particles->dt[pnum]);
  particle_dlat += (((((v1 + (2 * v2)) + (2 * v3)) + v4) / 6.0) * particles->dt[pnum]);
  particles->lon_nextloop[pnum] = particles->lon[pnum] + particle_dlon;
  particles->lat_nextloop[pnum] = particles->lat[pnum] + particle_dlat;
  particles->depth_nextloop[pnum] = particles->depth[pnum] + particle_ddepth;
  particles->time_nextloop[pnum] = particles->time[pnum] + particles->dt[pnum];
  return particles->state[pnum];
}

void particle_loop(int num_particles, JITParticlep *particles, double endtime, double dt, CField *U1, CField *V3, CField *U2, CField *V4)
{
  int pnum;
  double sign_dt;
  sign_dt = dt > 0 ? 1 : -1;
  for (pnum = 0; pnum < num_particles; ++pnum)
  {
    while ((particles->state[pnum] == EVALUATE || particles->state[pnum] == REPEAT))
    {
      double pre_dt;
      pre_dt = particles->dt[pnum];
      if (sign_dt*particles->time_nextloop[pnum] >= sign_dt*(endtime))
        break;
      if (fabs(endtime - particles->time_nextloop[pnum]) < fabs(particles->dt[pnum])-1e-6)
        particles->dt[pnum] = fabs(endtime - particles->time_nextloop[pnum]) * sign_dt;
      particles->state[pnum] = AdvectionRK4(particles, pnum, particles->time_nextloop[pnum], U1, V3, U2, V4);
      if (particles->state[pnum] == SUCCESS)
      {
        if (sign_dt*particles->time[pnum] < sign_dt*endtime)
        {
          particles->state[pnum] = EVALUATE;
        }
        else
        {
          particles->state[pnum] = SUCCESS;
        }
      }
      if (particles->state[pnum] == STOPALLEXECUTION)
        return;
      particles->dt[pnum] = pre_dt;
      if ((particles->state[pnum] == REPEAT || particles->state[pnum] == DELETE))
      {
        break;
      }
    }
  }
}