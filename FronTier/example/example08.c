/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


/*
*				example0.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a circle in a double vortex field. It demonstrates
*	the resolution of the front tracking method.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static float level_circle_func(POINTER,float*);
static int test_double_vortex_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,float*);

char *in_name,*restart_state_name,*restart_name,*out_name;
bool RestartRun;
int RestartStep;
bool binary = YES;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
	        /* equation for line is x^2/a^2 + y^2/b^2 = 1 */
        float x0,x1;
        float y0,y1;        
        float R; 
        float r;  
} CIRCLE_PARAMS;


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

typedef struct {
	float i1,i2;
        float cen1[2],cen2[2];
} DOUBLE_VORTEX_PARAMS;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	CIRCLE_PARAMS circle_params;	/* level function parameters */
	DOUBLE_VORTEX_PARAMS dv_params; /* velocity function parameters */
	Locstate  sl;

	f_basic.dim = 2;
	FrontInitStandardIO(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 100;	f_basic.gmax[1] = 100;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_name,"%s.ts%s",restart_name,right_flush(RestartStep,7));
#if defined(__MPI__)
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

	FrontStartUp(&front,&f_basic);

	if (!RestartRun)
	{
	    /* Initialize interface through level function */

	    circle_params.x0 = 0.5;
	    circle_params.y0 = 0.2;
	    circle_params.x1 = 0.5;
	    circle_params.y1 = 0.5;
	    circle_params.R = 0.3;
            circle_params.r = 0.3;
	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&circle_params;
	    level_func_pack.func = level_circle_func;

	    FrontInitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	dv_params.cen1[0] = 0.25;
	dv_params.cen1[1] = 0.50;
	dv_params.cen2[0] = 0.75;
	dv_params.cen2[1] = 0.50;
	dv_params.i1 = -0.5;
	dv_params.i2 =  0.5;

	velo_func_pack.func_params = (POINTER)&dv_params;
	velo_func_pack.func = test_double_vortex_vel;
	velo_func_pack.point_propagate = fourth_order_point_propagate;

	FrontInitVelo(&front,&velo_func_pack);

	/* Propagate the front */

	test_propagate(&front);

	clean_up(0);
	return 0;
}

static  void test_propagate(
        Front *front)
{
        int ip,im,status,count;
        Front *newfront;
        float dt,dt_frac,CFL;
        bool is_print_time, is_movie_time, time_limit_reached;
        char s[10];
        float fcrds[MAXD];
        int  dim = front->rect_grid->dim;

	front->max_time = 1.5;
	front->max_step = 10000;
	front->print_time_interval = 0.5;
	front->movie_frame_interval = 0.02;

        CFL = Time_step_factor(front);

        printf("dim = %d\n", dim);
	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun)
	{
            redistribute(front,YES,NO);

            front->time = 0.0;
            front->dt = 0.0;
            front->step = 0;

	    // Always output the initial interface.
	    FrontPrintOut(front,out_name);
            FrontMovieFrame(front,out_name,binary);
            ip = im = 1;

	    // This is a virtual propagation to get maximum front 
	    // speed to determine the first time step.

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            front->dt = CFL*FrontHypTimeStep(front); 
	}
	else
	{
	    ip = (int)(front->time/front->print_time_interval + 1.0);
            im = (int)(front->time/front->movie_frame_interval + 1.0);
	}
	front->dt = FrontOutputTimeControl(front,
			&is_movie_time,&is_print_time,
			&time_limit_reached,&im,&ip);

        for (;;)
        {
	    /* Propagating interface for time step dt */

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            assign_interface_and_free_front(front,newfront);

            ++front->step;
            front->time += front->dt;

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

            front->dt = CFL*FrontHypTimeStep(front); 

	    /* Output section */

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (is_print_time || time_limit_reached)
                print_front_output(front,out_name);
            if (is_movie_time || time_limit_reached)
                show_front_output(front,out_name,binary);

            if (time_limit_reached)
                    break;

	    front->dt = FrontOutputTimeControl(front,
			&is_movie_time,
			&is_print_time,
			&time_limit_reached,
			&im,&ip);
        }
        (void) delete_interface(front->interf);
}       /* end test_propagate */

/********************************************************************
 *	Sample (circle) level function for the initial interface    *
 ********************************************************************/

static float level_circle_func(
        POINTER func_params,
        float *coords)
{
	CIRCLE_PARAMS *circle_params = (CIRCLE_PARAMS*)func_params;
	float x0,y0,x1,y1,R,r,dist1,dist2;

	x0 = circle_params->x0;
	y0 = circle_params->y0;
	x1 = circle_params->x1;
	y1 = circle_params->y1;
	R  = circle_params->R;
        r = circle_params->r;
	dist1 = sqrt(sqr(coords[0] - x0) + sqr(coords[1] - y0)) - R;
	dist2 = sqrt(sqr(coords[0] - x1) + sqr(coords[1] - y1)) - r;
	if(dist1 > 0.0 && dist2 < 0.0) return -min(dist1,-dist2);
	else return min(fabs(dist1),fabs(dist2));

	
}	/* end level_circle_func */


/********************************************************************
 *	Sample (circle) velocity function for the front    *
 ********************************************************************/

static int test_double_vortex_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	float *vel)
{
	DOUBLE_VORTEX_PARAMS *dv_params = (DOUBLE_VORTEX_PARAMS*)params;
	float *coords = Coords(p);
	float d1,d2;
	float s1,s2;
	float *cen1 = dv_params->cen1;
	float *cen2 = dv_params->cen2;
	float dx1,dy1;
	float dx2,dy2;

	dx1 = coords[0] - cen1[0]; 
	dy1 = coords[1] - cen1[1];
	dx2 = coords[0] - cen2[0]; 
	dy2 = coords[1] - cen2[1];

	d1 = sqrt(sqr(dx1) + sqr(dy1));
	d2 = sqrt(sqr(dx2) + sqr(dy2));

	s1 = dv_params->i1/2.0/PI/d1;
	s2 = dv_params->i2/2.0/PI/d2;

	vel[0] =  s1*dy1/d1 + s2*dy2/d2;
	vel[1] = -s1*dx1/d1 - s2*dx2/d2;
}	/* end test_double_vortex_vel */
