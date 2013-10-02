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
*				gbstate.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routine for the evaluation of boundary states, in particular
*	the evaluation of states to be used in the imposition
*	of boundary conditions.
*
*/


#include <gdecs/gdecs.h>

/* LOCAL Function prototypes */
LOCAL	void	nearest_state_and_normal_on_hypersurface(float*,float*,
				 COMPONENT,HYPER_SURF*,Locstate,float*,Front*);
//#bjet2
LOCAL   void    connect_by_riemann_wave(Locstate, const Locstate, int, const float*);
LOCAL   void    nearest_interior_boundary_state(float*,HYPER_SURF*,Front*,POINTER,Locstate);


/*
*			flow_through_boundary_state():
*	
*	Assigns the boundary state at the exterior point coords to be the state
*	at the nearest interface point to coords..
*
*/

/*ARGSUSED*/
EXPORT	void flow_through_boundary_state(
	float		*coords,
	HYPER_SURF	*hs,
	Front		*front,
	POINTER		p2wave,
	Locstate	state)
{
	INTERFACE	*intfc = hs->interface;
	COMPONENT	ext_comp = exterior_component(intfc);
	HYPER_SURF	*hs_on = NULL;
	HYPER_SURF_ELEMENT *hse_on = NULL;
	bool pc_excluded = is_excluded_comp(positive_component(hs),intfc); 
	bool nc_excluded = is_excluded_comp(negative_component(hs),intfc); 
	float		coords_on[MAXD];
	float		t[MAXD];

	debug_print("flow_through","Entered flow_through_boundary_state()\n");

	if (pc_excluded && !nc_excluded)
	    ext_comp = positive_component(hs);
	else if (!pc_excluded && nc_excluded)
	    ext_comp = negative_component(hs);
	else
	    ext_comp = exterior_component(intfc);

	if (intfc->modified)
	{
	    HYPER_SURF *chs = NULL;
	    if ((intfc != front->interf) && (!front->interf->modified))
	    {
		chs = find_correspond_hyper_surface(hs,NULL,NULL,front,
			                            front->interf);
		if (chs != NULL) 
		{ 
		    if (!nearest_interface_point(coords,ext_comp,
		                    front->interf,INCLUDE_BOUNDARIES,
				    chs,coords_on,t,&hse_on,&hs_on))
		    	chs = NULL;
		    else if (!equivalent_comps(positive_component(hs),
		                          positive_component(chs),intfc) ||
		        !equivalent_comps(negative_component(hs),
		                          negative_component(chs),intfc))
			chs = NULL;
		}
	    }
	    if ((chs == NULL) && !long_nearest_interface_point(coords,ext_comp,
			                     intfc,INCLUDE_BOUNDARIES,hs,
					     coords_on,t,&hse_on,&hs_on))
	    {
	    	screen("ERROR in flow_through_boundary_state(), "
	    	       "long_nearest_interface_point() failed\n");
	    	clean_up(ERROR);
	    }
	}
	else 
	{
	    if (!nearest_interface_point(coords,ext_comp,intfc,
					INCLUDE_BOUNDARIES,hs,coords_on,
					t,&hse_on,&hs_on))
	    {
	    	screen("ERROR in flow_through_boundary_state(), "
	    	       "nearest_interface_point() failed\n");
	    	clean_up(ERROR);
	    }
	}
	if (is_excluded_comp(negative_component(hs),intfc))
	    state_along_hypersurface_element(positive_component(hs),
	    				     t,hse_on,hs_on,state);
	else if (is_excluded_comp(positive_component(hs),intfc))
	    state_along_hypersurface_element(negative_component(hs),
					     t,hse_on,hs_on,state);
	else 
	{
	    screen("ERROR: in flow_through_boundary_state\n"
	           "No boundary component on boundary\n");
	    clean_up(ERROR);
	}

/*	if (is_gravity() == YES)
	{
	    int    i, dim = front->rect_grid->dim;
	    float  v[MAXD], d, g, c, u, tmp, dens;
	    const float  eps = 10.0*MACH_EPS; //TOLERANCE

	    for (i = 0; i < dim; ++i)
		v[i] = coords[i] - coords_on[i];
	    d = mag_vector(v,dim);
	    g = scalar_product(gravity(coords,front->time),v,dim);
	    dens = Dens(state);
	    u = scalar_product(Mom(state),v,dim)/Dens(state);
	    c = d*sound_speed(state);
	    tmp = c*c - u*u;
	    if (fabs(tmp) < eps)
	    {
	    	(void) printf("WARNING in flow_through_boundary_state(), "
			      "velocity equals sound speed.\n");
	    	tmp = tmp > 0 ? eps : -eps;
	    }
	    dens += 3.0*d*d*dens*g/tmp;
	    for (i = 0; i < dim; ++i)
	    	v[i] = vel(i,state);
	    state_on_adiabat_with_dens(state,dens,state,GAS_STATE);
	    add_velocity_to_state(state,v);
	}
*/

	if (debugging("flow_through"))
	{
	    (void) printf("State returned by flow_through_boundary_state()\n");
	    fprint_raw_gas_data(stdout,state,front->rect_grid->dim);
	}

	debug_print("flow_through","Left flow_through_boundary_state()\n");
}		/*end flow_through_boundary_state*/

/*
*			constant_pressure_flow_through_boundary_state():
*	
*	Assigns the boundary state at the exterior point coords to be the state
*	at the nearest interface point to coords..
*
*/

/*ARGSUSED*/
EXPORT	void constant_pressure_flow_through_boundary_state(
	float		*coords,
	HYPER_SURF	*hs,
	Front		*front,
	POINTER		p2wave,
	Locstate	state)
{
	Locstate bst = boundary_state(hs);
	int      st_type;

	debug_print("flow_through",
		"Entered constant_pressure_flow_through_boundary_state()\n");
	flow_through_boundary_state(coords,hs,front,p2wave,state);
	st_type = state_type(state);
	set_state(state,TGAS_STATE,state);
	Dens(state) = Dens(bst);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int  i;
            if(Params(state)->n_comps != 1)
            {
                for(i = 0; i < Params(state)->n_comps; i++)
                    pdens(state)[i] = pdens(bst)[i];
            }
        }
	Press(state) = pressure(bst);
	set_state(state,st_type,state);
	if (debugging("flow_through"))
	{
	    (void) printf("State returned by ");
	    (void) printf("constant_pressure_flow_through_boundary_state()\n");
	    fprint_raw_gas_data(stdout,state,front->rect_grid->dim);
	}
	debug_print("flow_through",
		"Left constant_pressure_flow_through_boundary_state()\n");
}		/*end constant_pressure_flow_through_boundary_state*/

/*
*			time_dep_pressure_flow_through_boundary_state():
*	
*	Assigns the boundary state at the exterior point coords to be the state
*	at the nearest interface point to coords..
*
*/

/*ARGSUSED*/
/*	
 *	 time_dep_pressure_flow_through_boundary_state();
 *
 *      Return state whose pressure at the time:
 *
 *        __________         peak_pressure (pr2)
 *       /          \
 *      /            \
 *      --------------       base_pressure (pr1)
 *     0 t1        t2 t3
 *
 *     tw - 0  = time for warming-up   (tw)
 *     tp - tw = time for peak         (tp)
 *     tc - tp = time for cooling-down (tc)
 *
*/

void	get_constant_state(Locstate, int, float*, float);

EXPORT	void time_dep_pressure_flow_through_boundary_state(
	float		*coords,
	HYPER_SURF	*hs,
	Front		*front,
	POINTER		p2wave,
	Locstate	state)
{
	int      st_type;
	FD_DATA *fd_data = (FD_DATA *)boundary_state_data(hs);
	float	tr = fd_data->tr;
	float	tp = fd_data->tp;
	float	ts = fd_data->ts;
	float	pr_a = fd_data->pr_a;
	float	pr_p = fd_data->pr_p;
	float	time = front->time;
	Locstate ref_st = fd_data->state;
	float    v[MAXD], p;
	static   Locstate tmp = NULL;
	int      i, dim = front->rect_grid->dim;

	debug_print("time_dep",
		"Entered time_dep_pressure_flow_through_boundary_state()\n");
        

	//TMP
        //add_to_debug("time_dep");	
        if(debugging("time_dep_bdry"))
        {
            printf("ENtered time_dep_pressure_flow_through_boundary_state\n");
            printf(" ext_comp = %d,neg_comp, %d pos_comp %d\n",
                     exterior_component(hs->interface),
                     negative_component(hs),positive_component(hs));
        }

	if(!tmp)
	    scalar(&tmp,size_of_state(front->interf));
	//#bjet2 
	zero_scalar(tmp, size_of_state(front->interf));

	// flow_through_boundary_state(coords,hs,front,p2wave,tmp);
        nearest_interior_boundary_state(coords,hs,front,p2wave,tmp);
	st_type = state_type(tmp);
	
	Set_params(state,tmp);

	//TMP_TK	
	//if(debugging("const_bdry"))
	{
	    ft_assign(state, tmp, front->sizest);
	    get_constant_state(state, 2, coords, front->time);

	    debug_print("time_dep",
		"Left time_dep_pressure_flow_through_boundary_state()\n");
	    return;
	}

	set_type_of_state(state,TGAS_STATE);
	
	//TMP liuxt12
	if (debugging("strange_state"))
	    verbose_print_state("time_dep_pressure_flow_through_boundary_state",state);

	if (time < tr)	          p = (pr_p-pr_a)/tr*time + pr_a;
	else if (time <= tr+tp)   p = pr_p;
	else if (time < tr+tp+ts) p = -(pr_p-pr_a)/ts*(time-tr-tp) + pr_p;
	else                      p = pr_a;

        state_on_adiabat_with_pr(ref_st,p,state,TGAS_STATE);
	
        if (dim == 1)
        {
            float v_adj;
            float v = Vel(tmp)[0];
            float c = sound_speed(tmp);
                                                                                                                                             
            if (v >= c) // supersonic
              v_adj = 0.0;
            else
              v_adj = riemann_wave_curve(tmp,p);
            if(debugging("old_flowthrough"))
              v_adj = 0.0;
            Vel(state)[0] = Vel(tmp)[0] - v_adj;
        }
        else if (dim == 2)
        {
            float nor[MAXD];
            CURVE *c;
            c = Curve_of_hs(hs);
            if (is_bdry(c)) //Rectangular boundary
            {
                float t[MAXD], length;
                int i;
                for (i = 0; i < dim; ++i)
                  t[i] = Coords(c->end->posn)[i] - Coords(c->start->posn)[i];
                length = mag_vector(t,dim);
                for (i = 0; i < dim; ++i)
                  t[i] /= length;
                nor[0] = t[1];
                nor[1] = -t[0];
                // should be (0,-1,0) for jet run
            }
            else
            {
                printf("Not impelemented in %s(%d).\n",__FILE__,__LINE__);
                clean_up(ERROR);
            }
            connect_by_riemann_wave(state,tmp,dim,nor);
        }
        else
        {
            // NOTE: wired normal, should compute from hypersurf
            float nor[MAXD] = {0, 0, -1.0};
            // TMP
            float vn;
            if(debugging("time_dep_bdry"))
            {
                vn = scalar_product(VelocityVector(tmp,NULL),nor,dim);
                printf("time_dep_pressure_flow_through_boundary_state(),"
                   " flow_through state, press by time %10.8g\n", p);
                printf("normal vel %g\n",vn);
            }

            //here set tmp to TGAS_STATE, tmp is static, it will keep this type 
	    connect_by_riemann_wave(state,tmp,dim,nor);

            if(debugging("time_dep_bdry"))
            {
                verbose_print_state("flow through st", tmp);
                verbose_print_state("final boundary st", state);
                printf("coords %g %g %g\n", coords[0], coords[1], coords[2]);
            }
        }
        //st_type comes from tmp previously
	//set_state(state,st_type,state);
	
	//TMP as vel inflow
	/*
	set_state(state,TGAS_STATE,tmp);
	
	Dens(state) = 0.6608;
	Press(state) = 1.013;
	Vel(state)[0] = 0.0;
	Vel(state)[1] = 0.0;
	Vel(state)[2] = -p;
	set_state(state,st_type,state);
	*/

	if (debugging("time_dep"))
	{
	    if(front->rect_grid->dim == 3)	
	        printf("coords [%g %g %g] , Zvel %g ", 
			coords[0], coords[1], coords[2], Mom(state)[2]);	
	    (void) printf("time = %f  pressure = %f\n",
	    		time,pressure(state));
	}
	//remove_from_debug("time_dep");
	debug_print("time_dep",
		"Left time_dep_pressure_flow_through_boundary_state()\n");
}	/*end time_dep_pressure_flow_through_boundary_state */


/*
*			g_neumann_bdry_state():
*
*	Imposes a reflecting boundary condition at points laying outside of a 
*	neumann boundary.  The point (x,y) is reflected through the neumann 
*	boundary curve Nbdry, the state at this reflected point (xref, yref) 
*	is evaluated using hyp_solution and the state at (x, y) is then set to
*	be the state with the same density, pressure and tangential (wrt 
*	Nbdry) component of velocity as the state at (xref, yref) and opposite
*	normal component of velocity.  It is assumed that the reflected point 
*	should lie in the component int_comp.
*
*	The use of this function may be problem dependent and it is accessed
*	through the function pointer front->neumann_bdry_state.
*
*       The method mentioned above is one of several methods implemented
*       in this function, and it is the default.  The other methods are
*       accessible by setting an appropriate debug string.
*/

/*ARGSUSED*/
EXPORT	int g_neumann_bdry_state(
	float		*coords,
	COMPONENT	int_comp,
	POINT		*pt,
	HYPER_SURF	*Nbdry,
	Front		*front,
	POINTER		p2wave,
	Locstate	state)
{
	Wave		*wave = (Wave*)p2wave;
	float		nor[MAXD];
	float		coords_ref[MAXD];
	float		coords_tmp[MAXD];

	int		i, dim = front->rect_grid->dim;

	float           tmp_vel[MAXD], vn;

	/* coords_ref is the location of the reflection of the point coords
	   across the Neumann boundary. */

	if (debugging("no_nbdry") 
	    || !reflect_pt_about_Nbdry(coords,coords_ref,nor,
					  int_comp,Nbdry,front))
	    return NO;

	/* Get the state at coords_ref. */
	hyp_solution(coords_ref,int_comp,Nbdry,UNKNOWN_SIDE,
		     front,wave,state,NULL);

	set_state(state,TGAS_STATE,state);
	vn = 0.0;
	for (i = 0; i < dim; i ++)
	{
	    tmp_vel[i] = Vel(state)[i];
	    vn += tmp_vel[i] * nor[i];
	}

	if (no_slip(Nbdry))
        {
	    float alpha = adherence_coeff(Nbdry);
            for (i = 0; i < dim; i ++)
	    	Vel(state)[i] = (1-2*alpha)*Vel(state)[i];
	}

	zero_normal_velocity(state,nor,dim);
	for (i = 0; i < dim; i ++)
	    Vel(state)[i] += - vn * nor[i];
	set_state(state,GAS_STATE,state);
	return YES;

	/* multiple versions */

	if (debugging("nbdry1"))
	{
	    /* This version extrapolates the density and pressure of the wall
	       state at coords_tmp, but reflects the velocity of the state at
	       coords_ref. */

	    float v_dot_n;
	    static Locstate tmpst = NULL;
	    if (tmpst == NULL)
		alloc_state(front->interf,&tmpst,front->sizest);
	

	    /* Reflect the velocity of the state at coords_ref. */
	    
	    set_state(state,TGAS_STATE,state);

	    v_dot_n = 2.0*scalar_product(Vel(state),nor,dim);
	    for (i = 0; i < dim; i ++) 
		Vel(state)[i] -= v_dot_n*nor[i];

	    /* Get the wall state at coords_tmp, which is the midpoint of the
	       line segment connecting coords to coords_ref. */

	    for (i = 0; i < dim; ++i) 
		coords_tmp[i] = 0.5*(coords[i]+coords_ref[i]);

	    hyp_solution(coords_tmp,int_comp,Nbdry,UNKNOWN_SIDE,front,
			 wave,tmpst,NULL);

	    /* Copy the density and pressure of tmpst to state. */

	    Dens(state) = Dens(tmpst);
	    Press(state) = pressure(tmpst);
	    set_state(state,GAS_STATE,state);
	}
	else if (debugging("nbdry2"))
	{
	    /* This version linearly extrapolates the TGAS_STATEs at the
	       reflection point and boundary (coords_ref and coords_tmp,
	       respectively). */

	    /* Get the wall state at coords_tmp, which is the midpoint of the
	       line segment connecting coords to coords_ref. */

	    static Locstate tmpst = NULL;
	    if (tmpst == NULL)
		alloc_state(front->interf,&tmpst,front->sizest);

	    for (i = 0; i < dim; ++i) 
		coords_tmp[i] = 0.5*(coords[i]+coords_ref[i]);

	    hyp_solution(coords_tmp,int_comp,Nbdry,UNKNOWN_SIDE,front,
			 wave,tmpst,NULL);
	    
	    /* Linearly extrapolate the density, pressure, and velocity. */

	    set_state(state,TGAS_STATE,state);
	    set_state(tmpst,TGAS_STATE,tmpst);
	    
	    Dens(state) = 2.0*Dens(tmpst) - Dens(state);
	    Press(state) = 2.0*Press(tmpst) - Press(state);
	    for (i = 0; i < dim; ++i) 
		Vel(state)[i] = 2.0*Vel(tmpst)[i] - Vel(state)[i];
	    
	    set_state(state,GAS_STATE,state);
	}
	else if (debugging("nbdry3"))
	{
	    /* This version reflects the velocity of the state at coords_ref
	       and linearly extrapolates the density and pressure at coords_ref
	       and at one cell further in.  The state at the wall is not
	       used. */

	    float v_dot_n;
	    static Locstate tmpst = NULL;

	    if (tmpst == NULL)
		alloc_state(front->interf,&tmpst,front->sizest);

	    /* Reflect the velocity of the state at coords_ref. */
	    
	    set_state(state,TGAS_STATE,state);

	    v_dot_n = 2.0*scalar_product(Vel(state),nor,dim);
	    for (i = 0; i < dim; i ++) 
		Vel(state)[i] -= v_dot_n*nor[i];

	    /* Get the location and state one cell further in. */

	    for (i = 0; i < dim; ++i) 
		coords_tmp[i] = 2.0*coords_ref[i] - coords[i];
	    
	    hyp_solution(coords_tmp,int_comp,Nbdry,UNKNOWN_SIDE,front,
			 wave,tmpst,NULL);
	    
	    /* Linearly extrapolate the density and pressure. */

	    Dens(state) = 2.0*Dens(state) - Dens(tmpst);
	    Press(state) = 2.0*Press(state) - pressure(tmpst);
	    
	    set_state(state,GAS_STATE,state);
	}
	else if (debugging("nbdry4") || debugging("nbdry4a"))
	{
	    /* This version does a quadratic extrapolation using the states
	       at coords_ref and at one and two cells further in.  The state
	       at the wall is not used. */

	    static Locstate tmpst = NULL, tmpst1 = NULL;
	    if (tmpst1 == NULL)
	    {
		alloc_state(front->interf,&tmpst,front->sizest);
		alloc_state(front->interf,&tmpst1,front->sizest);
	    }

	    /* Get the location and state one cell further in. */

	    for (i = 0; i < dim; ++i) 
		coords_tmp[i] = 2.0*coords_ref[i] - coords[i];
	    
	    hyp_solution(coords_tmp,int_comp,Nbdry,UNKNOWN_SIDE,front,
			 wave,tmpst,NULL);
	    
	    /* Get the location and state two cells further in. */

	    for (i = 0; i < dim; ++i) 
		coords[i] = 2.0*coords_tmp[i] - coords_ref[i];
	    
	    hyp_solution(coords,int_comp,Nbdry,UNKNOWN_SIDE,front,
			 wave,tmpst1,NULL);
	    
	    /* Check to see if the densities are monotonic.
	    
	    if ((Dens(state) > Dens(tmpst) && Dens(tmpst) < Dens(tmpst1))
	     || (Dens(state) < Dens(tmpst) && Dens(tmpst) > Dens(tmpst1)))
	    {
		(void) printf("WARNING in g_neumann_bdry_state() --\n");
		(void) printf("Dens(state) = %g\n", Dens(state));
		(void) printf("Dens(tmpst) = %g\n", Dens(tmpst));
		(void) printf("Dens(tmpst1) = %g\n", Dens(tmpst1));
		(void) printf("Densities are not monotonic near wall.\n");
	    }
	    */

	    set_state(state,TGAS_STATE,state);

	    /* Yet another variant:  Reflect or extrapolate the velocity. */

	    if (debugging("nbdry4a"))  /* reflect the velocity */
	    {
		float v_dot_n = 2.0*scalar_product(Vel(state),nor,dim);
		for (i = 0; i < dim; i ++) 
		    Vel(state)[i] -= v_dot_n*nor[i];
	    }
	    else                       /* extrapolate the velocity */
	    {
		for (i = 0; i < dim; ++i)
		    Vel(state)[i] = vel(i,tmpst1) - 3.0*vel(i,tmpst)
			+ 3.0*Vel(state)[i];
	    }

	    /* Extrapolate the density and pressure. */

	    Dens(state) = Dens(tmpst1) - 3.0*Dens(tmpst) + 3.0*Dens(state);
	    Press(state) = pressure(tmpst1) - 3.0*pressure(tmpst)
		+ 3.0*Press(state);

	    set_state(state,GAS_STATE,state);
	}
	else if (debugging("nbdry5"))
	{
	    /* This version does a quadratic extrapolation using the states
	       at the wall, at coords_ref and at one cell further in. */

	    static Locstate tmpst = NULL, tmpst1 = NULL;
	    if (tmpst == NULL)
	    {
		alloc_state(front->interf,&tmpst,front->sizest);
		alloc_state(front->interf,&tmpst1,front->sizest);
	    }

	    /* Get the location and state one cell further in. */

	    for (i = 0; i < dim; ++i) 
		coords_tmp[i] = 2.0*coords_ref[i] - coords[i];
	    
	    hyp_solution(coords_tmp,int_comp,Nbdry,UNKNOWN_SIDE,front,
			 wave,tmpst1,NULL);
	    
	    /* Get the location and state at the wall. */

	    for (i = 0; i < dim; ++i) 
		coords_tmp[i] = 0.5*(coords[i]+coords_ref[i]);

	    hyp_solution(coords,int_comp,Nbdry,UNKNOWN_SIDE,front,
			 wave,tmpst,NULL);
	    
	    set_state(state,TGAS_STATE,state);

	    /* Extrapolate the density, pressure, and velocity. */

	    Dens(state) = (Dens(tmpst1) - 6.0*Dens(state) 
			   + 8.0*Dens(tmpst))/3.0;
	    Press(state) = (pressure(tmpst1) - 6.0*pressure(state)
			    + 8.0*pressure(tmpst))/3.0;
	    for (i = 0; i < dim; ++i)
		Vel(state)[i] = (vel(i,tmpst1) - 6.0*Vel(state)[i]
		    + 8.0*vel(i,tmpst))/3.0;
	    
	    set_state(state,GAS_STATE,state);
	}
	else
	{

	    /* Default method: Determine the state according to the comment above
	       this function.  If there is gravity, then reflect_state() calls
	       g_reflect_and_stratify_state(), which stratifies the density and
	       pressure with respect to the state being reflected. */
    
	    /* coords_tmp is the point on the reflection plane. */

	    for (i = 0; i < dim; ++i) 
	        coords_tmp[i] = 0.5*(coords[i]+coords_ref[i]);

	    reflect_state(state,front->interf,coords,coords_tmp,nor);

	}
	return YES;
}		/*end g_neumann_bdry_state*/



/*
*		g_fixed_boundary_state():
*	
*	Sets the state to the boundary_state of the given hypersurface
*	to the value of the fixed state boundary_state(hs).
*	If the boundary state params are different from params on
*	the hypersurface (an indication that a contact has broken
*	through on the boundary) then the boundary state is filtered
*	through a Riemann problem so the the params of the returned
*	state is consistent with those in the adjacent interior region.
*/

/*ARGSUSED*/
EXPORT	void g_fixed_boundary_state(
	float		*coords,
	HYPER_SURF	*hs,
	Front		*front,
	POINTER		p2wave,
	Locstate	state)
{
	Locstate		bs = boundary_state(hs);
	float			W[3];
	float			nor[3];
	float			coords_on[3];
	int			dim = front->rect_grid->dim;
	static	Locstate	st = NULL, tmpst = NULL;

	//TMP fixed states
	ft_assign(state,bs,front->sizest);
	return;

	if (is_exterior_comp(negative_component(hs),hs->interface))
	{
		switch (dim)
		{
#if defined(ONED)
	    case 1:
	    	if (Same_params(right_state(Point_of_hs(hs)),bs))
	    	{
	    		ft_assign(state,bs,front->sizest);
	    		return;
	    	}
	    	break;
#endif /* defined(ONED) */
#if defined(TWOD)
	    case 2:
	    	if (Same_params(right_start_state(Curve_of_hs(hs)),bs))
	    	{
	    	    ft_assign(state,bs,front->sizest);
	    	    return;
	    	}
		break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
	    {
	        SURFACE	*s = Surface_of_hs(hs);
	        Locstate	sl, sr;
	        TRI	*tri = first_tri(s);
	        slsr(Point_of_tri(tri)[0],Hyper_surf_element(tri),hs,&sl,&sr);
	        if (Same_params(sr,bs))
	        {
	            ft_assign(state,bs,front->sizest);
	            return;
	        }
	    }
	        break;
#endif /* defined(THREED) */
	    }
	    if (st == NULL)
	    {
	    	alloc_state(front->interf,&st,front->sizest);
	    	alloc_state(front->interf,&tmpst,front->sizest);
	    }
	    nearest_state_and_normal_on_hypersurface(coords,coords_on,
			positive_component(hs),hs,st,nor,front);
	    w_speed(coords_on,bs,st,tmpst,state,W,0.0,nor,CONTACT,front);
	    return;
	}
	else if (is_exterior_comp(positive_component(hs),hs->interface))
	{
	    switch (dim)
	    {
#if defined(ONED)
	    case 1:
	    	if (Same_params(left_state(Point_of_hs(hs)),bs))
	    	{
	    	    ft_assign(state,bs,front->sizest);
	    	    return;
	    	}
	    	break;
#endif /* defined(ONED) */
#if defined(TWOD)
	    case 2:
	    	if (Same_params(left_start_state(Curve_of_hs(hs)),bs))
		{
		    ft_assign(state,bs,front->sizest);
		    return;
		}
		break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
	    {
	        SURFACE	*s = Surface_of_hs(hs);
	    	Locstate	sl, sr;
	    	TRI	*tri = first_tri(s);
	    	slsr(Point_of_tri(tri)[0],Hyper_surf_element(tri),hs,&sl,&sr);
	    	if (Same_params(sl,bs))
	    	{
	    	    ft_assign(state,bs,front->sizest);
	    	    return;
	    	}
	    }
	    	break;
#endif /* defined(THREED) */
	    }
	    if (st == NULL)
	    {
	    	alloc_state(front->interf,&st,front->sizest);
	    	alloc_state(front->interf,&tmpst,front->sizest);
	    }
	    nearest_state_and_normal_on_hypersurface(coords,coords_on,
			negative_component(hs),hs,st,nor,front);
	    w_speed(coords_on,st,bs,state,tmpst,W,0.0,nor,CONTACT,front);
	    return;
	}
}		/*end g_fixed_boundary_state*/

/*
*		g_time_dependent_boundary_state():
*	
*	Evalutes the boundary state at the given time by interpolation
*       from a table of times, densities, pressures, and velocities.
*	If the boundary state params are different from params on
*	the hypersurface (an indication that a contact has broken
*	through on the boundary) then the boundary state is filtered
*	through a Riemann problem so that the params of the returned
*	state is consistent with those in the adjacent interior region.
*/

EXPORT	void g_time_dependent_boundary_state(
	float		*coords,
	HYPER_SURF	*hs,
	Front		*front,
	POINTER		p2wave,
	Locstate	state)
{
	/* variables copied from g_fixed_boundary_state() above */
	float			W[3];
	float			nor[3];
	float		        coords_on[MAXD];
	int			dim = front->rect_grid->dim;
	static	Locstate	st = NULL, tmpst = NULL;
        static  Locstate        bs = NULL;
	/* variables specific to g_time_dependent_boundary_state() */
        TD_BSTATE               *tbs = (TD_BSTATE*)boundary_state_data(hs);
        float                   time = front->time;
        int                     i, n;
        COMPONENT               comp;
        HYPER_SURF_ELEMENT      *hse;
        static POINT            *tmpp;
	float		        t[MAXD];
        if (bs == NULL)
        {
	    alloc_state(front->interf,&bs,front->sizest);
            tmpp = Static_point(front->interf);
        }
        for (i = 0; i < dim; ++i)
          Coords(tmpp)[i] = coords[i];

	if (is_excluded_comp(negative_component(hs),hs->interface))
            comp = positive_component(hs);
	else if (is_excluded_comp(positive_component(hs),hs->interface))
            comp = negative_component(hs);
        else
        {
           screen("ERROR in g_time_dependent_boundary_state() - no exterior "
                  "component\n");
           print_interface(hs->interface);
           clean_up(ERROR);
        }

	if (nearest_interface_point(coords,comp,front->interf,
                                    INCLUDE_BOUNDARIES,hs,
				    coords_on,t,&hse,&hs) != YES)
	{
	    screen("ERROR in nearest_interior_comp(), "
		   "nearest_interface_point() failed\n");
	    clean_up(ERROR);
	}
        normal(tmpp,hse,hs,nor,front);
	if (is_excluded_comp(positive_component(hs),hs->interface))
        {
            for (i = 0; i < dim; ++i)
               nor[i] = -nor[i];
        }

        if (time <= tbs->t[0])
        {
            Dens(bs) = tbs->rho[0];
            Press(bs) = tbs->p[0];
            for (i = 0; i < dim; ++i)
               Vel(bs)[i] = tbs->v[0]*nor[i];
            Params(bs) = tbs->params;
            set_type_of_state(bs,TGAS_STATE);
            set_state(bs,GAS_STATE,bs);
        }
        else if (tbs->t[tbs->n-1] <= time)
        {
            Dens(bs) = tbs->rho[tbs->n-1];
            Press(bs) = tbs->p[tbs->n-1];
            for (i = 0; i < dim; ++i)
               Vel(bs)[i] = tbs->v[tbs->n-1]*nor[i];
            Params(bs) = tbs->params;
            set_type_of_state(bs,TGAS_STATE);
            set_state(bs,GAS_STATE,bs);
        }
        else
        {
            for (n = 0; n < tbs->n; ++n)
            {
               if ((n < (tbs->n - 1)) && (tbs->t[n] <= time) && (time < tbs->t[n+1]))
               {
                  float alpha = (time - tbs->t[n])/(tbs->t[n+1]-tbs->t[n]);
                  float v;
                  Dens(bs) = (1.0-alpha)*tbs->rho[n]+alpha*tbs->rho[n+1];
                  Press(bs) = (1.0-alpha)*tbs->p[n]+alpha*tbs->p[n+1];
                  v = (1.0-alpha)*tbs->v[n] + alpha*tbs->v[n+1];
                  for (i = 0; i < dim; ++i)
                     Vel(bs)[i] = v*nor[i];
                  Params(bs) = tbs->params;
                  set_type_of_state(bs,TGAS_STATE);
                  set_state(bs,GAS_STATE,bs);
		  break;
               }
            }
         }
       

	if (is_excluded_comp(negative_component(hs),hs->interface))
	{
		switch (dim)
		{
#if defined(ONED)
	    case 1:
	    	if (Same_params(right_state(Point_of_hs(hs)),bs))
	    	{
	    		ft_assign(state,bs,front->sizest);
	    		return;
	    	}
	    	break;
#endif /* defined(ONED) */
#if defined(TWOD)
	    case 2:
	    	if (Same_params(right_start_state(Curve_of_hs(hs)),bs))
	    	{
	    	    ft_assign(state,bs,front->sizest);
	    	    return;
	    	}
		break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
	    {
	        SURFACE	*s = Surface_of_hs(hs);
	        Locstate	sl, sr;
	        TRI	*tri = first_tri(s);
	        slsr(Point_of_tri(tri)[0],Hyper_surf_element(tri),hs,&sl,&sr);
	        if (Same_params(sr,bs))
	        {
	            ft_assign(state,bs,front->sizest);
	            return;
	        }
	    }
	        break;
#endif /* defined(THREED) */
	    }
	    if (st == NULL)
	    {
	    	alloc_state(front->interf,&st,front->sizest);
	    	alloc_state(front->interf,&tmpst,front->sizest);
	    }
	    nearest_state_and_normal_on_hypersurface(coords,coords_on,
			positive_component(hs),hs,st,nor,front);
	    w_speed(coords_on,bs,st,tmpst,state,W,0.0,nor,CONTACT,front);
	    return;
	}
	else if (is_excluded_comp(positive_component(hs),hs->interface))
	{
	    switch (dim)
	    {
#if defined(ONED)
	    case 1:
	    	if (Same_params(left_state(Point_of_hs(hs)),bs))
	    	{
	    	    ft_assign(state,bs,front->sizest);
	    	    return;
	    	}
	    	break;
#endif /* defined(ONED) */
#if defined(TWOD)
	    case 2:
	    	if (Same_params(left_start_state(Curve_of_hs(hs)),bs))
		{
		    ft_assign(state,bs,front->sizest);
		    return;
		}
		break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
	    {
	        SURFACE	*s = Surface_of_hs(hs);
	    	Locstate	sl, sr;
	    	TRI	*tri = first_tri(s);
	    	slsr(Point_of_tri(tri)[0],Hyper_surf_element(tri),hs,&sl,&sr);
	    	if (Same_params(sl,bs))
	    	{
	    	    ft_assign(state,bs,front->sizest);
	    	    return;
	    	}
	    }
	    	break;
#endif /* defined(THREED) */
	    }
	    if (st == NULL)
	    {
	    	alloc_state(front->interf,&st,front->sizest);
	    	alloc_state(front->interf,&tmpst,front->sizest);
	    }
	    nearest_state_and_normal_on_hypersurface(coords,coords_on,
			negative_component(hs),hs,st,nor,front);
	    w_speed(coords_on,st,bs,state,tmpst,W,0.0,nor,CONTACT,front);
	    return;
	}
	ft_assign(state,bs,front->sizest);
}		/*end g_time_dependent_boundary_state*/

LOCAL	void	nearest_state_and_normal_on_hypersurface(
	float		*coords,
	float		*coords_on,
	COMPONENT	comp,
	HYPER_SURF	*hs,
	Locstate	state,
	float		*nor,
	Front		*front)
{
	HYPER_SURF_ELEMENT	*hse;
	HYPER_SURF		*hs_on;
	float			t[3];
	int			dim = front->rect_grid->dim;
	static	float		**nor_p = NULL;
#if defined(TWOD) || defined(THREED)
	float			mag;
#endif /* defined(TWOD) || defined(THREED) */

	if (nor_p == NULL)
		bi_array(&nor_p,3,3,FLOAT);
	if (nearest_interface_point(coords,comp,hs->interface,
				    INCLUDE_BOUNDARIES,hs,coords_on,t,&hse,
				    &hs_on) != YES)
	{
	   screen("ERROR in nearest_state_and_normal_on_hypersurface(), "
	          "nearest_interface_point() failed\n");
	   clean_up(ERROR);
	}
	switch (dim)
	{
#if defined(ONED)
	case 1:
	    normal(Point_of_hs(hs),hse,hs,nor,front);
	    break;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	    normal(Bond_of_hse(hse)->start,hse,hs,nor_p[0],front);
	    normal(Bond_of_hse(hse)->end,hse,hs,nor_p[1],front);
	    t[1] = 1.0 - t[0];
	    nor[0] = t[0]*nor_p[0][0] + t[1]*nor_p[1][0];
	    nor[1] = t[0]*nor_p[0][1] + t[1]*nor_p[1][1];
	    mag = mag_vector(nor,dim);
	    nor[0] /= mag;
	    nor[1] /= mag;
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	    normal(Point_of_tri(Tri_of_hse(hse))[0],hse,hs,nor_p[0],front);
	    normal(Point_of_tri(Tri_of_hse(hse))[1],hse,hs,nor_p[1],front);
	    normal(Point_of_tri(Tri_of_hse(hse))[2],hse,hs,nor_p[2],front);
	    nor[0] = t[0]*nor_p[0][0] + t[1]*nor_p[1][0] + t[2]*nor_p[2][0];
	    nor[1] = t[0]*nor_p[0][1] + t[1]*nor_p[1][1] + t[2]*nor_p[2][1];
	    nor[2] = t[0]*nor_p[0][2] + t[1]*nor_p[1][2] + t[2]*nor_p[2][2];
	    mag = mag_vector(nor,dim);
	    nor[0] /= mag;
	    nor[1] /= mag;
	    nor[2] /= mag;
	    break;
#endif /* defined(THREED) */
	}
	state_along_hypersurface_element(comp,t,hse,hs,state);
}		/*end nearest_state_and_normal_on_hypersurface*/

//#bjet2
LOCAL void connect_by_riemann_wave(
        Locstate        lstate,
        const Locstate  rstate,
        int             dim,
        const float     *nor)
{
        float v, vn, c;
        float v_adj=0.0;
        int i;

        set_state(lstate,TGAS_STATE,lstate);
        set_state(rstate,TGAS_STATE,rstate);

        vn = scalar_product(Vel(rstate),nor,dim);
        c = sound_speed(rstate);
        if( vn>=c ) // supersonic
        {
            v_adj = 0.0;
        }
        else
        {
            v_adj = riemann_wave_curve(rstate,Press(lstate));
        }
        v = vn + v_adj;
        for(i=0;i<dim;++i)
        {
            // Tangential velocity of lstate and rstate the same?
            // NEW code, preserve tangential velocity
            // Vel(lstate)[i] = (v-vn)*nor[i] + Vel(rstate)[i];
            /*
               equivalently,
               take the tangential velocity of rstate and
               the computed normal velocity
               float v_tan[MAXD];
               v_tan[i] = Vel(rstate)[i] - vn*nor[i];
               Vel(lstate)[i] = v*nor[i] + v_tan[i];
            */
            // Orignal
            Vel(lstate)[i] = v*nor[i]; // setting tangential velocity zero
        }
        Dens(lstate) = density(lstate);
}

/*
*                       nearest_interior_boundary_state():
*
*       Assigns the boundary state at the exterior point coords to be the state
*       at the nearest interior cell to coords..
*
*/


/*ARGSUSED*/
LOCAL   void nearest_interior_boundary_state(
        float           *coords,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         p2wave,
        Locstate        state)
{
        INTERFACE       *intfc = hs->interface;
        COMPONENT       ext_comp = exterior_component(intfc);
        HYPER_SURF      *hs_on = NULL;
        HYPER_SURF_ELEMENT *hse_on = NULL;
        bool pc_excluded = is_excluded_comp(positive_component(hs),intfc);
        bool nc_excluded = is_excluded_comp(negative_component(hs),intfc);
        float           coords_on[MAXD];
        float           t[MAXD], ncoords[MAXD];
        const  RECT_GRID *rgr = ((Wave*)p2wave)->rect_grid;
        int    i, icoords[MAXD];

        debug_print("flow_through","Entered nearest_interior_boundary_state()\n");

        rect_in_which(coords,icoords,rgr);
        for(i = 0; i < rgr->dim; ++i)
        {
            if (icoords[i] >= rgr->gmax[i]+rgr->ubuf[i])
                icoords[i] = rgr->gmax[i]+rgr->ubuf[i]-1;
            if (icoords[i] < -rgr->lbuf[i])
                icoords[i] = -rgr->lbuf[i];
        }
        for(i = 0; i < rgr->dim; i++)
            ncoords[i] = Rect_coords(icoords,(Wave*)p2wave)[i];
        
	//print_general_vector("crds= ", coords, 3, "\n");
	//print_general_vector("ncrds= ", ncoords, 3, "\n");
	//printf("%d %d %d\n", pc_excluded, 
	//		negative_component(hs), positive_component(hs));
	
	//#bjet2
	if(pc_excluded)
            hyp_solution(ncoords,negative_component(hs),hs,
                 UNKNOWN_SIDE,front,(Wave*)p2wave,state,NULL);
        else
            hyp_solution(ncoords,positive_component(hs),hs,
                 UNKNOWN_SIDE,front,(Wave*)p2wave,state,NULL);

	//print_general_vector("Bcrds= ", coords, 3, "\n");
	//print_gas_state(state);

        if(!is_obstacle_state(state))
            return;

	printf("WARNING in nearest_interior_boundary_state, state is obstacle.\n");

        if (pc_excluded && !nc_excluded)
            ext_comp = positive_component(hs);
        else if (!pc_excluded && nc_excluded)
            ext_comp = negative_component(hs);
        else
            ext_comp = exterior_component(intfc);

        if (intfc->modified)
        {
            HYPER_SURF *chs = NULL;
            if ((intfc != front->interf) && (!front->interf->modified))
            {
                chs = find_correspond_hyper_surface(hs,NULL,NULL,front,
                                                    front->interf);
                if (chs != NULL)
                {
                    if (!nearest_interface_point(coords,ext_comp,
                                    front->interf,INCLUDE_BOUNDARIES,
                                    chs,coords_on,t,&hse_on,&hs_on))
                        chs = NULL;
                    else if (!equivalent_comps(positive_component(hs),
                                          positive_component(chs),intfc) ||
                        !equivalent_comps(negative_component(hs),
                                          negative_component(chs),intfc))
                        chs = NULL;
                }
            }
            if ((chs == NULL) && !long_nearest_interface_point(coords,ext_comp,
                                             intfc,INCLUDE_BOUNDARIES,hs,
                                             coords_on,t,&hse_on,&hs_on))
            {
                screen("ERROR in flow_through_boundary_state(), "
                       "long_nearest_interface_point() failed\n");
                clean_up(ERROR);
            }
        }
        else
        {
            if (!nearest_interface_point(coords,ext_comp,intfc,
                                        INCLUDE_BOUNDARIES,hs,coords_on,
                                        t,&hse_on,&hs_on))
            {
                screen("ERROR in flow_through_boundary_state(), "
                       "nearest_interface_point() failed\n");
                clean_up(ERROR);
            }
        }
        if (is_excluded_comp(negative_component(hs),intfc))
        {
            //state_in_tri
            state_along_hypersurface_element(positive_component(hs),
                                             t,hse_on,hs_on,state);

        }
        else if (is_excluded_comp(positive_component(hs),intfc))
        {
            printf("TEST ERROR,pos_comp %d\n",positive_component(hs));
            state_along_hypersurface_element(negative_component(hs),
                                             t,hse_on,hs_on,state);
        }
        else
        {
            screen("ERROR: in nearest_interior_boundary_state\n"
                   "No boundary component on boundary\n");
            clean_up(ERROR);
        }
        if (is_gravity() == YES)
        {
            int    i, dim = front->rect_grid->dim;
            float  v[MAXD], d, g, c, u, tmp, dens;
            const float  eps = 10.0*MACH_EPS; /*TOLERANCE*/

            for (i = 0; i < dim; ++i)
                v[i] = coords[i] - coords_on[i];
            d = mag_vector(v,dim);
            g = scalar_product(gravity(coords,front->time),v,dim);
            dens = Dens(state);
            u = scalar_product(Mom(state),v,dim)/Dens(state);
            c = d*sound_speed(state);
            tmp = c*c - u*u;
            if (fabs(tmp) < eps)
            {
                (void) printf("WARNING in nearest_interior_boundary_state(), "
                              "velocity equals sound speed.\n");
                tmp = tmp > 0 ? eps : -eps;
            }
            dens += 3.0*d*d*dens*g/tmp;
            for (i = 0; i < dim; ++i)
                v[i] = vel(i,state);
            state_on_adiabat_with_dens(state,dens,state,GAS_STATE);
            add_velocity_to_state(state,v);
        }

        if (debugging("flow_through"))
        {
            (void) printf("State returned by nearest_interior_boundary_state()\n");
            fprint_raw_gas_data(stdout,state,front->rect_grid->dim);
        }

        debug_print("flow_through","Left nearest_interior_boundary_state()\n");
}               /*end nearest_interior_boundary_state*/


/* This file outlines the new flow boundary condition. The icndition uses a boundary state, however, the boundary state is not necessarily used in the stencil. When the flow is supersonic, the boundary state is used unaltered. When the flow is subsonic, a partial Riemann problem is solved fixing some of the fluid variables. There will be functions to allow any 2 variables to be fixed, and solving the third such that no wave enters the domain.
 */

// using gamma law gas

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define RK_PRINT(...) {printf("%s:%d: ",__FILE__,__LINE__);printf(__VA_ARGS__);}

#define dot(a,b) (((a)[0]*(b)[0])+((a)[1]*(b)[1])+((a)[2]*(b)[2]))

#if 1
//TODO: COMMENTS!!!
#include <geos/poly.h>

enum _CRX_SIDE {
        PREV = -1,
        NEXT =  1
};
typedef enum _CRX_SIDE CRX_SIDE;

void interp_flow_bdry_states(
        int      start,
        HYPER_SURF *hs,
        int      side,
        int      endpt,
        int      new_comp,
        Stencil  *sten,
        CRXING   *cross)
{
        int             i, j, indx, in_one;
        size_t          sizest = sten->fr->sizest;
        Locstate        s1, s2;
        int             st_type;
        static Locstate tmpst = NULL;
        if(tmpst == NULL)
            alloc_state(sten->fr->interf,&tmpst,sten->fr->sizest);

        indx = start*side;

        float dir[3] = {0.,0.,0.}; // TODO: impl this
        // TMP DEBUG !!!! RK
	// if x coords is less than 0, then must be inflow bdry
	// else must be outflow boundary
        if (Coords(sten->p[indx])[0] < 0.f)
        {
            dir[0] = -1.;
            s1 = sten->st[indx+1];
        }
        else
        {
            dir[0] = 1.;
            s1 = sten->st[indx-1];
        }

        for (j = 0; j < 3; ++j)
                Coords(sten->p[indx])[j] = Coords(cross->pt)[j];

        s2 = sten->st[indx] = sten->worksp[indx];

        if(Coords(sten->p[indx])[2] > 3.9)
        {
	 	get_constant_state(sten->st[indx], 2, Coords(sten->p[indx]), sten->fr->time);
	//RK_PRINT("Crds indx-i= %d %f %f %f\n", (indx-i),Coords(sten->p[indx-i])[0],Coords(sten->p[indx-i])[1],Coords(sten->p[indx-i])[2]);
	//verbose_print_state("TK indx-i",sten->st[indx-i]);

//	RK_PRINT("Crds indx-1= %f %f %f\n", Coords(sten->p[indx-1])[0],Coords(sten->p[indx-1])[1],Coords(sten->p[indx-1])[2]);
//	RK_PRINT("Crds indx= %f %f %f\n", Coords(sten->p[indx])[0],Coords(sten->p[indx])[1],Coords(sten->p[indx])[2]);
//	    verbose_print_state("TK indx-1",sten->st[indx-1]);
//	    verbose_print_state("TK indx",sten->st[indx]);
	  
/*
            s1 = sten->st[indx-1];
            dir[0] = 0.;
            dir[1] = 0.;
            dir[2] = 1.;
            st_type = state_type(s2);
            set_state(s2, TGAS_STATE, s2);
            set_state(s1, TGAS_STATE, s1);

            set_state(tmpst, TGAS_STATE, s2);

            rk_new_open_bdry_a(Coords(sten->p[indx]),
                        s1,//simst 
                        s2,//bdryst 
                        dir,//dir 
                        tmpst//stenst 
                        );

            set_state(s2, st_type, tmpst);
            set_state(s1, st_type, s1);
*/
            ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
            ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);

            for (i = start+1; i <= endpt; ++i)
            {
                indx = i*side;
                in_one = indx - side;
                sten->st[indx] = sten->st[in_one];
                for (j = 0; j < 3; ++j)
                    Coords(sten->p[indx])[j] = Coords(sten->p[in_one])[j];
                ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
                ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
            }
            return;
        }

        //the inflow dirichlet bdry velocity
        if(new_comp == 2)
        {
            get_constant_state(s2, 2, Coords(sten->p[indx]), sten->fr->time);
            st_type = state_type(s2);
            set_state(s2, TGAS_STATE, s2);
            //no this for the run
            //Dens(s2) = density(s1);
            //Press(s2) = pressure(s1);
            set_state(s2, st_type, s2);
        }

        //crossflow bdry
        if(new_comp == 3)
        {
            st_type = state_type(s2);
            set_state(s2, TGAS_STATE, s2);

            set_state(s1, TGAS_STATE, s1);
            set_state(s2, TGAS_STATE, s2);
            //assert(Dens(s2) == density(s1));
            //assert(Press(s2) == pressure(s1));
            //Dens(s2) = density(s1);
            //Press(s2) = pressure(s1);
            //set_state(s2, st_type, s2);

            set_state(tmpst, TGAS_STATE, s2); 

            // RK NEW !!!!!! VVVVVV
            rk_new_open_bdry(Coords(sten->p[indx]),
                        /*simst*/ s1,
                        /*bdryst*/ s2,
                        /*dir*/ dir,
                        /*stenst*/ tmpst
                        );

            set_state(s2, st_type, tmpst);
            set_state(s1, st_type, s1);
        }

        if(NO && new_comp == 2)
        {
            st_type = state_type(s2);
            set_state(s2, TGAS_STATE, s2);
            Dens(s2) = density(s1);
            Press(s2) = 40.0;
            set_state(s2, st_type, s2);
        }
        ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
        ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);

        for (i = start+1; i <= endpt; ++i)
        {
            indx = i*side;
            in_one = indx - side;
            sten->st[indx] = sten->st[in_one];
            for (j = 0; j < 3; ++j)
                Coords(sten->p[indx])[j] = Coords(sten->p[in_one])[j];
            ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
            ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
        }

        return;
}

void    set_flow_bdry_sten_states(
        int        i,
        HYPER_SURF *hs,
        Stencil    *sten,
        CRX_SIDE   side,
        int endpt)
{
        Wave            *wave = sten->wave;
        Front           *fr = sten->fr;
        int             indx; // == ??
        size_t          sizest = fr->sizest;

        Locstate simst = sten->st[indx-1];
        float dir[3] = {0.,0.,0.}; // TODO: impl this
        // TMP DEBUG !!!! RK
        if (Coords(sten->p[i*side])[0] < 0)
            dir[0] = -1.;
        else
            dir[0] = 1.;


        for (; i <= endpt; ++i)
        {
            indx = i * side;
            sten->st[indx] = sten->worksp[indx];
            // FIXME: this says dirichlet, but is used for flow
            Locstate bdryst; // the theoretical bdry state
            copy_state(bdryst, sten->st[indx]);
            Locstate stenst;  // the state to be used in the stencil

            evaluate_dirichlet_boundary_state(Coords(sten->p[indx]),hs,fr,wave,
                                             bdryst);
                                             //sten->st[indx]);

            rk_new_open_bdry(Coords(sten->p[indx]), simst, bdryst, dir, stenst);            ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
            ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
            simst = stenst;
        }
}
static void flow_bdry_rp_P( double *norm, // out of domain
                                                Locstate st1,
                                                Locstate st4,
                                                Locstate st3)
{
        double T1  = Temperature(st1);
        double P1  = pressure(st1);
        double* U1 = Vel(st1);
        double V1  = -dot(U1,norm);
        double d1  = Dens(st1);

        double T4  = Temperature(st4);
        double P4  = pressure(st4);
        double* U4 = Vel(st4);
        double V4  = -dot(U4,norm);
        double d4  = Dens(st4);

        double T3=T4;
        double U3[3]; memcpy(U3, U4, 3*sizeof(double));
        double V3 = V4;
        double d3 = Dens(st3);

        double P3;

        double gam;
        gam = MGamma(st4);

        // detect shock or raref

        if(0 /*shock*/)
        {
            P3 = P4 + pow(U3 - U4,3)/(1./d4-1./d3);
        }
        else // rarefaction
        {
            // Right rarefaction
            P3= pow(sqr((gam-1.)/2. * (V3-V4) + sound_speed(st4))*d4/(gam*pow(P4,1.0/gam)), 1./(1.-1./gam));

        }

        // copy to st3
        //pressure(st3) = P3;
        //pressure(st3) = P3;
        Press(st3) = P3;


        // set st3 tangential velocity
        return;
}

void rk_new_open_bdry(
        double *crds,
        Locstate sim_state,
        Locstate bdry_state,
        float *dir, // out of domain direction
        /*out*/ Locstate sten_state
)
{
        //RK_PRINT("enter rk_new_open_bdry\n");
        //RK_PRINT("sim_state\n");
        //verbose_print_state("sim_state",sim_state);
        //RK_PRINT("-------------------------\n)");
        //RK_PRINT("bdry_state\n");
        //verbose_print_state("bdry_state", bdry_state);

        double *sim_vel = Vel(sim_state);
        double sim_p = pressure(sim_state);
        double sim_T = Temperature(sim_state);

        double *bdry_vel = Vel(bdry_state); // FIXME???
        double bdry_p = pressure(bdry_state);
        double bdry_T = Temperature(bdry_state);

	double abs_vel[3];
	for(int i=0; i < 3; ++i)
		abs_vel[i]=fabs(sim_vel[i]);
        double m = mach_number(sim_state, sim_vel);
        //RK_PRINT("m=%f\n",m);
        if(m>1)
        {
                // use bdry state.
                //RK_PRINT("using bdry_state\n");
                sten_state = bdry_state;
                return;
        }
/*
	// set in average_p_if_mach
	if(dir[0]==1)
		Press(bdry_state)=wall_pressure_avg; // set to the precomputed global
*/
        double v = dot(sim_vel, dir);
        //RK_PRINT("v=%f\n",v);
        if(v>0) // outflow
        {
//              sten_state = NULL;
                //RK_PRINT("outflow!!\n");

                flow_bdry_rp_P(dir, sim_state, bdry_state, sten_state);
                // ...
        }
        else if (v<=0) // inflow
	{
                //RK_PRINT("inflow!!\n");
//              sten_state = NULL;
                flow_bdry_rp_P(dir, bdry_state, sim_state, sten_state);

        }
        /*
        else // v==0; anomaly
        {
                assert(0);
        }
        */

        //sten_state = bdry_state; // TMP DEBUG!!!!!!

        //RK_PRINT("Done!!\n");
        return;


}
#endif //0
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void rk_new_open_bdry_a(
        double *crds,
        Locstate sim_state,
        Locstate bdry_state,
        float *dir, // out of domain direction
        /*out*/ Locstate sten_state
)
{
        //RK_PRINT("enter rk_new_open_bdry\n");
        //RK_PRINT("sim_state\n");
        //verbose_print_state("sim_state",sim_state);
        //RK_PRINT("-------------------------\n)");
        //RK_PRINT("bdry_state\n");
        //verbose_print_state("bdry_state", bdry_state);

        double *sim_vel = Vel(sim_state);
        double sim_p = pressure(sim_state);
        double sim_T = Temperature(sim_state);

        double *bdry_vel = Vel(bdry_state); // FIXME???
        double bdry_p = pressure(bdry_state);
        double bdry_T = Temperature(bdry_state);

        double m = mach_number(sim_state, sim_vel);
        if(m>1)
        {
                // use bdry state.
                //RK_PRINT("using bdry_state\n");
                sten_state = bdry_state;
                return;
        }

        double v = dot(sim_vel, dir);
        RK_PRINT("v=%f\n",v);
        if(v>0) // outflow
        {
//              sten_state = NULL;
                //RK_PRINT("outflow!!\n");

                flow_bdry_rp_P(dir, bdry_state, sim_state, sten_state);
                // ...
        }
        else if (v<=0) // inflow
        {
                //RK_PRINT("inflow!!\n");
//              sten_state = NULL;
                flow_bdry_rp_P(dir, sim_state, bdry_state, sten_state);

        }
        /*
        else // v==0; anomaly
        {
                assert(0);
        }
        */

        //sten_state = bdry_state; // TMP DEBUG!!!!!!

        //RK_PRINT("Done!!\n");
        return;


}

