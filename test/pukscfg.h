/**
 * Do everything needed to load in the configuration variables
 */

/* Read command line options */
while( ( optc = getopt( argc, argv, "C:p:q:r:e:t:o:x:" ) ) != -1 )
{
        switch( optc )
        {
		case 'C':
			/* Initialize and load the symbol table */
			symtab_init( &sym, SYMTAB_INIT_SIZE );
			ret = cfgread_load_symbols_f( optarg, &sym );
			if( ret == 0 )
			{
				b_config_loaded = 1;
				fprintf( stderr, " * Loaded following symbols from \"%s\":\n\n", optarg );
				symtab_print( &sym );
				fprintf( stderr, "\n" );
			}
			else
				fprintf( stderr, " * ERROR: Failed loading variables from file \"%s\": Returned %d\n\n", optarg, ret );
			break;
		case 'p':
			/* Load points, number, and dimension */
			load_points( optarg, &dim, &npts, &pts, 0 );
			b_points_loaded = 1;

			/* Output notification */
			fprintf( stderr, " * Loaded %d points in %d dimensions from \"%s\"\n\n", npts, dim, optarg );
			break;
		case 'q':
			/* Load quadrature points from file */
			load_quadrature( optarg, &quadn, &qpts, &qwts, &ret );

			/* Make sure nothing went wrong and output what happens */
			if( ret == 0 )
			{
				b_quad_loaded = 1; /* Only mark loaded if there is no trouble */
				fprintf( stderr, " * Loaded quadrature with %d points from \"%s\"\n\n", quadn, optarg );
			}
			else
				fprintf( stderr, " * ERROR: Problem loading quadrature rule from \"%s\": Return value %d\n\n", optarg, ret );
			break;
		case 'r':
			/* Load radial neighbors number */
			rdn = atoi( optarg );
			b_rdn_set = 1;
			fprintf( stderr, " * Setting radial neighbor number to %d\n\n", rdn );
			break;
		case 'e':
			/* Load the number of eigenvalues to use */
			neig = atoi( optarg );
			if( neig > 0 )
			{
				b_neig_set = 1;
				fprintf( stderr, " * Setting the number of eigenvalues to %d\n\n", neig );
			}
			else
				fprintf( stderr, " * ERROR: Invalid number of eigenvalues: %d\n\n", neig );
			break;
		case 't':
			/* Load the name of the routine to use to solve the tridiagonal eigenvalue problem */
			strncpy( s_treig_routine, optarg, TOKEN_BUFFER_LENGTH );
			b_treig_routine_set = 1;
			fprintf( stderr, " * Using \"%s\" to solve the tridiagonal eigenvalue problem\n\n", optarg );
			break;
		case 'o':
			strcpy( buf, optarg );
			nfr = 0;
			n = parse_stokenize( buf, tok, "," );
			for(i=0;i<n;i++)
			{
				strcpy( bug, tok[i] );
				m = parse_stokenize( bug, tol, "-" );
				if( m > 2 ) /* Only want definitive ranges like 127-439 */
					continue;
				if( m == 1 ) /* Then print single page */
				{
					rng[nfr][0] = atol( tol[0] );
					rng[nfr][1] = rng[nfr][0]; /* Same indicating a single frame */
					++nfr;
				}
				else
				{
					rng[nfr][0] = atol( tol[0] );
					rng[nfr][1] = atol( tol[1] );
					if( rng[nfr][1] < rng[nfr][0] )
						fprintf( stderr, " * ERROR: Bad range given; ignoring\n\n" );
					else
						++nfr;
				}
			}
			b_out_eig_range_set = 1; /* Output eigenvector range is specified on command line */
			fprintf( stderr, " * Eigenvector output range specified: " );
			for(i=0;i<nfr;i++)
			{
				if( rng[i][1] == rng[i][0] )
					fprintf( stderr, "%d ", rng[i][0] );
				else
					fprintf( stderr, "%d-%d ", rng[i][0], rng[i][1] );
			}
			fprintf( stderr, "\n\n" );
			break;
		case 'x': /* Load the external potential */
			ret = load_nuclei( optarg, &nuc, 0 );
			if( ret != 0 )
				fprintf( stderr, " ---> ERROR: Failed to load external potential from file \"%s\"\n", optarg );
			else
			{
				fprintf( stderr, " ---> Loaded external potential from file \"%s\"\n\n", optarg );
				nuclei_print( &nuc );
				b_nuc_loaded = 1;
			}
			break;
		case '?':
		default:
			fprintf( stderr, "I don't understand the jibberish you are spouting.\n" );
			return 0;
	}
}

/* Process the configuration variables */
if( b_config_loaded )
{
	/* Print that searching for variables in configuration file */
	fprintf( stderr, " * Reading and executing configuration variables:\n\n" );

	/* Look for points file first if not loaded already */
	if( !b_points_loaded )
	{
		s = symtab_lookup( &sym, "pointCloud" );
		if( s != NULL )
		{
			if( load_points( (char*) s->data, &dim, &npts, &pts, 0 ) != 0 )
				fprintf( stderr, " ---> ERROR: Failed to load configuration file \"%s\"\n", (char*) s->data );
			else
			{
				fprintf( stderr, " ---> Loaded %d points in %d dimensions from \"%s\"\n", npts, dim, (char*) s->data );
				b_points_loaded = 1;
			}
		}
	}
	if( !b_quad_loaded )
	{
		s = symtab_lookup( &sym, "quadrature" );
		if( s != NULL )
		{
			load_quadrature( (char*) s->data, &quadn, &qpts, &qwts, &ret );
			if( ret == 0 )
			{
				fprintf( stderr, " ---> Loaded quadrature rule with %d points from \"%s\"\n", quadn, (char*) s->data );
				b_quad_loaded = 1;
			}
			else
				fprintf( stderr, " ---> ERROR: Problem loading quadurature file in \"%s\": Return value %d\n", (char*) s->data, ret );
		}
	}
	if( !b_rdn_set )
	{
		s = symtab_lookup( &sym, "radialNeighbors" );
		if( s != NULL )
		{
			rdn = atoi( (char*) s->data );
			b_rdn_set = 1;
			fprintf( stderr, " ---> Loaded radial neighbors number %d\n", rdn );
		}
	}
	if( !b_domain_set )
	{
		s = symtab_lookup( &sym, "domain" );
		if( s != NULL )
		{
			load_domain( (char*) s->data, &dm, &ret );
			if( ret == 0 )
			{
				fprintf( stderr, " ---> Loaded domain information\n" );
				b_domain_set = 1;
			}
		}
	}
	s = symtab_lookup( &sym, "boundaryTolerance" );
	if( s != NULL )
	{
		f_bndry_tol = atof( (char*) s->data );
		fprintf( stderr, " ---> Loaded boundary tolerance of %15.7f\n", f_bndry_tol );
	}
	if( !b_neig_set )
	{
		s = symtab_lookup( &sym, "numberEigenvalues" );
		if( s != NULL )
		{
			neig = atoi( (char*) s->data );
			if( neig > 0 )
			{
				b_neig_set = 1;
				fprintf( stderr, " ---> Loaded number of eigenvalues to %d\n", neig );
			}
			else
				fprintf( stderr, " ---> ERROR: Invalid number of eigenvalues: %d\n", neig );
		}
	}
	s = symtab_lookup( &sym, "overlFileName" );
        if( s != NULL )
        {
                strncpy( s_overl_fname, (char*) s->data, TOKEN_BUFFER_LENGTH );
                fprintf( stderr, " ---> File name of overlap matrix set to \"%s\"\n", s_overl_fname );
        }
        s = symtab_lookup( &sym, "stiffFileName" );
        if( s != NULL )
        {
                strncpy( s_stiff_fname, (char*) s->data, TOKEN_BUFFER_LENGTH );
                fprintf( stderr, " ---> File name of stiffness matrix set to \"%s\"\n", s_stiff_fname );
        }
	s = symtab_lookup( &sym, "projFileName" );
	if( s != NULL )
	{
		strncpy( s_proj_fname, (char*) s->data, TOKEN_BUFFER_LENGTH );
		fprintf( stderr, " ---> File name of projection matrix set to \"%s\"\n", s_proj_fname );
	}
	s = symtab_lookup( &sym, "loadOverlMat" );
        if( s != NULL )
        {
                if( strcmp( (char*) s->data, "yes" ) == 0 )
                {
                        b_load_overl_mat = 1;
                        fprintf( stderr, " ---> Loading overlap matrix from file \"%s\"\n", s_overl_fname );
                }
        }
	s = symtab_lookup( &sym, "loadStiffMat" );
	if( s != NULL )
	{
		if( strcmp( (char*) s->data, "yes" ) == 0 )
		{
			b_load_stiff_mat = 1;
			fprintf( stderr, " ---> Loading stiffness matrix from file \"%s\"\n", s_stiff_fname );
		}
	}
	s = symtab_lookup( &sym, "saveOverlMat" );
        if( s != NULL )
        {
                if( strcmp( (char*) s->data, "yes" ) == 0 )
                {
                        b_save_overl_mat = 1;
                        fprintf( stderr, " ---> Saving overlap matrix to file \"%s\"\n", s_overl_fname );
                }
        }
        s = symtab_lookup( &sym, "saveStiffMat" );
        if( s != NULL )
        {
                if( strcmp( (char*) s->data, "yes" ) == 0 )
                {
                        b_save_stiff_mat = 1;
                        fprintf( stderr, " ---> Saving stiffness matrix to file \"%s\"\n", s_stiff_fname );
                }
        }
	s = symtab_lookup( &sym, "saveProjMat" );
	if( s != NULL )
	{
		if( strcmp( (char*) s->data, "yes" ) == 0 )
		{
			b_save_proj_mat = 1;
			fprintf( stderr, " ---> Saving projection matrix to file \"%s\"\n", s_proj_fname );
		}
	}
	s = symtab_lookup( &sym, "boundaryProject" );
	if( s != NULL )
	{
		if( strcmp( (char*) s->data, "yes" ) == 0 )
		{
			b_boundary_proj = 1;
			fprintf( stderr, " ---> Using projection method to enforce boundary conditions\n" );
		}
	}
	s = symtab_lookup( &sym, "singularityOrder" );
	if( s != NULL )
	{
		i_sing_order = atoi( (char*) s->data );
		fprintf( stderr, " ---> Using singularity order %d for singular particles\n", i_sing_order );
	}
	if( !b_treig_routine_set )
	{
		s = symtab_lookup( &sym, "treigRoutine" );
		if( s!= NULL )
		{
			strncpy( s_treig_routine, (char*) s->data, TOKEN_BUFFER_LENGTH );
			b_treig_routine_set = 1;
			fprintf( stderr, " ---> Using routine \"%s\" for tridiagonal eigenvalue solver\n", s_treig_routine );
		}
	}
	s = symtab_lookup( &sym, "maxLanczosSteps" );
	if( s != NULL )
	{
		i_max_lanczos_steps = atoi( (char*) s->data );
		fprintf( stderr, " ---> Using maximum of %d Lanczos steps\n", i_max_lanczos_steps );
	}
	s = symtab_lookup( &sym, "useSingular" );
	if( s != NULL )
	{
		if( strncmp( (char*) s->data, "no", TOKEN_BUFFER_LENGTH ) == 0 )
		{
			b_use_singular = 0;
			fprintf( stderr, " ---> Deactivating singular nodes\n" );
		}
		else
		{
			b_use_singular = 1;
			fprintf( stderr, " ---> Using singular nodes\n" );
		}
	}
	if( !b_points_loaded )
		fprintf( stderr, " ---> ERROR: Cannot load grid without loading points; need dimension\n" );
	else
	{
		dx = (int*) malloc( dim * sizeof(int) );
		x0 = (double*) malloc( dim * sizeof(double) );
		wx = (double*) malloc( dim * sizeof(double) );
		dd = (double*) malloc( dim * sizeof(double) );
		s = symtab_lookup( &sym, "gridDivision" );
		if( s != NULL )
		{
			strncpy( buf, (char*) s->data, TOKEN_BUFFER_LENGTH );
			n = parse_stokenize( buf, tok, "," );
			if( n < dim )
			{
				fprintf( stderr, " ---> ERROR: Too few entries for gridDivision\n" );
				for(i=0;i<n;i++)
					dx[i] = atoi( tok[i] );
				for(i=n;i<dim;i++)
					dx[i] = dx[n-1]; /* Fill it in with the last entry given */
			}
			else
				for(i=0;i<dim;i++)
					dx[i] = atoi( tok[i] );
			b_grid_set_dx = 1;
			fprintf( stderr, " ---> Set grid: dx = " );
			for(i=0;i<dim;i++)
				fprintf( stderr, "%5d", dx[i] );
			fprintf( stderr, "\n" );
		}
		s = symtab_lookup( &sym, "gridOrigin" );
		if( s != NULL )
		{
			strncpy( buf, (char*) s->data, TOKEN_BUFFER_LENGTH );
			n = parse_stokenize( buf, tok, "," );
			if( n < dim )
			{
				fprintf( stderr, " ---> ERROR: Too few entries for gridOrigin\n" );
				for(i=0;i<n;i++)
					x0[i] = atof( tok[i] );
				for(i=n;i<dim;i++)
					x0[i] = x0[n-1];
			}
			else
				for(i=0;i<dim;i++)
					x0[i] = atof( tok[i] );
			b_grid_set_x0 = 1;
			fprintf( stderr, " ---> Set grid: x0 = " );
			for(i=0;i<dim;i++)
				fprintf( stderr, "%12.5f", x0[i] );
			fprintf( stderr, "\n" );
		}
		s = symtab_lookup( &sym, "gridDimensions" );
		if( s != NULL )
		{
			strncpy( buf, (char*) s->data, TOKEN_BUFFER_LENGTH );
			n = parse_stokenize( buf, tok, "," );
			if( n < dim )
			{
				fprintf( stderr, " ---> ERROR: Too few entries for gridDimensions\n" );
				for(i=0;i<n;i++)
					wx[i] = atof( tok[i] );
				for(i=n;i<dim;i++)
					wx[i] = wx[n-1];
			}
			else
				for(i=0;i<dim;i++)
					wx[i] = atof( tok[i] );
			b_grid_set_wx = 1;
			fprintf( stderr, " ---> Set grid: wx = " );
			for(i=0;i<dim;i++)
				fprintf( stderr, "%12.5f", wx[i] );
			fprintf( stderr, "\n" );
		}
		for(i=0;i<dim;i++)
			dd[i] = wx[i] / (double) ( dx[i] - 1 );
	}
	s = symtab_lookup( &sym, "useExternalPotential" );
	if( s != NULL )
	{
		if( strncmp( (char*) s->data, "no", TOKEN_BUFFER_LENGTH ) == 0 )
		{
			fprintf( stderr, " ---> Turning off external potential for this run\n" );
			b_use_external_potential = 0;
		}
	}
	if( b_use_external_potential && !b_nuc_loaded )
	{
		s = symtab_lookup( &sym, "nucleiFile" );
		if( s != NULL )
		{
			ret = load_nuclei( (char*) s->data, &nuc, 0 );
			if( ret != 0 )
				fprintf( stderr, " ---> ERROR: Failed to load nuclei file from \"%s\"\n", (char*) s->data );
			else
			{
				fprintf( stderr, " ---> Loaded nuclear potential from file \"%s\"\n\n", (char*) s->data );
				nuclei_print( &nuc );
				b_nuc_loaded = 1;
			}
		}
	}
	
	fprintf( stderr, "\n" );
}

