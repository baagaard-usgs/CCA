/**
 * @file cca.c
 * @brief Main file for CCA library.
 *
 * @section DESCRIPTION
 *
 * Delivers the prototype CCA model which consists of En-Jui Lee's full 3D
 * tomographic results for central California.
 *
 */

#include "cca.h"

/**
 * Initializes the CCA plugin model within the UCVM framework. In order to initialize
 * the model, we must provide the UCVM install path and optionally a place in memory
 * where the model already exists.
 *
 * @param dir The directory in which UCVM has been installed.
 * @param label A unique identifier for the velocity model.
 * @return CCA_CODE_SUCCESS if initialization was successful, otherwise CCA_CODE_ERROR.
 */
int cca_init(const char *dir, const char *label) {
	int tempVal = 0;
	char configbuf[512];
	double north_height_m = 0, east_width_m = 0, rotation_angle = 0;

	// Initialize variables.
	cca_configuration = calloc(1, sizeof(cca_configuration_t));
	cca_velocity_model = calloc(1, sizeof(cca_model_t));
        cca_vs30_map = calloc(1, sizeof(cca_vs30_map_config_t));

	// Configuration file location.
	sprintf(configbuf, "%s/%s/config", dir, label);

        // Set up model directories.
        sprintf(cca_vs30_etree_file, "%s/ucvm/ucvm.e", dir);

	// Read the cca_configuration file.
	if (cca_read_configuration(configbuf, cca_configuration) != CCA_CODE_SUCCESS)
		return(CCA_CODE_ERROR);

	// Set up the iteration directory.
	sprintf(cca_iteration_directory, "%s/%s/%s", dir, label, cca_configuration->model_dir);

	// Can we allocate the model, or parts of it, to memory. If so, we do.
	tempVal = cca_try_reading_model(cca_velocity_model);

	if (tempVal == CCA_CODE_SUCCESS) {
		fprintf(stderr, "WARNING: Could not load model into memory. Reading the model from the\n");
		fprintf(stderr, "hard disk may result in slow performance.");
	} else if (tempVal == CCA_CODE_ERROR) {
		cca_print_error("No model file was found to read from.");
		return(CCA_CODE_ERROR);
	}

        if (cca_read_vs30_map(cca_vs30_etree_file, cca_vs30_map) != CCA_CODE_SUCCESS) {
                cca_print_error("Could not read the Vs30 map data from UCVM.");
                return(CCA_CODE_ERROR);
        }

	// We need to convert the point from lat, lon to UTM, let's set it up.
	if (!(cca_lonlat = pj_init_plus("+proj=lonlat +datum=WGS84"))) {
		cca_print_error("Could not set up latitude and longitude projection.");
		return(CCA_CODE_ERROR);
	}
	if (!(cca_utm = pj_init_plus("+proj=utm +zone=10 +datum=NAD27 +units=m +no_defs"))) {
		cca_print_error("Could not set up UTM projection.");
		return(CCA_CODE_ERROR);
	}
        if (!(cca_aeqd = pj_init_plus(cca_vs30_map->projection))) {
                cca_print_error("Could not set up AEQD projection.");
                return(CCA_CODE_ERROR);
        }

	// In order to simplify our calculations in the query, we want to rotate the box so that the bottom-left
	// corner is at (0m,0m). Our box's height is cca_total_height_m and cca_total_width_m. We then rotate the
	// point so that is is somewhere between (0,0) and (cca_total_width_m, cca_total_height_m). How far along
	// the X and Y axis determines which grid points we use for the interpolation routine.

	// Calculate the rotation angle of the box.
	north_height_m = cca_configuration->top_left_corner_n - cca_configuration->bottom_left_corner_n;
	east_width_m = cca_configuration->top_left_corner_e - cca_configuration->bottom_left_corner_e;

	// Rotation angle. Cos, sin, and tan are expensive computationally, so calculate once.
	rotation_angle = atan(east_width_m / north_height_m);

	cca_cos_rotation_angle = cos(rotation_angle);
	cca_sin_rotation_angle = sin(rotation_angle);

	cca_total_height_m = sqrt(pow(cca_configuration->top_left_corner_n - cca_configuration->bottom_left_corner_n, 2.0f) +
						  pow(cca_configuration->top_left_corner_e - cca_configuration->bottom_left_corner_e, 2.0f));
	cca_total_width_m  = sqrt(pow(cca_configuration->top_right_corner_n - cca_configuration->top_left_corner_n, 2.0f) +
						  pow(cca_configuration->top_right_corner_e - cca_configuration->top_left_corner_e, 2.0f));

        // Get the cos and sin for the Vs30 map rotation.
        cca_cos_vs30_rotation_angle = cos(cca_vs30_map->rotation * DEG_TO_RAD);
        cca_sin_vs30_rotation_angle = sin(cca_vs30_map->rotation * DEG_TO_RAD);

	// Let everyone know that we are initialized and ready for business.
	cca_is_initialized = 1;

	return CCA_CODE_SUCCESS;
}

/**
 * Called when the model is being discarded. Free all variables.
 *
 * @return CCA_CODE_SUCCESS
 */
int cca_finalize() {
	pj_free(cca_lonlat);
	pj_free(cca_utm);

	if (cca_velocity_model) free(cca_velocity_model);
	if (cca_configuration) free(cca_configuration);
        if (cca_vs30_map) free(cca_vs30_map);

	return CCA_CODE_SUCCESS;
}

/**
 * Returns the version information.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int cca_version(char *ver, int len)
{
  int verlen;
  verlen = strlen(cca_version_string);
  if (verlen > len - 1) {
    verlen = len - 1;
  }
  memset(ver, 0, len);
  strncpy(ver, cca_version_string, verlen);
  return 0;
}

/**
 * Set model parameter.
 * 
 * @param[in] name Name of parameter.
 * @param[in] value Value of parameter. 
 */
int cca_set_param(const char *name, const char *value)
{
  return(CCA_CODE_SUCCESS);
}

/**
 * Queries CCA at the given points and returns the data that it finds.
 *
 * @param points The points at which the queries will be made.
 * @param data The data that will be returned (Vp, Vs, density, Qs, and/or Qp).
 * @param numpoints The total number of points to query.
 * @return CCA_CODE_SUCCESS or CCA_CODE_ERROR.
 */
int cca_query(cca_point_t *points, cca_properties_t *data, int numpoints) {
	int i = 0, err=0;
	double point_u=0, point_v=0;
	double point_x = 0, point_y = 0;
	int load_x_coord = 0, load_y_coord = 0, load_z_coord = 0;
	double x_frac = 0, y_frac = 0, z_frac = 0;
	cca_properties_t surrounding_points[8];

	for (i = 0; i < numpoints; i++) {

		// We need to be below the surface to service this query.
		if (points[i].depth < 0) {
			data[i].vp = -1;
			data[i].vs = -1;
			data[i].rho = -1;
			data[i].qp = -1;
			data[i].qs = -1;
			continue;
		}

		point_u = points[i].longitude * DEG_TO_RAD;
		point_v = points[i].latitude * DEG_TO_RAD;
		err = pj_transform(cca_lonlat, cca_utm, 1, 1, &point_u, &point_v, NULL);
		if (err) {
			fprintf(stderr, "Error during coordinate transformation using Proj. %s\n", pj_strerrno(err));
			return(CCA_CODE_ERROR);
		}

		// Point within rectangle.
		point_u -= cca_configuration->bottom_left_corner_e;
		point_v -= cca_configuration->bottom_left_corner_n;

		// We need to rotate that point, the number of degrees we calculated above.
		point_x = cca_cos_rotation_angle * point_u - cca_sin_rotation_angle * point_v;
		point_y = cca_sin_rotation_angle * point_u + cca_cos_rotation_angle * point_v;

		// Which point base point does that correspond to?
		load_x_coord = floor(point_x / cca_total_width_m * (cca_configuration->nx - 1));
		load_y_coord = floor(point_y / cca_total_height_m * (cca_configuration->ny - 1));

		// And on the Z-axis?
		load_z_coord = (cca_configuration->depth / cca_configuration->depth_interval - 1) -
					   floor(points[i].depth / cca_configuration->depth_interval);

		// Are we outside the model's X and Y boundaries?
		if (load_x_coord > cca_configuration->nx - 2 || load_y_coord > cca_configuration->ny - 2 || load_x_coord < 0 || load_y_coord < 0) {
			data[i].vp = -1;
			data[i].vs = -1;
			data[i].rho = -1;
			data[i].qp = -1;
			data[i].qs = -1;
			continue;
		}

		// Get the X, Y, and Z fractions for the bilinear or cca_trilinear interpolation below.
	        double x_interval=(cca_configuration->nx > 1) ?
                     cca_total_width_m / (cca_configuration->nx-1) : cca_total_width_m;
                double y_interval=(cca_configuration->ny > 1) ?
                     cca_total_height_m / (cca_configuration->ny-1) : cca_total_height_m;


		x_frac = fmod(point_x, x_interval) / (x_interval);
		y_frac = fmod(point_y, y_interval) / (y_interval);
		z_frac = fmod(points[i].depth, cca_configuration->depth_interval) / cca_configuration->depth_interval;

		if (load_z_coord < 1) {
			// We're below the model boundaries. Bilinearly interpolate the bottom plane and use that value.
			data[i].vp = -1;
			data[i].vs = -1;
			data[i].rho = -1;
			data[i].qp = -1;
			data[i].qs = -1;

		        continue;
                } else {
if ((points[i].depth < cca_configuration->depth_interval) && (cca_configuration->gtl == 1)) {
                           cca_get_vs30_based_gtl(&(points[i]), &(data[i]));
                           data[i].rho=cca_calculate_density(data[i].vs);

                        } else {
			    // Read all the surrounding point properties.
			    cca_read_properties(load_x_coord,     load_y_coord,     load_z_coord,     &(surrounding_points[0]));	// Orgin.
			    cca_read_properties(load_x_coord + 1, load_y_coord,     load_z_coord,     &(surrounding_points[1]));	// Orgin + 1x
			    cca_read_properties(load_x_coord,     load_y_coord + 1, load_z_coord,     &(surrounding_points[2]));	// Orgin + 1y
			    cca_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord,     &(surrounding_points[3]));	// Orgin + x + y, forms top plane.
			    cca_read_properties(load_x_coord,     load_y_coord,     load_z_coord - 1, &(surrounding_points[4]));	// Bottom plane origin
			    cca_read_properties(load_x_coord + 1, load_y_coord,     load_z_coord - 1, &(surrounding_points[5]));	// +1x
			    cca_read_properties(load_x_coord,     load_y_coord + 1, load_z_coord - 1, &(surrounding_points[6]));	// +1y
			    cca_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[7]));	// +x +y, forms bottom plane.

			    cca_trilinear_interpolation(x_frac, y_frac, z_frac, surrounding_points, &(data[i]));
		}
}

		// Calculate Qp and Qs.
		if (data[i].vs < 1500)
			data[i].qs = data[i].vs * 0.02;
		else
			data[i].qs = data[i].vs * 0.10;

		data[i].qp = data[i].qs * 1.5;
	}

	return CCA_CODE_SUCCESS;
}

/**
 * Retrieves the material properties (whatever is available) for the given data point, expressed
 * in x, y, and z co-ordinates.
 *
 * @param x The x coordinate of the data point.
 * @param y The y coordinate of the data point.
 * @param z The z coordinate of the data point.
 * @param data The properties struct to which the material properties will be written.
 */
void cca_read_properties(int x, int y, int z, cca_properties_t *data) {
	// Set everything to -1 to indicate not found.
	data->vp = -1;
	data->vs = -1;
	data->rho = -1;
	data->qp = -1;
	data->qs = -1;
	float *ptr = NULL;
	FILE *fp = NULL;
	int location = z * cca_configuration->nx * cca_configuration->ny + y * cca_configuration->nx + x;

	// Check our loaded components of the model.
	if (cca_velocity_model->vs_status == 2) {
		// Read from memory.
		ptr = (float *)cca_velocity_model->vs;
		data->vs = ptr[location];
	} else if (cca_velocity_model->vs_status == 1) {
		// Read from file.
		fp = (FILE *)cca_velocity_model->vs;
		fseek(fp, location * sizeof(float), SEEK_SET);
		fread(&(data->vs), sizeof(float), 1, fp);
	}

	// Check our loaded components of the model.
	if (cca_velocity_model->vp_status == 2) {
		// Read from memory.
		ptr = (float *)cca_velocity_model->vp;
		data->vp = ptr[location];
	} else if (cca_velocity_model->vp_status == 1) {
		// Read from file.
		fseek(fp, location * sizeof(float), SEEK_SET);
		fread(&(data->vp), sizeof(float), 1, fp);
	}

	// Check our loaded components of the model.
	if (cca_velocity_model->rho_status == 2) {
		// Read from memory.
		ptr = (float *)cca_velocity_model->rho;
		data->rho = ptr[location];
	} else if (cca_velocity_model->rho_status == 1) {
		// Read from file.
		fseek(fp, location * sizeof(float), SEEK_SET);
		fread(&(data->rho), sizeof(float), 1, fp);
	}
}

/**
 * Trilinearly interpolates given a x fraction, y fraction, z fraction and a cube of
 * data properties in top origin format (top plane first, bottom plane second).
 *
 * @param x_frac X fraction
 * @param y_frac Y fraction
 * @param z_frac Z fraction
 * @param eight_points Eight surrounding data properties
 * @param ret_properties Returned data properties
 */
void cca_trilinear_interpolation(double x_frac, double y_frac, double z_frac,
							 cca_properties_t *eight_points, cca_properties_t *ret_properties) {
	cca_properties_t *temp_array = calloc(2, sizeof(cca_properties_t));
	cca_properties_t *four_points = eight_points;

	cca_bilinear_interpolation(x_frac, y_frac, four_points, &temp_array[0]);

	// Now advance the pointer four "cca_properties_t" spaces.
	four_points += 4;

	// Another interpolation.
	cca_bilinear_interpolation(x_frac, y_frac, four_points, &temp_array[1]);

	// Now linearly interpolate between the two.
	cca_linear_interpolation(z_frac, &temp_array[0], &temp_array[1], ret_properties);

	free(temp_array);
}

/**
 * Bilinearly interpolates given a x fraction, y fraction, and a plane of data properties in
 * origin, bottom-right, top-left, top-right format.
 *
 * @param x_frac X fraction.
 * @param y_frac Y fraction.
 * @param four_points Data property plane.
 * @param ret_properties Returned data properties.
 */
void cca_bilinear_interpolation(double x_frac, double y_frac, cca_properties_t *four_points, cca_properties_t *ret_properties) {
	cca_properties_t *temp_array = calloc(2, sizeof(cca_properties_t));
	cca_linear_interpolation(x_frac, &four_points[0], &four_points[1], &temp_array[0]);
	cca_linear_interpolation(x_frac, &four_points[2], &four_points[3], &temp_array[1]);
	cca_linear_interpolation(y_frac, &temp_array[0], &temp_array[1], ret_properties);
	free(temp_array);
}

/**
 * Linearly interpolates given a fraction from x0 to x1, a data point at x0, and a data point at x1.
 *
 * @param frac fraction of the way from x0 to x1 (from 0 to 1 interval).
 * @param x0 Data point at x0.
 * @param x1 Data point at x1.
 * @param ret_properties Resulting data properties.
 */
void cca_linear_interpolation(double frac, cca_properties_t *x0, cca_properties_t *x1, cca_properties_t *ret_properties) {
	ret_properties->vp  = (1.0 - frac) * x0->vp  + frac * x1->vp;
	ret_properties->vs  = (1.0 - frac) * x0->vs  + frac * x1->vs;
	ret_properties->rho = (1.0 - frac) * x0->rho + frac * x1->rho;
	ret_properties->qp  = (1.0 - frac) * x0->qp  + frac * x1->qp;
	ret_properties->qs  = (1.0 - frac) * x0->qs  + frac * x1->qs;
}

/**
 * Reads the cca_configuration file describing the various properties of CVM-S5 and populates
 * the cca_configuration struct. This assumes cca_configuration has been "calloc'ed" and validates
 * that each value is not zero at the end.
 *
 * @param file The cca_configuration file location on disk to read.
 * @param config The cca_configuration struct to which the data should be written.
 * @return CCA_CODE_SUCCESS or CCA_CODE_ERRORure, depending on if file was read CCA_CODE_SUCCESSfully.
 */
int cca_read_configuration(char *file, cca_configuration_t *config) {
	FILE *fp = fopen(file, "r");
	char key[40];
	char value[80];
	char line_holder[128];

	// If our file pointer is null, an error has occurred. Return(CCA_CODE_ERROR).
	if (fp == NULL) {
		cca_print_error("Could not open the cca_configuration file.");
		return(CCA_CODE_ERROR);
	}

	// Read the lines in the cca_configuration file.
	while (fgets(line_holder, sizeof(line_holder), fp) != NULL) {
		if (line_holder[0] != '#' && line_holder[0] != ' ' && line_holder[0] != '\n') {
			sscanf(line_holder, "%s = %s", key, value);

			// Which variable are we editing?
			if (strcmp(key, "utm_zone") == 0) 				config->utm_zone = atoi(value);
			if (strcmp(key, "model_dir") == 0)				sprintf(config->model_dir, "%s", value);
			if (strcmp(key, "nx") == 0) 	  				config->nx = atoi(value);
			if (strcmp(key, "ny") == 0) 	  			 	config->ny = atoi(value);
			if (strcmp(key, "nz") == 0) 	  			 	config->nz = atoi(value);
			if (strcmp(key, "depth") == 0) 	  			 	config->depth = atof(value);
			if (strcmp(key, "top_left_corner_e") == 0) 		config->top_left_corner_e = atof(value);
			if (strcmp(key, "top_left_corner_n") == 0)		config->top_left_corner_n = atof(value);
			if (strcmp(key, "top_right_corner_e") == 0)		config->top_right_corner_e = atof(value);
			if (strcmp(key, "top_right_corner_n") == 0)		config->top_right_corner_n = atof(value);
			if (strcmp(key, "bottom_left_corner_e") == 0)	config->bottom_left_corner_e = atof(value);
			if (strcmp(key, "bottom_left_corner_n") == 0)	config->bottom_left_corner_n = atof(value);
			if (strcmp(key, "bottom_right_corner_e") == 0)	config->bottom_right_corner_e = atof(value);
			if (strcmp(key, "bottom_right_corner_n") == 0)	config->bottom_right_corner_n = atof(value);
			if (strcmp(key, "depth_interval") == 0)			config->depth_interval = atof(value);
                        if (strcmp(key, "p0") == 0)                                             config->p0 = atof(value);
                        if (strcmp(key, "p1") == 0)                                             config->p1 = atof(value);
                        if (strcmp(key, "p2") == 0)                                             config->p2 = atof(value);
                        if (strcmp(key, "p3") == 0)                                             config->p3 = atof(value);
                        if (strcmp(key, "p4") == 0)                                             config->p4 = atof(value);
                        if (strcmp(key, "p5") == 0)                                             config->p5 = atof(value);
                        if (strcmp(key, "gtl") == 0) {
                                if (strcmp(value, "on") == 0) config->gtl = 1;
                                else config->gtl = 0;
                        }
// anything else, just ignore
		}
	}

	// Have we set up all cca_configuration parameters?
	if (config->utm_zone == 0 || config->nx == 0 || config->ny == 0 || config->nz == 0 || config->model_dir[0] == '\0' ||
		config->top_left_corner_e == 0 || config->top_left_corner_n == 0 || config->top_right_corner_e == 0 ||
		config->top_right_corner_n == 0 || config->bottom_left_corner_e == 0 || config->bottom_left_corner_n == 0 ||
		config->bottom_right_corner_e == 0 || config->bottom_right_corner_n == 0 || config->depth == 0 ||
		config->depth_interval == 0) {
		cca_print_error("One cca_configuration parameter not specified. Please check your cca_configuration file.");
		return(CCA_CODE_ERROR);
	}

	fclose(fp);

	return CCA_CODE_SUCCESS;
}

/**
 *  * Calculates the density based off of Vs. Based on Nafe-Drake scaling relationship.
 *   *
 *    * @param vs The Vs value off which to scale.
 *     * @return Density, in g/m^3.
 *      */
double cca_calculate_density(double vs) {
        double retVal;
        vs = vs / 1000;
        retVal = cca_configuration->p0 + cca_configuration->p1 * vs + cca_configuration->p2 * pow(vs, 2) +
                         cca_configuration->p3 * pow(vs, 3) + cca_configuration->p4 * pow(vs, 4) + cca_configuration->p5 * pow(vs, 5);
        retVal = retVal * 1000;
        return retVal;
}


/**
 * Prints the error string provided.
 *
 * @param err The error string to print out to stderr.
 */
void cca_print_error(char *err) {
	fprintf(stderr, "An error has occurred while executing CCA. The error was:\n\n");
	fprintf(stderr, "%s", err);
	fprintf(stderr, "\n\nPlease contact software@scec.org and describe both the error and a bit\n");
	fprintf(stderr, "about the computer you are running CCA on (Linux, Mac, etc.).\n");
}

/**
 * Tries to read the model into memory.
 *
 * @param model The model parameter struct which will hold the pointers to the data either on disk or in memory.
 * @return 2 if all files are read to memory, CCA_CODE_SUCCESS if file is found but at least 1
 * is not in memory, CCA_CODE_ERROR if no file found.
 */
int cca_try_reading_model(cca_model_t *model) {
	double base_malloc = cca_configuration->nx * cca_configuration->ny * cca_configuration->nz * sizeof(float);
	int file_count = 0;
	int all_read_to_memory = 1;
	char current_file[128];
	FILE *fp;

	// Let's see what data we actually have.
	sprintf(current_file, "%s/vp.dat", cca_iteration_directory);
	if (access(current_file, R_OK) == 0) {
		model->vp = malloc(base_malloc);
		if (model->vp != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->vp, 1, base_malloc, fp);
			fclose(fp);
			model->vp_status = 2;
		} else {
			all_read_to_memory = 0;
			model->vp = fopen(current_file, "rb");
			model->vp_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/vs.dat", cca_iteration_directory);
	if (access(current_file, R_OK) == 0) {
		model->vs = malloc(base_malloc);
		if (model->vs != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->vs, 1, base_malloc, fp);
			fclose(fp);
			model->vs_status = 2;
		} else {
			all_read_to_memory = 0;
			model->vs = fopen(current_file, "rb");
			model->vs_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/density.dat", cca_iteration_directory);
	if (access(current_file, R_OK) == 0) {
		model->rho = malloc(base_malloc);
		if (model->rho != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->rho, 1, base_malloc, fp);
			fclose(fp);
			model->rho_status = 2;
		} else {
			all_read_to_memory = 0;
			model->rho = fopen(current_file, "rb");
			model->rho_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/qp.dat", cca_iteration_directory);
	if (access(current_file, R_OK) == 0) {
		model->qp = malloc(base_malloc);
		if (model->qp != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->qp, 1, base_malloc, fp);
			fclose(fp);
			model->qp_status = 2;
		} else {
			all_read_to_memory = 0;
			model->qp = fopen(current_file, "rb");
			model->qp_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/qs.dat", cca_iteration_directory);
	if (access(current_file, R_OK) == 0) {
		model->qs = malloc(base_malloc);
		if (model->qs != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->qs, 1, base_malloc, fp);
			fclose(fp);
			model->qs_status = 2;
		} else {
			all_read_to_memory = 0;
			model->qs = fopen(current_file, "rb");
			model->qs_status = 1;
		}
		file_count++;
	}

	if (file_count == 0)
		return(CCA_CODE_ERROR);
	else if (file_count > 0 && all_read_to_memory == 0)
		return CCA_CODE_SUCCESS;
	else
		return 2;
}


/**
 * Reads the format of the Vs30 data e-tree. This file location is typically specified
 * in the cca_configuration file of the model.
 *
 * @param filename The e-tree's file location from which to read.
 * @param map The outputted map cca_configuration structure.
 */
int cca_read_vs30_map(char *filename, cca_vs30_map_config_t *map) {
	char appmeta[512];
	char *token;
	int index = 0, retVal = 0;
	map->vs30_map = etree_open(filename, O_RDONLY, 64, 0, 3);
	retVal = snprintf(appmeta, sizeof(appmeta), "%s", etree_getappmeta(map->vs30_map));

	if (retVal >= 0 && retVal < 128) {
		return(CCA_CODE_ERROR);
	}

	// Now we need to parse the map cca_configuration.
	index = 0;
	token = strtok(appmeta, "|");

	while (token != NULL) {
		switch (index) {
	    case 0:
	    	snprintf(map->type, sizeof(map->type), "%s", token);
	    	break;
	    case 1:
	    	snprintf(map->description, sizeof(map->description), "%s", token);
	    	break;
	    case 2:
	    	snprintf(map->author, sizeof(map->author), "%s", token);
	    	break;
	    case 3:
	    	snprintf(map->date, sizeof(map->date), "%s", token);
	    	break;
	    case 4:
	    	sscanf(token, "%lf", &(map->spacing));
	    	break;
	    case 5:
	    	snprintf(map->schema, sizeof(map->schema), "%s", token);
	    	break;
	    case 6:
	    	snprintf(map->projection, sizeof(map->projection), "%s", token);
	    	break;
	    case 7:
	    	sscanf(token, "%lf,%lf,%lf", &(map->origin_point.longitude), &(map->origin_point.latitude),
	    			&(map->origin_point.depth));
	    	break;
	    case 8:
	    	sscanf(token, "%lf", &(map->rotation));
	    	break;
	    case 9:
	    	sscanf(token, "%lf,%lf,%lf", &(map->x_dimension), &(map->y_dimension), &(map->z_dimension));
	    	break;
	    case 10:
	    	sscanf(token, "%u,%u,%u", &(map->x_ticks), &(map->y_ticks), &(map->z_ticks));
	    	break;
	    default:
	    	fprintf(stderr, "Unexpected metadata. Please check your Vs30 e-tree within UCVM.\n");
	    	return(CCA_CODE_ERROR);
	    	break;
		}
	    index++;
	    token = strtok(NULL, "|");
	}

	return(CCA_CODE_SUCCESS);

}

/**
 * Given a latitude and longitude in WGS84 co-ordinates, we find the corresponding e-tree octant
 * in the Vs30 map e-tree and read the value as well as interpolate bilinearly.
 *
 * @param longitude The longitude in WGS84 format.
 * @param latitude The latitude in WGS84 format.
 * @param map The Vs30 map structure as defined during the initialization procedure.
 * @return The Vs30 value at that point, or -1 if outside the boundaries.
 */
double cca_get_vs30_value(double longitude, double latitude, cca_vs30_map_config_t *map) {
	// Convert both points to UTM.
	double longitude_utm_e = longitude * DEG_TO_RAD;
	double latitude_utm_n = latitude * DEG_TO_RAD;
	double vs30_long_utm_e = map->origin_point.longitude * DEG_TO_RAD;
	double vs30_lat_utm_n = map->origin_point.latitude * DEG_TO_RAD;
	double temp_rotated_point_n = 0.0, temp_rotated_point_e = 0.0;
	double rotated_point_n = 0.0, rotated_point_e = 0.0;
	double percent = 0.0;
	int loc_x = 0, loc_y = 0;
	etree_addr_t addr;
	cca_vs30_mpayload_t vs30_payload[4];

	int max_level = ceil(log(map->x_dimension / map->spacing) / log(2.0));

	etree_tick_t edgetics = (etree_tick_t)1 << (ETREE_MAXLEVEL - max_level);
	double map_edgesize = map->x_dimension / (double)((etree_tick_t)1<<max_level);

	pj_transform(cca_lonlat, cca_aeqd, 1, 1, &longitude_utm_e, &latitude_utm_n, NULL);
	pj_transform(cca_lonlat, cca_aeqd, 1, 1, &vs30_long_utm_e, &vs30_lat_utm_n, NULL);

	// Now that both are in UTM, we can subtract and rotate.
	temp_rotated_point_e = longitude_utm_e - vs30_long_utm_e;
	temp_rotated_point_n = latitude_utm_n - vs30_lat_utm_n;

	rotated_point_e = cca_cos_vs30_rotation_angle * temp_rotated_point_e - cca_sin_vs30_rotation_angle * temp_rotated_point_n;
	rotated_point_n = cca_sin_vs30_rotation_angle * temp_rotated_point_e + cca_cos_vs30_rotation_angle * temp_rotated_point_n;

	// Are we within the box?
	if (rotated_point_e < 0 || rotated_point_n < 0 || rotated_point_e > map->x_dimension ||
		rotated_point_n > map->y_dimension) return -1;

	// Get the integer location of the grid point within the map.
	loc_x = floor(rotated_point_e / map_edgesize);
	loc_y = floor(rotated_point_n / map_edgesize);

	// We need the four surrounding points for bilinear interpolation.
	addr.level = ETREE_MAXLEVEL;
	addr.x = loc_x * edgetics; addr.y = loc_y * edgetics; addr.z = 0;
    /* Adjust addresses for edges of grid */
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[0]));
	addr.x = (loc_x + 1) * edgetics; addr.y = loc_y * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[1]));
	addr.x = loc_x * edgetics; addr.y = (loc_y + 1) * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[2]));
	addr.x = (loc_x + 1) * edgetics; addr.y = (loc_y + 1) * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[3]));

	percent = fmod(rotated_point_e / map->spacing, map->spacing) / map->spacing;
	vs30_payload[0].vs30 = percent * vs30_payload[0].vs30 + (1 - percent) * vs30_payload[1].vs30;
	vs30_payload[1].vs30 = percent * vs30_payload[2].vs30 + (1 - percent) * vs30_payload[3].vs30;

	return vs30_payload[0].vs30;
}

/**
 * Gets the GTL value using the Wills and Wald dataset, given a latitude, longitude and depth.
 *
 * @param point The point at which to retrieve the property. Note, depth is ignored.
 * @param data The material properties at the point specified, or -1 if not found.
 * @return CCA_CODE_SUCCESS or CCA_CODE_ERRORure.
 */
int cca_get_vs30_based_gtl(cca_point_t *point, cca_properties_t *data) {
        double a = 0.5, b = 0.6, c = 0.5;
	double percent_z = point->depth / cca_configuration->depth_interval;
	double f = 0.0, g = 0.0;
	double vs30 = 0.0, vp30 = 0.0;

	// Double check that we're above the first layer.
	if (percent_z > 1) return(CCA_CODE_ERROR);

	// Query for the point at depth_interval.
	cca_point_t *pt = calloc(1, sizeof(cca_point_t));
	cca_properties_t *dt = calloc(1, sizeof(cca_properties_t));

	pt->latitude = point->latitude;
	pt->longitude = point->longitude;
	pt->depth = cca_configuration->depth_interval;

	if (cca_query(pt, dt, 1) != CCA_CODE_SUCCESS) return(CCA_CODE_ERROR);

	// Now we need the Vs30 data value.
	vs30 = cca_get_vs30_value(point->longitude, point->latitude, cca_vs30_map);

	if (vs30 == -1) {
		data->vp = -1;
		data->vs = -1;
	} else {
		// Get the point's material properties within the GTL.
		f = percent_z + b * (percent_z - pow(percent_z, 2.0f));
		g = a - a * percent_z + c * (pow(percent_z, 2.0f) + 2.0 * sqrt(percent_z) - 3.0 * percent_z);
		data->vs = f * dt->vs + g * vs30;
//fprintf(stderr,"XXX f %f and g %f\n", f, g);
		vs30 = vs30 / 1000;
		vp30 = 0.9409 + 2.0947 * vs30 - 0.8206 * pow(vs30, 2.0f) + 0.2683 * pow(vs30, 3.0f) - 0.0251 * pow(vs30, 4.0f);
		vp30 = vp30 * 1000;
		data->vp = f * dt->vp + g * vp30;
	}

	free(pt);
	free(dt);

	return(CCA_CODE_SUCCESS);
}


// The following functions are for dynamic library mode. If we are compiling
// a static library, these functions must be disabled to avoid conflicts.
#if defined(BUILD_SHARED_LIBRARY)

/**
 * Init function loaded and called by the UCVM library. Calls cca_init.
 *
 * @param dir The directory in which UCVM is installed.
 * @return CCA_CODE_SUCCESS or CCA_CODE_ERRORure.
 */
int model_init(const char *dir, const char *label) {
	return cca_init(dir, label);
}

/**
 * Finalize function loaded and called by the UCVM library. Calls cca_finalize.
 *
 * @return CCA_CODE_SUCCESS
 */
int model_finalize() {
	return cca_finalize();
}

/**
 * Version function loaded and called by the UCVM library. Calls cca_version.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int model_version(char *ver, int len) {
	return cca_version(ver, len);
}

/**
 * Set model parameter.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int model_set_param(const char *name, const char *value) {
	return cca_set_param(name, value);
}

/**
 * Query function loaded and called by the UCVM library. Calls cca_query.
 *
 * @param points The basic_point_t array containing the points.
 * @param data The basic_properties_t array containing the material properties returned.
 * @param numpoints The number of points in the array.
 * @return CCA_CODE_SUCCESS or CCA_CODE_ERROR.
 */
int model_query(cca_point_t *points, cca_properties_t *data, int numpoints) {
	return cca_query(points, data, numpoints);
}

int (*get_model_init())(const char *, const char *) {
        return &cca_init;
}
int (*get_model_finalize())() {
         return &cca_finalize;
}
int (*get_model_version())(char *, int) {
         return &cca_version;
}
int (*get_model_set_param())(const char *, const char*) {
         return &cca_set_param;
}
int (*get_model_query())(cca_point_t *, cca_properties_t *, int) {
         return &cca_query;
}

#endif
