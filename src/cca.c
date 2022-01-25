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

#include "proj.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif
#define DEG_TO_RAD M_PI / 180.0

const char *cca_version_string = "CCA";

// Variables

/** Location of the ucvm.e e-tree file. */
char cca_vs30_etree_file[256];
/** Location of En-Jui's latest iteration files. */
char cca_iteration_directory[256];

/** Configuration parameters. */
cca_configuration_t *cca_configuration = NULL;

/** Holds the configuration parameters for the Vs30 map. */
cca_vs30_map_config_t *cca_vs30_map = NULL;

/** Holds pointers to the velocity model data OR indicates it can be read from file. */
cca_model_t *cca_velocity_model = NULL;

/** Proj coordinate transformation objects. */
PJ *cca_geo2utm = NULL;
PJ *cca_geo2aeqd = NULL;

/** The cosine of the rotation angle used to rotate the box and point around the bottom-left corner. */
double cca_cos_rotation_angle = 0;
/** The sine of the rotation angle used to rotate the box and point around the bottom-left corner. */
double cca_sin_rotation_angle = 0;

/** The height of this model's region, in meters. */
double cca_total_height_m = 0;
/** The width of this model's region, in meters. */
double cca_total_width_m = 0;

/** The cosine of the Vs30 map's rotation. */
double cca_cos_vs30_rotation_angle = 0;
/** The sine of the Vs30 map's rotation. */
double cca_sin_vs30_rotation_angle = 0;

/**
 * Creates the CCA plugin model within the UCVM framework. In order to initialize
 * the model, we must provide the UCVM install path and optionally a place in memory
 * where the model already exists.
 *
 * @param models_dir The directory containing the UCVM models ($datarootdir/$package/model).
 * @param label A unique identifier for the velocity model.
 * @return CCA_CODE_SUCCESS if initialization was successful, otherwise CCA_CODE_ERROR.
 */
int
cca_create(const char *models_dir,
           const char *label) {
    char config_filename[512];
    int err = CCA_CODE_SUCCESS;

    // Read the cca_configuration file.
    cca_configuration = calloc(1, sizeof(cca_configuration_t));
    sprintf(config_filename, "%s/%s/config", models_dir, label);
    if (cca_read_configuration(config_filename, cca_configuration) != CCA_CODE_SUCCESS) {
        return (CCA_CODE_ERROR);
    }

    // Can we allocate the model, or parts of it, to memory. If so, we do.
    sprintf(cca_iteration_directory, "%s/%s/%s", models_dir, label, cca_configuration->model_dir);
    cca_velocity_model = calloc(1, sizeof(cca_model_t));
    err = cca_try_reading_model(cca_velocity_model);
    if (err == CCA_CODE_SUCCESS) {
        fprintf(stderr, "WARNING: Could not load model into memory. Reading the model from the\n");
        fprintf(stderr, "hard disk may result in slow performance.");
    } else if (err == CCA_CODE_ERROR) {
        cca_print_error("No model file was found to read from.");
        return (CCA_CODE_ERROR);
    }

    // Read Vs30 model.
    sprintf(cca_vs30_etree_file, "%s/ucvm/ucvm.e", models_dir);
    cca_vs30_map = calloc(1, sizeof(cca_vs30_map_config_t));
    if (cca_read_vs30_map(cca_vs30_etree_file, cca_vs30_map) != CCA_CODE_SUCCESS) {
        cca_print_error("Could not read the Vs30 map data from UCVM.");
        return (CCA_CODE_ERROR);
    }

    return CCA_CODE_SUCCESS;
}


/**
 * Initializes the CCA plugin model within the UCVM framework.
 *
 * @return CCA_CODE_SUCCESS if initialization was successful, otherwise CCA_CODE_ERROR.
 */
int
cca_initialize(void) {
    char cca_projstr[64];
    double north_height_m = 0, east_width_m = 0, rotation_angle = 0;

    // We need to convert the point from lat, lon to UTM, let's set it up.
    snprintf(cca_projstr, 64, "+proj=utm +zone=%d +datum=NAD27 +units=m +no_defs", cca_configuration->utm_zone);
    if (!(cca_geo2utm = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326", cca_projstr, NULL))) {
        cca_print_error("Could not set up Proj transformation from EPSG:4325 to UTM.");
        cca_print_error(proj_context_errno_string(PJ_DEFAULT_CTX, proj_context_errno(PJ_DEFAULT_CTX)));
        return (CCA_CODE_ERROR);
    }
    assert(cca_vs30_map);
    if (!(cca_geo2aeqd = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326", cca_vs30_map->projection, NULL))) {
        cca_print_error("Could not set up Proj transformation from EPSG:4326 to AEQD projection.");
        cca_print_error(proj_context_errno_string(PJ_DEFAULT_CTX, proj_context_errno(PJ_DEFAULT_CTX)));
        return (CCA_CODE_ERROR);
    }

    // In order to simplify our calculations in the query, we want to rotate the box so that the bottom-left
    // corner is at (0m,0m). Our box's height is cca_total_height_m and cca_total_width_m. We then rotate the
    // point so that is is somewhere between (0,0) and (cca_total_width_m, cca_total_height_m). How far along
    // the X and Y axis determines which grid points we use for the interpolation routine.

    // Calculate the rotation angle of the box.
    assert(cca_configuration);
    north_height_m = cca_configuration->top_left_corner_n - cca_configuration->bottom_left_corner_n;
    east_width_m = cca_configuration->top_left_corner_e - cca_configuration->bottom_left_corner_e;

    // Rotation angle. Cos, sin, and tan are expensive computationally, so calculate once.
    rotation_angle = atan(east_width_m / north_height_m);

    cca_cos_rotation_angle = cos(rotation_angle);
    cca_sin_rotation_angle = sin(rotation_angle);

    cca_total_height_m = sqrt(pow(cca_configuration->top_left_corner_n - cca_configuration->bottom_left_corner_n, 2.0f) +
                              pow(cca_configuration->top_left_corner_e - cca_configuration->bottom_left_corner_e, 2.0f));
    cca_total_width_m = sqrt(pow(cca_configuration->top_right_corner_n - cca_configuration->top_left_corner_n, 2.0f) +
                             pow(cca_configuration->top_right_corner_e - cca_configuration->top_left_corner_e, 2.0f));

    // Get the cos and sin for the Vs30 map rotation.
    cca_cos_vs30_rotation_angle = cos(cca_vs30_map->rotation * DEG_TO_RAD);
    cca_sin_vs30_rotation_angle = sin(cca_vs30_map->rotation * DEG_TO_RAD);

    return CCA_CODE_SUCCESS;
}


/**
 * Called when the model is being discarded. Free all variables.
 *
 * @return CCA_CODE_SUCCESS on success or CCA_CODE_ERROR on failure.
 */
int
cca_finalize(void) {
    proj_destroy(cca_geo2utm);cca_geo2utm = NULL;
    proj_destroy(cca_geo2aeqd);cca_geo2aeqd = NULL;

    if (cca_velocity_model) { free(cca_velocity_model);cca_velocity_model = NULL; }
    if (cca_configuration) { free(cca_configuration);cca_configuration = NULL; }
    if (cca_vs30_map) { free(cca_vs30_map);cca_vs30_map = NULL; }

    return CCA_CODE_SUCCESS;
}


/**
 * Returns the version information.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return CCA_CODE_SUCCESS on success or CCA_CODE_ERROR on failure.
 */
int
cca_version(char *ver,
            int len) {
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
 * @return CCA_CODE_SUCCESS on success or CCA_CODE_ERROR on failure.
 */
int
cca_set_parameter(const char *name,
                  const char *value) {
    if (strcasecmp(name, "use_gtl") == 0) {
        cca_configuration->use_gtl = strcasecmp(value, "true") == 0 ? 1 : 0;
    } else {
        fprintf(stderr, "Unknown parameter %s=%s for CCA model.\n", name, value);
        return CCA_CODE_ERROR;
    }

    return CCA_CODE_SUCCESS;
}


/**
 * Queries CCA at the given points and returns the data that it finds.
 *
 * @param points The points at which the queries will be made.
 * @param data The data that will be returned (Vp, Vs, density, Qs, and/or Qp).
 * @param numpoints The total number of points to query.
 * @return CCA_CODE_SUCCESS on success or CCA_CODE_ERROR on failure.
 */
int
cca_query(cca_point_t *points,
          cca_properties_t *data,
          int numpoints,
          cca_query_flags_t *qflags) {
    int i = 0, err = CCA_CODE_SUCCESS;
    double point_u = 0, point_v = 0;
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

        PJ_COORD xyzSrc = proj_coord(points[i].latitude, points[i].longitude, 0.0, HUGE_VAL);
        PJ_COORD xyzDest = proj_trans(cca_geo2utm, PJ_FWD, xyzSrc);
        err = proj_context_errno(PJ_DEFAULT_CTX);
        if (err) {
            fprintf(stderr, "Error occurred while transforming latitude=%.4f, longitude=%.4f to UTM.\n",
                    points[i].latitude, points[i].longitude);
            fprintf(stderr, "Proj error: %s\n", proj_context_errno_string(PJ_DEFAULT_CTX, err));
            return CCA_CODE_ERROR;
        }

        point_u = xyzDest.xyzt.x;
        point_v = xyzDest.xyzt.y;

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
        if ((load_x_coord > cca_configuration->nx - 2) || (load_y_coord > cca_configuration->ny - 2) || (load_x_coord < 0) || (load_y_coord < 0)) {
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;
            continue;
        }

        // Get the X, Y, and Z fractions for the bilinear or cca_trilinear interpolation below.
        double x_interval = (cca_configuration->nx > 1) ?
                            cca_total_width_m / (cca_configuration->nx-1) : cca_total_width_m;
        double y_interval = (cca_configuration->ny > 1) ?
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
            if ((points[i].depth < cca_configuration->depth_interval) && cca_configuration->use_gtl) {
                cca_get_vs30_based_gtl(&(points[i]), &(data[i]));
                data[i].rho = cca_calculate_density(data[i].vs);

            } else {
                // Read all the surrounding point properties.

                /* Top plane */
                cca_read_properties(load_x_coord,     load_y_coord,     load_z_coord,     &(surrounding_points[0]));
                cca_read_properties(load_x_coord + 1, load_y_coord,     load_z_coord,     &(surrounding_points[1]));
                cca_read_properties(load_x_coord,     load_y_coord + 1, load_z_coord,     &(surrounding_points[2]));
                cca_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord,     &(surrounding_points[3]));

                /* Bottom plane */
                cca_read_properties(load_x_coord,     load_y_coord,     load_z_coord - 1, &(surrounding_points[4]));
                cca_read_properties(load_x_coord + 1, load_y_coord,     load_z_coord - 1, &(surrounding_points[5]));
                cca_read_properties(load_x_coord,     load_y_coord + 1, load_z_coord - 1, &(surrounding_points[6]));
                cca_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[7]));

                cca_trilinear_interpolation(x_frac, y_frac, z_frac, surrounding_points, &(data[i]));
            }
        }

        // Calculate Qp and Qs.
        if (data[i].vs < 1500) {
            data[i].qs = data[i].vs * 0.02;
        } else {
            data[i].qs = data[i].vs * 0.10;
        }

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
void
cca_read_properties(int x,
                    int y,
                    int z,
                    cca_properties_t *data) {
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
void
cca_trilinear_interpolation(double x_frac,
                            double y_frac,
                            double z_frac,
                            cca_properties_t *eight_points,
                            cca_properties_t *ret_properties) {
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
void
cca_bilinear_interpolation(double x_frac,
                           double y_frac,
                           cca_properties_t *four_points,
                           cca_properties_t *ret_properties) {
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
void
cca_linear_interpolation(double frac,
                         cca_properties_t *x0,
                         cca_properties_t *x1,
                         cca_properties_t *ret_properties) {
    ret_properties->vp = (1.0 - frac) * x0->vp  + frac * x1->vp;
    ret_properties->vs = (1.0 - frac) * x0->vs  + frac * x1->vs;
    ret_properties->rho = (1.0 - frac) * x0->rho + frac * x1->rho;
    ret_properties->qp = (1.0 - frac) * x0->qp  + frac * x1->qp;
    ret_properties->qs = (1.0 - frac) * x0->qs  + frac * x1->qs;
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
int
cca_read_configuration(char *file,
                       cca_configuration_t *config) {
    FILE *fp = fopen(file, "r");
    char key[40];
    char value[80];
    char line_holder[128];

    // If our file pointer is null, an error has occurred. Return(CCA_CODE_ERROR).
    if (fp == NULL) {
        cca_print_error("Could not open the cca_configuration file.");
        return (CCA_CODE_ERROR);
    }

    // Read the lines in the cca_configuration file.
    while (fgets(line_holder, sizeof(line_holder), fp) != NULL) {
        if ((line_holder[0] != '#') && (line_holder[0] != ' ') && (line_holder[0] != '\n')) {
            sscanf(line_holder, "%s = %s", key, value);

            if (strcmp(key, "utm_zone") == 0) {config->utm_zone = atoi(value);}
            if (strcmp(key, "model_dir") == 0) {sprintf(config->model_dir, "%s", value);}
            if (strcmp(key, "nx") == 0) {config->nx = atoi(value);}
            if (strcmp(key, "ny") == 0) {config->ny = atoi(value);}
            if (strcmp(key, "nz") == 0) {config->nz = atoi(value);}
            if (strcmp(key, "depth") == 0) {config->depth = atof(value);}
            if (strcmp(key, "top_left_corner_e") == 0) {config->top_left_corner_e = atof(value);}
            if (strcmp(key, "top_left_corner_n") == 0) {config->top_left_corner_n = atof(value);}
            if (strcmp(key, "top_right_corner_e") == 0) {config->top_right_corner_e = atof(value);}
            if (strcmp(key, "top_right_corner_n") == 0) {config->top_right_corner_n = atof(value);}
            if (strcmp(key, "bottom_left_corner_e") == 0) {config->bottom_left_corner_e = atof(value);}
            if (strcmp(key, "bottom_left_corner_n") == 0) {config->bottom_left_corner_n = atof(value);}
            if (strcmp(key, "bottom_right_corner_e") == 0) {config->bottom_right_corner_e = atof(value);}
            if (strcmp(key, "bottom_right_corner_n") == 0) {config->bottom_right_corner_n = atof(value);}
            if (strcmp(key, "depth_interval") == 0) {config->depth_interval = atof(value);}
            if (strcmp(key, "p0") == 0) { config->p0 = atof(value);}
            if (strcmp(key, "p1") == 0) { config->p1 = atof(value);}
            if (strcmp(key, "p2") == 0) { config->p2 = atof(value);}
            if (strcmp(key, "p3") == 0) { config->p3 = atof(value);}
            if (strcmp(key, "p4") == 0) { config->p4 = atof(value);}
            if (strcmp(key, "p5") == 0) { config->p5 = atof(value);}
            if (strcmp(key, "use_gtl") == 0) {
                config->use_gtl = strcasecmp(value, "true") == 0 ? 1 : 0;
            }
            // anything else, just ignore
        }
    }

    // Have we set up all cca_configuration parameters?
    if ((config->utm_zone == 0) || (config->nx == 0) || (config->ny == 0) || (config->nz == 0) || (config->model_dir[0] == '\0') ||
        (config->top_left_corner_e == 0) || (config->top_left_corner_n == 0) || (config->top_right_corner_e == 0) ||
        (config->top_right_corner_n == 0) || (config->bottom_left_corner_e == 0) || (config->bottom_left_corner_n == 0) ||
        (config->bottom_right_corner_e == 0) || (config->bottom_right_corner_n == 0) || (config->depth == 0) ||
        (config->depth_interval == 0)) {
        cca_print_error("One cca_configuration parameter not specified. Please check your cca_configuration file.");
        return (CCA_CODE_ERROR);
    }

    fclose(fp);

    return CCA_CODE_SUCCESS;
}


/**
 * Calculates the density based off of Vs. Based on Nafe-Drake scaling relationship.
 *
 * @param vs The Vs value off which to scale.
 * @return Density, in g/m^3.
 */
double
cca_calculate_density(double vs) {
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
void
cca_print_error(const char *err) {
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
int
cca_try_reading_model(cca_model_t *model) {
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

    if (file_count == 0) {
        return (CCA_CODE_ERROR);
    } else if ((file_count > 0) && (all_read_to_memory == 0)) {
        return CCA_CODE_SUCCESS;
    } else {
        return 2;
    }
}


/**
 * Reads the format of the Vs30 data e-tree. This file location is typically specified
 * in the cca_configuration file of the model.
 *
 * @param filename The e-tree's file location from which to read.
 * @param map The outputted map cca_configuration structure.
 */
int
cca_read_vs30_map(char *filename,
                  cca_vs30_map_config_t *map) {
    char appmeta[512];
    char *token;
    int index = 0, retVal = 0;
    map->vs30_map = etree_open(filename, O_RDONLY, 64, 0, 3);
    retVal = snprintf(appmeta, sizeof(appmeta), "%s", etree_getappmeta(map->vs30_map));

    if ((retVal >= 0) && (retVal < 128)) {
        return (CCA_CODE_ERROR);
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
            return (CCA_CODE_ERROR);
            break;
        }
        index++;
        token = strtok(NULL, "|");
    }

    return CCA_CODE_SUCCESS;

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
double
cca_get_vs30_value(double longitude,
                   double latitude,
                   cca_vs30_map_config_t *map) {
    // Convert both points to UTM.
    double point_x, point_y, origin_x, origin_y;
    double temp_rotated_point_x = 0.0, temp_rotated_point_y = 0.0;
    double rotated_point_x = 0.0, rotated_point_y = 0.0;
    double percent = 0.0;
    int loc_x = 0, loc_y = 0;
    etree_addr_t addr;
    cca_vs30_mpayload_t vs30_payload[4];

    int max_level = ceil(log(map->x_dimension / map->spacing) / log(2.0));

    etree_tick_t edgetics = (etree_tick_t)1 << (ETREE_MAXLEVEL - max_level);
    double map_edgesize = map->x_dimension / (double)((etree_tick_t)1<<max_level);

    PJ_COORD xyzSrc = proj_coord(latitude, longitude, 0.0, HUGE_VAL);
    PJ_COORD xyzDest = proj_trans(cca_geo2aeqd, PJ_FWD, xyzSrc);
    point_x = xyzDest.xyzt.x;
    point_y = xyzDest.xyzt.y;

    xyzSrc = proj_coord(map->origin_point.latitude, map->origin_point.longitude, 0.0, HUGE_VAL);
    xyzDest = proj_trans(cca_geo2aeqd, PJ_FWD, xyzSrc);
    origin_x = xyzDest.xyzt.x;
    origin_y = xyzDest.xyzt.y;

    // Now that both are in UTM, we can subtract and rotate.
    temp_rotated_point_x = point_x - origin_x;
    temp_rotated_point_y = point_y - origin_y;

    rotated_point_x = cca_cos_vs30_rotation_angle * temp_rotated_point_x - cca_sin_vs30_rotation_angle * temp_rotated_point_y;
    rotated_point_y = cca_sin_vs30_rotation_angle * temp_rotated_point_x + cca_cos_vs30_rotation_angle * temp_rotated_point_y;

    // Are we within the box?
    if ((rotated_point_x < 0) || (rotated_point_y < 0) || (rotated_point_x > map->x_dimension) ||
        (rotated_point_y > map->y_dimension) ) { return -1;}

    // Get the integer location of the grid point within the map.
    loc_x = floor(rotated_point_x / map_edgesize);
    loc_y = floor(rotated_point_y / map_edgesize);

    // We need the four surrounding points for bilinear interpolation.
    addr.level = ETREE_MAXLEVEL;
    addr.x = loc_x * edgetics;addr.y = loc_y * edgetics;addr.z = 0;
    /* Adjust addresses for edges of grid */
    if (addr.x >= map->x_ticks) {addr.x = map->x_ticks - edgetics;}
    if (addr.y >= map->y_ticks) {addr.y = map->y_ticks - edgetics;}
    etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[0]));
    addr.x = (loc_x + 1) * edgetics;addr.y = loc_y * edgetics;
    if (addr.x >= map->x_ticks) {addr.x = map->x_ticks - edgetics;}
    if (addr.y >= map->y_ticks) {addr.y = map->y_ticks - edgetics;}
    etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[1]));
    addr.x = loc_x * edgetics;addr.y = (loc_y + 1) * edgetics;
    if (addr.x >= map->x_ticks) {addr.x = map->x_ticks - edgetics;}
    if (addr.y >= map->y_ticks) {addr.y = map->y_ticks - edgetics;}
    etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[2]));
    addr.x = (loc_x + 1) * edgetics;addr.y = (loc_y + 1) * edgetics;
    if (addr.x >= map->x_ticks) {addr.x = map->x_ticks - edgetics;}
    if (addr.y >= map->y_ticks) {addr.y = map->y_ticks - edgetics;}
    etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[3]));

    percent = fmod(rotated_point_x / map->spacing, map->spacing) / map->spacing;
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
int
cca_get_vs30_based_gtl(cca_point_t *point,
                       cca_properties_t *data) {
    cca_point_t vm_point;
    cca_properties_t vm_data;

    double a = 0.5, b = 0.6, c = 0.5;
    double percent_z = point->depth / cca_configuration->depth_interval;
    double f = 0.0, g = 0.0;
    double vs30 = 0.0, vp30 = 0.0;

    // Double check that we're above the first layer.
    if (percent_z > 1) {return (CCA_CODE_ERROR);}

    vm_point.latitude = point->latitude;
    vm_point.longitude = point->longitude;
    vm_point.depth = cca_configuration->depth_interval;
    memset(&vm_data, 0, sizeof(vm_data));
    if (cca_query(&vm_point, &vm_data, 1, NULL) != CCA_CODE_SUCCESS) {return (CCA_CODE_ERROR);}

    // Now we need the Vs30 data value.
    vs30 = cca_get_vs30_value(point->longitude, point->latitude, cca_vs30_map);

    if (vs30 == -1) {
        data->vp = -1;
        data->vs = -1;
    } else {
        // Get the point's material properties within the GTL.
        f = percent_z + b * (percent_z - pow(percent_z, 2.0f));
        g = a - a * percent_z + c * (pow(percent_z, 2.0f) + 2.0 * sqrt(percent_z) - 3.0 * percent_z);
        data->vs = f * vm_data.vs + g * vs30;

        vs30 = vs30 / 1000;
        vp30 = 0.9409 + 2.0947 * vs30 - 0.8206 * pow(vs30, 2.0f) + 0.2683 * pow(vs30, 3.0f) - 0.0251 * pow(vs30, 4.0f);
        vp30 = vp30 * 1000;
        data->vp = f * vm_data.vp + g * vp30;
    }

    return CCA_CODE_SUCCESS;
}


// The following functions are for dynamic library mode. If we are compiling
// a static library, these functions must be disabled to avoid conflicts.
#if defined(BUILD_SHARED_LIBRARY)

/**
 * Create function loaded and called by the UCVM library.
 *
 * @param dir The directory in which UCVM is installed.
 * @return CCA_CODE_SUCCESS or CCA_CODE_ERRORure.
 */
int
ucvmapi_model_create(const char *dir,
                     const char *label) {
    return cca_create(dir, label);
}


/**
 * Initialize function loaded and called by the UCVM library.
 *
 * @return CCA_CODE_SUCCESS or CCA_CODE_ERRORure.
 */
int
ucvmapi_model_initialize(void) {
    return cca_initialize();
}


/**
 * Finalize function loaded and called by the UCVM library. Calls cca_finalize.
 *
 * @return CCA_CODE_SUCCESS
 */
int
ucvmapi_model_finalize(void) {
    return cca_finalize();
}


/**
 * Version function loaded and called by the UCVM library.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int
ucvmapi_model_version(char *ver,
                      int len) {
    return cca_version(ver, len);
}


/**
 * Set model parameter.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int
ucvmapi_model_set_parameter(const char *name,
                            const char *value) {
    return cca_set_parameter(name, value);
}


/**
 * Query function loaded and called by the UCVM library. Calls cca_query.
 *
 * @param points The basic_point_t array containing the points.
 * @param data The basic_properties_t array containing the material properties returned.
 * @param numpoints The number of points in the array.
 * @return CCA_CODE_SUCCESS or CCA_CODE_ERROR.
 */
int
ucvmapi_model_query(cca_point_t *points,
                    cca_properties_t *data,
                    int numpoints,
                    cca_query_flags_t *qflags) {
    return cca_query(points, data, numpoints, qflags);
}


#endif
