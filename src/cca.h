/**
 * @file cca.h
 * @brief Main header file for CCA library.
 *
 * Delivers the CCA model which consists of En-Jui Lee's full 3D
 * tomographic results for central California.
 *
 */

#include "euclid/etree.h"

/** Defines a return value of success */
#define CCA_CODE_SUCCESS 0
/** Defines a return value for an error */
#define CCA_CODE_ERROR 1

typedef enum { CCA_COORD_GEO_DEPTH=0,
               CCA_COORD_GEO_ELEV } cca_ctype_t;

// Structures

/** Defines a point (latitude, longitude, and depth) in WGS84 format */
typedef struct cca_point_t {
    double longitude; /* longitude, degrees */
    double latitude; /* latitude, degrees */
    double depth; /* depth, meters */
} cca_point_t;

/** Defines the material properties this model will retrieve. */
typedef struct cca_properties_t {
    double vp; /* P wave velocity, meters per second */
    double vs; /* S wave velocity, meters per second */
    double rho; /* density, g/m**3 */
    double qp;
    double qs;
} cca_properties_t;

/** The CVM-S5 configuration structure. */
typedef struct cca_configuration_t {
    /** The zone of UTM projection */
    int utm_zone;
    /** The model directory */
    char model_dir[128];
    /** GTL on or off (1 or 0) */
    int gtl;
    /** Number of x points */
    int nx;
    /** Number of y points */
    int ny;
    /** Number of z points */
    int nz;
    /** Depth in meters */
    double depth;
    /** Top left corner easting in UTM projection */
    double top_left_corner_e;
    /** Top left corner northing in UTM projection */
    double top_left_corner_n;
    /** Top right corner easting in UTM projection */
    double top_right_corner_e;
    /** Top right corner northing in UTM projection */
    double top_right_corner_n;
    /** Bottom left corner easting in UTM projection */
    double bottom_left_corner_e;
    /** Bottom left corner northing in UTM projection */
    double bottom_left_corner_n;
    /** Bottom right corner easting in UTM projection */
    double bottom_right_corner_e;
    /** Bottom right corner northing in UTM projection */
    double bottom_right_corner_n;
    /** Z interval for the data */
    double depth_interval;
    /** Brocher 2005 scaling polynomial coefficient 10^0 */
    double p0;
    /** Brocher 2005 scaling polynomial coefficient 10^1 */
    double p1;
    /** Brocher 2005 scaling polynomial coefficient 10^2 */
    double p2;
    /** Brocher 2005 scaling polynomial coefficient 10^3 */
    double p3;
    /** Brocher 2005 scaling polynomial coefficient 10^4 */
    double p4;
    /** Brocher 2005 scaling polynomial coefficient 10^5 */
    double p5;
} cca_configuration_t;

/** The configuration structure for the Vs30 map. */
typedef struct cca_vs30_map_config_t {
    /** Pointer to the e-tree file */
    etree_t *vs30_map;
    /** The type of map */
    char type[20];
    /** A description of the map */
    char description[50];
    /** The map's author */
    char author[30];
    /** The date the map was created */
    char date[10];
    /** The spacing in meters */
    double spacing;
    /** The map's schema */
    char schema[50];
    /** The projection string in Proj.4 format */
    char projection[128];
    /** The origin point */
    cca_point_t origin_point;
    /** The number of degrees the map was rotated around origin */
    double rotation;
    /** The X dimension of the map */
    double x_dimension;
    /** The Y dimension of the map */
    double y_dimension;
    /** The Z dimension of the map */
    double z_dimension;
    /** Number of e-tree ticks in the X direction */
    int x_ticks;
    /** Number of e-tree ticks in the Y direction */
    int y_ticks;
    /** Number of e-tree ticks in the Z direction */
    int z_ticks;
} cca_vs30_map_config_t;

/** The model structure which points to available portions of the model. */
typedef struct cca_model_t {
    /** A pointer to the Vs data either in memory or disk. Null if does not exist. */
    void *vs;
    /** Vs status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
    int vs_status;
    /** A pointer to the Vp data either in memory or disk. Null if does not exist. */
    void *vp;
    /** Vp status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
    int vp_status;
    /** A pointer to the rho data either in memory or disk. Null if does not exist. */
    void *rho;
    /** Rho status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
    int rho_status;
    /** A pointer to the Qp data either in memory or disk. Null if does not exist. */
    void *qp;
    /** Qp status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
    int qp_status;
    /** A pointer to the Qs data either in memory or disk. Null if does not exist. */
    void *qs;
    /** Qs status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
    int qs_status;
} cca_model_t;

/** Contains the Vs30 and surface values from the UCVM map. */
typedef struct cca_vs30_mpayload_t {
    float surf; /* Surface elevation, meters */
    float vs30; /* Vs30 from Wills and Wald, units? */
} cca_vs30_mpayload_t;

// Constants

// UCVM API Required Functions

#if defined(BUILD_SHARED_LIBRARY)

/** Creates the model */
int ucvmapi_model_create(const char *models_dir,
                         const char *label);

/** Initializes the model */
int ucvmapi_model_initialize();

/** Cleans up the model (frees memory, etc.) */
int ucvmapi_model_finalize();

/** Returns version information */
int ucvmapi_model_version(char *ver,
                          int len);

/* Set model user parameter */
int ucvmapi_model_set_param(const char* name,
                            const char* value);

/** Queries the model */
int ucvmapi_model_query(cca_point_t *points,
                        cca_properties_t *data,
                        int numpts);

#endif

// CCA Related Functions

/** Creates the model */
int cca_create(const char *models_dir,
               const char *label);

/** Initializes the model */
int cca_initialize();

/** Cleans up the model (frees memory, etc.) */
int cca_finalize();

/** Returns version information */
int cca_version(char *ver,
                int len);

/* Set model user parameter */
int cca_set_param(const char* name,
                  const char* value);

/** Queries the model */
int cca_query(cca_point_t *points,
              cca_properties_t *data,
              int numpts);

// Non-UCVM Helper Functions
/** Reads the configuration file. */
int cca_read_configuration(char *file,
                           cca_configuration_t *config);

/** Prints out the error string. */
void cca_print_error(char *err);

/** Retrieves the value at a specified grid point in the model. */
void cca_read_properties(int x,
                         int y,
                         int z,
                         cca_properties_t *data);

/** Attempts to malloc the model size in memory and read it in. */
int cca_try_reading_model(cca_model_t *model);

/** Calculates density from Vs. */
double cca_calculate_density(double vs);

// GTL related
/** Retrieves the vs30 value for a given point. */
int cca_get_vs30_based_gtl(cca_point_t *point,
                           cca_properties_t *data);

/** Reads the specified Vs30 map from UCVM. */
int cca_read_vs30_map(char *filename,
                      cca_vs30_map_config_t *map);

/** Gets the Vs30 value at a point */
double cca_get_vs30_value(double longitude,
                          double latitude,
                          cca_vs30_map_config_t *map);

// Interpolation Functions
/** Linearly interpolates two cca_properties_t structures */
void cca_linear_interpolation(double percent,
                              cca_properties_t *x0,
                              cca_properties_t *x1,
                              cca_properties_t *ret_properties);

/** Bilinearly interpolates the properties. */
void cca_bilinear_interpolation(double x_percent,
                                double y_percent,
                                cca_properties_t *four_points,
                                cca_properties_t *ret_properties);

/** Trilinearly interpolates the properties. */
void cca_trilinear_interpolation(double x_percent,
                                 double y_percent,
                                 double z_percent,
                                 cca_properties_t *eight_points,
                                 cca_properties_t *ret_properties);
