/** @file DynamicPropsCmd.hpp
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.23
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt */

#ifndef DYNAMICPROPSCMD_H
#define DYNAMICPROPSCMD_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "DynamicProps"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "DynamicProps"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION ""
#endif

enum enum_privilegedAxis { privilegedAxis__NULL = -1, privilegedAxis_arg_x = 0, privilegedAxis_arg_y, privilegedAxis_arg_z };
enum enum_selectionMode { selectionMode__NULL = -1, selectionMode_arg_survival = 0, selectionMode_arg_restart };

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * input_arg;	/**< @brief input dump file.  */
  char * input_orig;	/**< @brief input dump file original value given at command line.  */
  const char *input_help; /**< @brief input dump file help description.  */
  char * output_arg;	/**< @brief output file name.  */
  char * output_orig;	/**< @brief output file name original value given at command line.  */
  const char *output_help; /**< @brief output file name help description.  */
  char * sele1_arg;	/**< @brief select first stuntdouble set.  */
  char * sele1_orig;	/**< @brief select first stuntdouble set original value given at command line.  */
  const char *sele1_help; /**< @brief select first stuntdouble set help description.  */
  char * sele2_arg;	/**< @brief select second stuntdouble set (if sele2 is not set, use script from sele1).  */
  char * sele2_orig;	/**< @brief select second stuntdouble set (if sele2 is not set, use script from sele1) original value given at command line.  */
  const char *sele2_help; /**< @brief select second stuntdouble set (if sele2 is not set, use script from sele1) help description.  */
  char * sele3_arg;	/**< @brief select third stuntdouble set.  */
  char * sele3_orig;	/**< @brief select third stuntdouble set original value given at command line.  */
  const char *sele3_help; /**< @brief select third stuntdouble set help description.  */
  int seleoffset_arg;	/**< @brief global index offset for a second object (used to define a vector between sites in molecule).  */
  char * seleoffset_orig;	/**< @brief global index offset for a second object (used to define a vector between sites in molecule) original value given at command line.  */
  const char *seleoffset_help; /**< @brief global index offset for a second object (used to define a vector between sites in molecule) help description.  */
  int order_arg;	/**< @brief Lengendre Polynomial Order.  */
  char * order_orig;	/**< @brief Lengendre Polynomial Order original value given at command line.  */
  const char *order_help; /**< @brief Lengendre Polynomial Order help description.  */
  int nbins_arg;	/**< @brief Number of bins (default='100').  */
  char * nbins_orig;	/**< @brief Number of bins original value given at command line.  */
  const char *nbins_help; /**< @brief Number of bins help description.  */
  int nzbins_arg;	/**< @brief Number of Z bins (default='100').  */
  char * nzbins_orig;	/**< @brief Number of Z bins original value given at command line.  */
  const char *nzbins_help; /**< @brief Number of Z bins help description.  */
  double rcut_arg;	/**< @brief cutoff radius (angstroms).  */
  char * rcut_orig;	/**< @brief cutoff radius (angstroms) original value given at command line.  */
  const char *rcut_help; /**< @brief cutoff radius (angstroms) help description.  */
  double OOcut_arg;	/**< @brief Oxygen-Oxygen cutoff radius (angstroms) (default='3.5').  */
  char * OOcut_orig;	/**< @brief Oxygen-Oxygen cutoff radius (angstroms) original value given at command line.  */
  const char *OOcut_help; /**< @brief Oxygen-Oxygen cutoff radius (angstroms) help description.  */
  double thetacut_arg;	/**< @brief HOO cutoff angle (degrees) (default='30').  */
  char * thetacut_orig;	/**< @brief HOO cutoff angle (degrees) original value given at command line.  */
  const char *thetacut_help; /**< @brief HOO cutoff angle (degrees) help description.  */
  double OHcut_arg;	/**< @brief Oxygen-Hydrogen cutoff radius (angstroms) (default='2.45').  */
  char * OHcut_orig;	/**< @brief Oxygen-Hydrogen cutoff radius (angstroms) original value given at command line.  */
  const char *OHcut_help; /**< @brief Oxygen-Hydrogen cutoff radius (angstroms) help description.  */
  enum enum_privilegedAxis privilegedAxis_arg;	/**< @brief which axis is special for spatial analysis (default = z axis) (default='z').  */
  char * privilegedAxis_orig;	/**< @brief which axis is special for spatial analysis (default = z axis) original value given at command line.  */
  const char *privilegedAxis_help; /**< @brief which axis is special for spatial analysis (default = z axis) help description.  */
  double length_arg;	/**< @brief maximum length (default='100').  */
  char * length_orig;	/**< @brief maximum length original value given at command line.  */
  const char *length_help; /**< @brief maximum length help description.  */
  double dipoleX_arg;	/**< @brief X-component of the dipole with respect to body frame (default='0.0').  */
  char * dipoleX_orig;	/**< @brief X-component of the dipole with respect to body frame original value given at command line.  */
  const char *dipoleX_help; /**< @brief X-component of the dipole with respect to body frame help description.  */
  double dipoleY_arg;	/**< @brief Y-component of the dipole with respect to body frame (default='0.0').  */
  char * dipoleY_orig;	/**< @brief Y-component of the dipole with respect to body frame original value given at command line.  */
  const char *dipoleY_help; /**< @brief Y-component of the dipole with respect to body frame help description.  */
  double dipoleZ_arg;	/**< @brief Z-component of the dipole with respect to body frame (default='-1.0').  */
  char * dipoleZ_orig;	/**< @brief Z-component of the dipole with respect to body frame original value given at command line.  */
  const char *dipoleZ_help; /**< @brief Z-component of the dipole with respect to body frame help description.  */
  enum enum_selectionMode selectionMode_arg;	/**< @brief How to treat objects which leave a dynamic selection and then return later (default = survival) (default='survival').  */
  char * selectionMode_orig;	/**< @brief How to treat objects which leave a dynamic selection and then return later (default = survival) original value given at command line.  */
  const char *selectionMode_help; /**< @brief How to treat objects which leave a dynamic selection and then return later (default = survival) help description.  */
  const char *selecorr_help; /**< @brief selection correlation function help description.  */
  const char *rcorr_help; /**< @brief mean squared displacement help description.  */
  const char *rcorrZ_help; /**< @brief mean squared displacement binned by Z help description.  */
  const char *vcorr_help; /**< @brief velocity correlation function help description.  */
  const char *vcorrZ_help; /**< @brief velocity correlation function along z-axis help description.  */
  const char *vcorrR_help; /**< @brief velocity correlation function projected radially help description.  */
  const char *vaOutProdcorr_help; /**< @brief Velocity - Velocity auto outer product correlation function help description.  */
  const char *waOutProdcorr_help; /**< @brief Angular Velocity - Angular Velocity auto outer product correlation function help description.  */
  const char *vwOutProdcorr_help; /**< @brief Velocity - Angular Velocity outer product correlation function help description.  */
  const char *wvOutProdcorr_help; /**< @brief Angular Velocity - Velocity outer product correlation function help description.  */
  const char *wcorr_help; /**< @brief charge velocity correlation function help description.  */
  const char *dcorr_help; /**< @brief dipole correlation function help description.  */
  const char *lcorr_help; /**< @brief Lengendre correlation function help description.  */
  const char *lcorrZ_help; /**< @brief Lengendre correlation function binned by Z help description.  */
  const char *cohZ_help; /**< @brief Lengendre correlation function for OH bond vectors binned by Z help description.  */
  const char *sdcorr_help; /**< @brief System dipole correlation function help description.  */
  const char *r_rcorr_help; /**< @brief Radial msd help description.  */
  const char *thetacorr_help; /**< @brief Angular msd help description.  */
  const char *drcorr_help; /**< @brief Directional msd for particles with unit vectors help description.  */
  const char *stresscorr_help; /**< @brief Stress tensor correlation function help description.  */
  const char *bondcorr_help; /**< @brief Bond extension correlation function help description.  */
  const char *freqfluccorr_help; /**< @brief Frequency Fluctuation correlation function help description.  */
  const char *jumptime_help; /**< @brief Hydrogen bond jump time correlation function help description.  */
  const char *jumptimeZ_help; /**< @brief Hydrogen bond jump time correlation function binned by Z help description.  */
  const char *jumptimeR_help; /**< @brief Hydrogen bond jump time correlation function binned by R around a third selection help description.  */
  const char *persistence_help; /**< @brief Hydrogen bond persistence correlation function help description.  */
  const char *pjcorr_help; /**< @brief Momentum - Angular Momentum cross correlation function help description.  */
  const char *ftcorr_help; /**< @brief Force - Torque cross correlation function help description.  */
  const char *ckcorr_help; /**< @brief Charge - Kinetic energy cross correlation function help description.  */
  const char *cscorr_help; /**< @brief Charge - Orientation order parameter (Cos\theta) cross correlation function help description.  */
  const char *facorr_help; /**< @brief Force - Force auto correlation function help description.  */
  const char *tfcorr_help; /**< @brief Torque - Force Cross correlation function help description.  */
  const char *tacorr_help; /**< @brief Torque auto correlation function help description.  */
  const char *disp_help; /**< @brief Displacement correlation function help description.  */
  const char *dispZ_help; /**< @brief Displacement correlation function binned by Z help description.  */
  const char *current_help; /**< @brief Current density auto correlation function help description.  */
  const char *onsager_help; /**< @brief Onsager coefficient correlation functions help description.  */
  const char *ddisp_help; /**< @brief Collective Dipole displacement function (Helfand moment of Current Density) help description.  */
  const char *rotAngleDisp_help; /**< @brief Displacement correlation function for rotation angles help description.  */
  const char *meandisp_help; /**< @brief mean displacement help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int input_given ;	/**< @brief Whether input was given.  */
  unsigned int output_given ;	/**< @brief Whether output was given.  */
  unsigned int sele1_given ;	/**< @brief Whether sele1 was given.  */
  unsigned int sele2_given ;	/**< @brief Whether sele2 was given.  */
  unsigned int sele3_given ;	/**< @brief Whether sele3 was given.  */
  unsigned int seleoffset_given ;	/**< @brief Whether seleoffset was given.  */
  unsigned int order_given ;	/**< @brief Whether order was given.  */
  unsigned int nbins_given ;	/**< @brief Whether nbins was given.  */
  unsigned int nzbins_given ;	/**< @brief Whether nzbins was given.  */
  unsigned int rcut_given ;	/**< @brief Whether rcut was given.  */
  unsigned int OOcut_given ;	/**< @brief Whether OOcut was given.  */
  unsigned int thetacut_given ;	/**< @brief Whether thetacut was given.  */
  unsigned int OHcut_given ;	/**< @brief Whether OHcut was given.  */
  unsigned int privilegedAxis_given ;	/**< @brief Whether privilegedAxis was given.  */
  unsigned int length_given ;	/**< @brief Whether length was given.  */
  unsigned int dipoleX_given ;	/**< @brief Whether dipoleX was given.  */
  unsigned int dipoleY_given ;	/**< @brief Whether dipoleY was given.  */
  unsigned int dipoleZ_given ;	/**< @brief Whether dipoleZ was given.  */
  unsigned int selectionMode_given ;	/**< @brief Whether selectionMode was given.  */
  unsigned int selecorr_given ;	/**< @brief Whether selecorr was given.  */
  unsigned int rcorr_given ;	/**< @brief Whether rcorr was given.  */
  unsigned int rcorrZ_given ;	/**< @brief Whether rcorrZ was given.  */
  unsigned int vcorr_given ;	/**< @brief Whether vcorr was given.  */
  unsigned int vcorrZ_given ;	/**< @brief Whether vcorrZ was given.  */
  unsigned int vcorrR_given ;	/**< @brief Whether vcorrR was given.  */
  unsigned int vaOutProdcorr_given ;	/**< @brief Whether vaOutProdcorr was given.  */
  unsigned int waOutProdcorr_given ;	/**< @brief Whether waOutProdcorr was given.  */
  unsigned int vwOutProdcorr_given ;	/**< @brief Whether vwOutProdcorr was given.  */
  unsigned int wvOutProdcorr_given ;	/**< @brief Whether wvOutProdcorr was given.  */
  unsigned int wcorr_given ;	/**< @brief Whether wcorr was given.  */
  unsigned int dcorr_given ;	/**< @brief Whether dcorr was given.  */
  unsigned int lcorr_given ;	/**< @brief Whether lcorr was given.  */
  unsigned int lcorrZ_given ;	/**< @brief Whether lcorrZ was given.  */
  unsigned int cohZ_given ;	/**< @brief Whether cohZ was given.  */
  unsigned int sdcorr_given ;	/**< @brief Whether sdcorr was given.  */
  unsigned int r_rcorr_given ;	/**< @brief Whether r_rcorr was given.  */
  unsigned int thetacorr_given ;	/**< @brief Whether thetacorr was given.  */
  unsigned int drcorr_given ;	/**< @brief Whether drcorr was given.  */
  unsigned int stresscorr_given ;	/**< @brief Whether stresscorr was given.  */
  unsigned int bondcorr_given ;	/**< @brief Whether bondcorr was given.  */
  unsigned int freqfluccorr_given ;	/**< @brief Whether freqfluccorr was given.  */
  unsigned int jumptime_given ;	/**< @brief Whether jumptime was given.  */
  unsigned int jumptimeZ_given ;	/**< @brief Whether jumptimeZ was given.  */
  unsigned int jumptimeR_given ;	/**< @brief Whether jumptimeR was given.  */
  unsigned int persistence_given ;	/**< @brief Whether persistence was given.  */
  unsigned int pjcorr_given ;	/**< @brief Whether pjcorr was given.  */
  unsigned int ftcorr_given ;	/**< @brief Whether ftcorr was given.  */
  unsigned int ckcorr_given ;	/**< @brief Whether ckcorr was given.  */
  unsigned int cscorr_given ;	/**< @brief Whether cscorr was given.  */
  unsigned int facorr_given ;	/**< @brief Whether facorr was given.  */
  unsigned int tfcorr_given ;	/**< @brief Whether tfcorr was given.  */
  unsigned int tacorr_given ;	/**< @brief Whether tacorr was given.  */
  unsigned int disp_given ;	/**< @brief Whether disp was given.  */
  unsigned int dispZ_given ;	/**< @brief Whether dispZ was given.  */
  unsigned int current_given ;	/**< @brief Whether current was given.  */
  unsigned int onsager_given ;	/**< @brief Whether onsager was given.  */
  unsigned int ddisp_given ;	/**< @brief Whether ddisp was given.  */
  unsigned int rotAngleDisp_given ;	/**< @brief Whether rotAngleDisp was given.  */
  unsigned int meandisp_given ;	/**< @brief Whether meandisp was given.  */

  char **inputs ; /**< @brief unnamed options (options without names) */
  unsigned inputs_num ; /**< @brief unnamed options number */
  int correlation_function_group_counter; /**< @brief Counter for group correlation_function */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);

extern const char *cmdline_parser_privilegedAxis_values[];  /**< @brief Possible values for privilegedAxis. */
extern const char *cmdline_parser_selectionMode_values[];  /**< @brief Possible values for selectionMode. */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* DYNAMICPROPSCMD_H */
