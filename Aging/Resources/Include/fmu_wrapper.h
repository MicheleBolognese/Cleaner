/*
 Copyright 2021 Modelon AB
*/
#ifndef FMU_WRAPPER_HEADER
#define FMU_WRAPPER_HEADER
#include <stdio.h>

#if !defined(Ext_Obj_Export)
    #if defined _WIN32 || defined __CYGWIN__
        #define Ext_Obj_Export __declspec(dllexport)
    #else
        #if defined __GNUC__
            #define Ext_Obj_Export __attribute__ ((visibility("default")))
        #else
            #define Ext_Obj_Export
        #endif
    #endif
#endif


typedef void (ModelicaErrorTYPE)(const char *msg);

Ext_Obj_Export void *initialize_wrapper(const char *instance_name,
                                        const char *identifier,
                                        const char *ext_lib_base_path,
                                        const char *GUID,
                                        ModelicaErrorTYPE *ModelicaError);
Ext_Obj_Export void wrapper_get_real(void *fmu_wr, const unsigned int *value_refs,
                                     size_t n_values, double *output_values,
                                     ModelicaErrorTYPE *ModelicaError);

Ext_Obj_Export void wrapper_set_real(void *fmu_wr, const unsigned int *value_refs,
                                     size_t n_values, double *input_values,
                                     ModelicaErrorTYPE *ModelicaError);

Ext_Obj_Export void wrapper_get_integer(void *fmu_wr, const unsigned int *value_refs,
                                        size_t n_values, int *output_values,
                                        ModelicaErrorTYPE *ModelicaError);

Ext_Obj_Export void wrapper_set_integer(void *fmu_wr, const unsigned int *value_refs,
                                        size_t n_values, int *input_values,
                                        ModelicaErrorTYPE *ModelicaError);

Ext_Obj_Export void wrapper_get_boolean(void *fmu_wr, const unsigned int *value_refs,
                                        size_t n_values, int *output_values,
                                        ModelicaErrorTYPE *ModelicaError);

Ext_Obj_Export void wrapper_set_boolean(void *fmu_wr, const unsigned int *value_refs,
                                        size_t n_values, int *input_values,
                                        ModelicaErrorTYPE *ModelicaError);

Ext_Obj_Export void wrapper_get_string(void *fmu_wr, const unsigned int *value_refs,
                                        size_t n_values, char **output_values,
                                        ModelicaErrorTYPE *ModelicaError);

Ext_Obj_Export void wrapper_set_string(void *fmu_wr, const unsigned int *value_refs,
                                        size_t n_values, char **input_values,
                                        ModelicaErrorTYPE *ModelicaError);

Ext_Obj_Export void wrapper_setup_experiment(
    void *fmu_wr,
    int tolerance_defined,
    double tolerance,
    double start_time,
    int stop_time_defined,
    double stop_time,
    ModelicaErrorTYPE *ModelicaError
);
Ext_Obj_Export void wrapper_enter_initialization_mode(
    void *fmu_wr,
    ModelicaErrorTYPE *ModelicaError
);
Ext_Obj_Export void wrapper_exit_initialization_mode(
    void *fmu_wr,
    ModelicaErrorTYPE *ModelicaError
);
Ext_Obj_Export void wrapper_do_step(
    void *obj,
    double current_time,
    double step_size,
    int no_prior_fmu_state,
    ModelicaErrorTYPE *ModelicaError
);
Ext_Obj_Export void wrapper_set_real_input_der(
    void *fmu_wr,
    const int *value_refs,
    size_t n_values,
    int *order,
    double *input_values,
    ModelicaErrorTYPE *ModelicaError
);
Ext_Obj_Export void free_wrapper(void *fmu_wrapper);

#endif
