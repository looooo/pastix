# include(RulesPrecisions)
# # reset variables
# set(generated_files "")
# set(generated_headers "")

# set(HEADERS
#   z_spm.h
# )

# ### generate the dsparse_cores headers for all possible precisions
# precisions_rules_py(generated_headers
#                     "${HEADERS}"
#                     PRECISIONS "p;s;d;c;z")

# add_custom_target(spm_headers
#   DEPENDS ${generated_headers} )

