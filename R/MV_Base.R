# MV_Base.R
#
# This file previously contained make_spd, mv_mixture_fit, and mv_noise_ks_score.
# Those functions have been extracted to non-legacy core modules:
#
#   - make_spd          -> R/core_spd.R
#   - mv_mixture_fit    -> R/mv_fit.R
#   - mv_noise_ks_score -> R/mv_fit.R
#
# All MV and RIMLE callers now resolve from non-legacy core files only.
# This file is retained as a historical placeholder.
