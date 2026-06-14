# ====================================================================================== #
# Executes the entire analysis pipeline
#
# ⚠️ NOTE
# - PMC runs will take multiple hours to execute
# - Selection of PMoAs for cross-validation and validation is done manually in scripts
# ====================================================================================== #

# ---- 2,4-D

include("calibration_24D.jl") # fit 2,4-D model
include("calibration_24D_summaries.jl") # additional summary output for 2,4-D calibrations
include("cross_validation_24D.jl") # cross-validation for 2,4-D

# ---- Flupyradifurone

include("calibration_flupyradifurone.jl") # fit Flupyradifurone model
include("calibration_flupyradifurone_summaries.jl") # additional summary output for Flup. calibs
include("cross_validation_flupyradifurone.jl") # cross-validation for Flupyradifurone

# ---- mixture

include("validation_mixure.jl") 

