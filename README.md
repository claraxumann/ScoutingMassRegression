# ScoutingMassRegression

fulltraining and fulltraining20epochs contain histograms, mass resolution, and efficiency plots for various trainings that I've run. Each training should have more or less the same plots. 

lossfunctions has plots comparing different loss functions. 
 - Loss 1: RatioSmoothL1Loss
 - Loss 2: SymmetricRatioSmoothL1Loss
 - Loss 3: LogCoshLoss (default)
 - Loss 4: SmoothL1Loss

The input files are pred.root files stored in various folders in `/eos/user/c/clxu/training_output`.
