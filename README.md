# ScoutingMassRegression

fulltraining and fulltraining20epochs contain histograms, mass resolution, and efficiency plots for various trainings that I've run. Each training should have more or less the same plots. 

lossfunctions has plots comparing different loss functions. 
 - Loss 1: RatioSmoothL1Loss
 - Loss 2: SymmetricRatioSmoothL1Loss
 - Loss 3: LogCoshLoss (default)
 - Loss 4: SmoothL1Loss

The input files are pred.root files stored in various folders in `/eos/user/c/clxu/training_output`.

When running, make sure `events` is using the correct file for the training you want to evaluate. For the mass resolution and efficiency plots, I'd recommend running those cells in order; some of the variables have the same names and otherwise they might get mixed together.

datacorrections contains cuts and histograms to make corrections to the simulation data. It can take in a list of root files and separates them into TTToSemiLeptonic, data, and other. 
 - To be put in the right category, the file name must contain the correct substring
    - a file with "TTToSemiLeptonic" contained in the filename will have cuts made on it and be split into different histograms
    - a file with "data" in the name will be put into a single histogram, and similarly for "other"
 - The last cell of the notebook creates the histograms and currently can plot mass and pT, but can be changed to plot other variables as well.
