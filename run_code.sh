# Argument	Type	Description
# ispO	bool	Use p–O collision settings
# isMC	bool	Run on Monte Carlo sample
# caseName	string	Tag for input/output files, useful for systematics
# remakeDS	bool	Rebuild the RooDataSet from reduced tables
# fitMass1D	bool	Perform 1D mass fit
# fitTauz1D	bool	Perform 1D τz fit
# fit2D	bool	Perform 2D mass–τz fit
# plotResults	bool	Produce and save summary plots

root -l 'InputToResults.C(false, false, "fineBins", true, true, true, true, true)'
