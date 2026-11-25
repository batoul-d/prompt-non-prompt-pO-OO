#!/bin/bash

# ============================================================
# CONFIGURATION
# ============================================================

# TODO: add sweight plots in a different part of the html file

# Run analysis
RUN_ANALYSIS=true
DO_QC=true
SAVE_LOG=false
SPLOT_QA=true

# Collision systems
PO=true
OO=false

# Case names
# PO_CASE="fb_test"
PO_CASE="fineBins_copy"
PO_CASE_MC="fineBins_copy"
OO_CASE="fineBins_VWG_4Gaussians"
OO_CASE_MC="fineBins_VWG_4Gaussians"

# Signal extraction method
METHOD="sPlot"

# Fit options
REMAKE_DS=true
FIT_MASS_1D=true
FIT_TAUZ_1D=true
FIT_TAUZ_NP_1D=true
FIT_2D=false
PLOT_RESULTS=false

# ============================================================
# RUN ANALYSIS
# ============================================================

if [ "$RUN_ANALYSIS" = true ]; then

    mkdir -p plots

    LOG_FILE_PO="plots/log_pO.txt"
    LOG_FILE_OO="plots/log_OO.txt"

    if [ "$SAVE_LOG" = true ]; then
        > "$LOG_FILE_PO"
        > "$LOG_FILE_OO"
    fi

    if [ "$PO" = true ]; then
        echo "Running pO analysis: $PO_CASE"

        if [ "$SAVE_LOG" = true ]; then
            echo "Running pO analysis: $PO_CASE" >> "$LOG_FILE_PO"
            root -l -q "InputToResults.C(true, false, \"$PO_CASE\", \"$METHOD\", $REMAKE_DS, $FIT_MASS_1D, $FIT_TAUZ_1D, $FIT_TAUZ_NP_1D, $FIT_2D, $PLOT_RESULTS)" >> "$LOG_FILE_PO" 2>&1
        else
            root -l -q "InputToResults.C(true, false, \"$PO_CASE\", \"$METHOD\", $REMAKE_DS, $FIT_MASS_1D, $FIT_TAUZ_1D, $FIT_TAUZ_NP_1D, $FIT_2D, $PLOT_RESULTS)"
        fi
    fi

    if [ "$OO" = true ]; then
        echo "Running OO analysis: $OO_CASE"

        if [ "$SAVE_LOG" = true ]; then
            echo "Running OO analysis: $OO_CASE" >> "$LOG_FILE_OO"
            root -l -q "InputToResults.C(false, false, \"$OO_CASE\", \"$METHOD\", $REMAKE_DS, $FIT_MASS_1D, $FIT_TAUZ_1D, $FIT_TAUZ_NP_1D, $FIT_2D, $PLOT_RESULTS)" >> "$LOG_FILE_OO" 2>&1
        else
            root -l -q "InputToResults.C(false, false, \"$OO_CASE\", \"$METHOD\", $REMAKE_DS, $FIT_MASS_1D, $FIT_TAUZ_1D, $FIT_TAUZ_NP_1D, $FIT_2D, $PLOT_RESULTS)"
        fi
    fi

fi

# ============================================================
# COMPARISON PLOTS
# ============================================================

if [ "$DO_QC" = true ]; then

    root -l -q "compareDataAndMCfitParams.C(\"$PO_CASE\", \"$PO_CASE_MC\", \"$OO_CASE\", \"$OO_CASE_MC\")"
    python3 drawAllPlotsHtml.py "$PO_CASE" "$OO_CASE" "$SPLOT_QA"
    open plots/comparison_report.html

fi