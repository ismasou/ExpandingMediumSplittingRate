#!/bin/env bash

BOLDRED="\033[1;31m"
RESET="\033[0m"
# This script generates the data to compare with Caron-Huot Gale's paper

# Compile the program
if ! make Expanding; then
    echo -e "${BOLDRED}Compilation failed${RESET}"
    exit 1
fi

# mkdir
mkdir -p Output
mkdir -p Output/Full
mkdir -p Output/Expanding
mkdir -p Output/Opacity

# Run the program
echo -e "${BOLDRED}Running first set of simulations${RESET}"
./Expanding.out -P 16 -T 0.4 -z 0.01 -a 1.0 -t0 0.1 -tmax 15 -Process G_GG
echo -e "${BOLDRED}Running second set of simulations${RESET}"
./Expanding.out -P 16 -T 0.4 -z 0.01 -a 1.0 -t0 0.3 -tmax 15 -Process G_GG
echo -e "${BOLDRED}Running third set of simulations${RESET}"
./Expanding.out -P 16 -T 0.4 -z 0.5  -a 1.0 -t0 0.1 -tmax 15 -Process G_GG
echo -e "${BOLDRED}Running fourth set of simulations${RESET}"
./Expanding.out -P 16 -T 0.4 -z 0.5  -a 1.0 -t0 0.3 -tmax 15 -Process G_GG

# Copy results to CH_G_Comparison
echo -e "${BOLDRED}Copying results to Plot Folder${RESET}"
cp -f Output/Expanding/rate_P16_z0.01_T0.4_alpha1_t00.1_t10.dat ExpandingPlot/Expanding/rate_P16_z0.01_T0.4_alpha1_t00.1.dat
cp -f Output/Expanding/rate_P16_z0.01_T0.4_alpha1_t00.3_t10.dat ExpandingPlot/Expanding/rate_P16_z0.01_T0.4_alpha1_t00.3.dat
cp -f Output/Expanding/rate_P16_z0.5_T0.4_alpha1_t00.1_t10.dat  ExpandingPlot/Expanding/rate_P16_z0.5_T0.4_alpha1_t00.1.dat
cp -f Output/Expanding/rate_P16_z0.5_T0.4_alpha1_t00.3_t10.dat  ExpandingPlot/Expanding/rate_P16_z0.5_T0.4_alpha1_t00.3.dat

# Compile Opacity Expansion
if ! make ExpandingOp; then
    echo -e "${BOLDRED}Compilation failed${RESET}"
    exit 1
fi

# Run Opacity Expansion
echo -e "${BOLDRED}Running Opacity Expansion${RESET}"
./OEExp.out -P 16 -T 0.4 -z 0.01 -a 1.0 -t0 0.1 -tmax 15 -Process G_GG
./OEExp.out -P 16 -T 0.4 -z 0.01 -a 1.0 -t0 0.3 -tmax 15 -Process G_GG
./OEExp.out -P 16 -T 0.4 -z 0.5  -a 1.0 -t0 0.1 -tmax 15 -Process G_GG
./OEExp.out -P 16 -T 0.4 -z 0.5  -a 1.0 -t0 0.3 -tmax 15 -Process G_GG

# Copy results to CH_G_Comparison
echo -e "${BOLDRED}Copying OE results to Plot Folder${RESET}"
cp -f Output/Opacity/bjorken_P16_z0.01_T0.4_alpha1_t00.1_t10.dat ExpandingPlot/Opacity/bjorken_P16_z0.01_T0.4_alpha1_t00.1.dat
cp -f Output/Opacity/bjorken_P16_z0.01_T0.4_alpha1_t00.3_t10.dat ExpandingPlot/Opacity/bjorken_P16_z0.01_T0.4_alpha1_t00.3.dat
cp -f Output/Opacity/bjorken_P16_z0.5_T0.4_alpha1_t00.1_t10.dat  ExpandingPlot/Opacity/bjorken_P16_z0.5_T0.4_alpha1_t00.1.dat
cp -f Output/Opacity/bjorken_P16_z0.5_T0.4_alpha1_t00.3_t10.dat  ExpandingPlot/Opacity/bjorken_P16_z0.5_T0.4_alpha1_t00.3.dat

# Generate Plots
echo -e "${BOLDRED}Generating plots${RESET}"
cd CH_G_Comparison || (echo -e "Directory not found" && exit 1)
gnuplot plot.gp

cd .. || echo -e "${BOLDRED}Directory not found${RESET}" && exit 1
echo -e "Plots generated in ExpandingPlot/"
