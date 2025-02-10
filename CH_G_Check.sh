#!/bin/env bash

BOLDRED="\033[1;31m"
RESET="\033[0m"
# This script generates the data to compare with Caron-Huot Gale's paper

# Compile the program
if ! make Static; then
    echo -e "${BOLDRED}Compilation failed${RESET}"
    exit 1
fi

# mkdir
mkdir -p Output
mkdir -p Output/Full
mkdir -p Output/Evolution
mkdir -p Output/Opacity

# Run the program
echo -e "${BOLDRED}Running first set of simulations${RESET}"
./Static.out -P 16 -T 0.2 -z 0.1875 -Process Q_GQ
echo -e "${BOLDRED}Running second set of simulations${RESET}"
./Static.out -P 16 -T 0.2 -z 0.5 -Process Q_GQ
echo -e "${BOLDRED}Running third set of simulations${RESET}"
./Static.out -P 16 -T 0.4 -z 0.1875 -Process Q_GQ
echo -e "${BOLDRED}Running fourth set of simulations${RESET}"
./Static.out -P 16 -T 0.4 -z 0.5 -Process Q_GQ

# Copy results to CH_G_Comparison
echo -e "${BOLDRED}Copying results to CH_G_Comparison${RESET}"
cp -f Output/Full/rate_P16_z0.1875_T0.2.dat CH_G_Comparison/rate-0.2_3.dat
cp -f Output/Full/rate_P16_z0.1875_T0.4.dat CH_G_Comparison/rate-0.4_3.dat
cp -f Output/Full/rate_P16_z0.5_T0.2.dat CH_G_Comparison/rate-0.2_8.dat
cp -f Output/Full/rate_P16_z0.5_T0.4.dat CH_G_Comparison/rate-0.4_8.dat

# Compile Opacity Expansion
if ! make Opacity; then
    echo -e "${BOLDRED}Compilation failed${RESET}"
    exit 1
fi

# Run Opacity Expansion
echo -e "${BOLDRED}Running Opacity Expansion${RESET}"
./Opacity.out -P 16 -T 0.2 -z 0.1875 -Process Q_GQ
./Opacity.out -P 16 -T 0.2 -z 0.5 -Process Q_GQ
./Opacity.out -P 16 -T 0.4 -z 0.1875 -Process Q_GQ
./Opacity.out -P 16 -T 0.4 -z 0.5 -Process Q_GQ

# Copy results to CH_G_Comparison
echo -e "${BOLDRED}Copying OE results to CH_G_Comparison${RESET}"
cp -f Output/Opacity/opacity_P16_z0.1875_T0.2.dat CH_G_Comparison/rate-opacity-0.2_3.dat
cp -f Output/Opacity/opacity_P16_z0.1875_T0.4.dat CH_G_Comparison/rate-opacity-0.4_3.dat
cp -f Output/Opacity/opacity_P16_z0.5_T0.2.dat CH_G_Comparison/rate-opacity-0.2_8.dat
cp -f Output/Opacity/opacity_P16_z0.5_T0.4.dat CH_G_Comparison/rate-opacity-0.4_8.dat

# Generate Plots
echo -e "${BOLDRED}Generating plots${RESET}"
cd CH_G_Comparison || (echo -e "Directory not found" && exit 1)
gnuplot plot.gp

cd .. || echo -e "${BOLDRED}Directory not found${RESET}" && exit 1
echo -e "Plots generated in CH_G_Comparison"
