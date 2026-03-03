#!/usr/bin/env python3

chi2_min = 5
chi2_max = 75   # inclusive

pt_bins = [
    "0-1",
    "1-2",
    "2-5",
    "5-15",
]

y_range = "2.5-3.6"
cent_range = "0-100"
fitStat = "todo"
doIterativeFit = "true"

modelSig = "extCB"
mean_mass = "[3.096,3,3.2]"
sigma_mass = "[0.05,0.04,0.07]"
alpha0_mass = "[0.763]"
n0_mass = "[4.963]"
alpha1_mass = "[2.125]"
n1_mass = "[0.889]"

modelBkg = "ChevPol3"
c0_mass = "[-1,-2,2]"
c1_mass = "[0.5,-2,2]"
c2_mass = "[-0.1,-0.4,0.4]"
c3_mass = "[0.5,-0.1,0.1]"
# c4_mass = "[0,-0.1,0.1]"
# c5_mass = "[0,-0.1,0.1]"
# c6_mass = "[0,-0.1,0.1]"

output_file = "chi2_pO_scan_table.txt"

with open(output_file, "w") as f:
    # Header
    f.write(
        "pt;y;cent;chi2;fitStat;doIterativeFit;modelSig_mass;mean_mass;"
        "sigma_mass;alpha0_mass;n0_mass;alpha1_mass;n1_mass;"
        "modelBkg_mass;c0_mass;c1_mass;c2_mass;c3_mass;\n"
    )

    # Content
    for chi2 in range(chi2_min, chi2_max + 1):
        chi2_bin = f"0-{chi2}"
        for pt in pt_bins:
            line = (
                f"{pt};{y_range};{cent_range};{chi2_bin};{fitStat};{doIterativeFit};"
                f"{modelSig};{mean_mass};{sigma_mass};{alpha0_mass};"
                f"{n0_mass};{alpha1_mass};{n1_mass};"
                f"{modelBkg};{c0_mass};{c1_mass};{c2_mass};{c3_mass};\n"
            )
            f.write(line)

print(f"Written {4 * (chi2_max - chi2_min + 1)} rows to {output_file}")

