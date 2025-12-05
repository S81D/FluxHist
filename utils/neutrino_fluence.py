import uproot
import numpy as np
import argparse

# Load Flux histogram, determine total neutrino flux per cm^2 (central values calculated from gsimple files)

'''
usage: neutrino_fluence.py [-h] [-i INPUT] {numu,nue,numubar,nuebar} ...

Compute integrated neutrino fluxes from the output ROOT file.

positional arguments:
  {numu,nue,numubar,nuebar}
                        Neutrino flavors to integrate.

options:
  -i, --input INPUT     Path to the output flux ROOT file.
'''


# flavor map to histogram name
FLAVOR_MAP = {
    "numu":    "numu_cv",
    "nue":     "nue_cv",
    "numubar": "numubar_cv",
    "nuebar":  "nuebar_cv"
}


# return integral of histogram
def integrate_flux(hist):
    values, edges = hist.to_numpy()
    return float(np.sum(values))


def main():
    parser = argparse.ArgumentParser(
        description="Compute integrated neutrino fluxes from output flux ROOT file."
    )

    parser.add_argument(
        "flavors",
        nargs="+",
        choices=FLAVOR_MAP.keys(),
        help="Neutrino flavors to integrate (numu, nue, numubar, nuebar)."
    )

    parser.add_argument(
        "-i", "--input",
        default='ANNIE_FLUX_full_detector.root`,
        help='Path to the output flux ROOT file.'
    )

    args = parser.parse_args()

    with uproot.open(args.input) as root_file:
        total_neutrino_flux = 0.0
        total_antineutrino_flux = 0.0

        print()
        for flavor in args.flavors:
            hist_name = FLAVOR_MAP[flavor]
            hist = root_file[hist_name]

            flux = integrate_flux(hist)
            print(f"Integrated {flavor} flux: {flux:.6g} per POT per cm^2")

            # group flavors into total integrated (anti)neutrino fluxes
            if "bar" in flavor:   # numubar or nuebar
                total_antineutrino_flux += flux
            else:
                total_neutrino_flux += flux

        print("\n--- Totals ---")
        if total_neutrino_flux > 0:
            print(f"Total neutrino flux: {total_neutrino_flux:.6g}")
        if total_antineutrino_flux > 0:
            print(f"Total antineutrino flux: {total_antineutrino_flux:.6g}")

    print("\ndone\n")


if __name__ == "__main__":
    main()

# done
