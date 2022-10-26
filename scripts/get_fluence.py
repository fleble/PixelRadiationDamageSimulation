import ROOT as rt
import argparse

class GetFlunce(object):
    def __init__(self, flunce_root_file = "FLUKA/fluence_field.root"):
        self.__fl_root_file__ = rt.TFile.Open(flunce_root_file)
        self.__fl_histo__ = self.__fl_root_file__.Get("fluence_allpart_6500GeV_phase1")

    def fluence(self, r, z):
        bin_x = self.__fl_histo__.GetXaxis().FindBin(r)
        bin_y = self.__fl_histo__.GetYaxis().FindBin(z)
        fl = self.__fl_histo__.GetBinContent(bin_x, bin_y)
        return fl

def main(args):
    gf = GetFlunce()
    print gf.fluence(float(args.r), float(args.z))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', dest='r', action='store', help='radius')
    parser.add_argument('-z', dest='z', action='store', help='z')
    args = parser.parse_args()
    main(args)

# 5698
