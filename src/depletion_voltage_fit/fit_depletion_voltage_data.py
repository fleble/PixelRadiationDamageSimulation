import argparse
from pathlib import Path
import datetime as dt

from array import array as pyroot_array
import ROOT

from utils.pythonUtils import import_module_from


ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)


def __get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input_simulation_file_name",
        help="Path to ROOT file with leakage current and depletion voltage simulation",
        required=True,
    )
    parser.add_argument(
        "-dv", "--depletion_voltage_data_file_name",
        help="Path to file with depletion voltage measurement data",
        required=True,
    )
    parser.add_argument(
        "-rog", "--readout_group",
        help="Read out group name, e.g. BPix_BmI_SEC4_LYR1",
    )
    # parser.add_argument(
    #     "-f", "--fluence_file",
    #     help='Input fluence ROOT file',
    #     default="data/fluence/fluence_field_phase1_6500GeV.root",
    # )
    parser.add_argument(
        "-od", "--output_directory",
        help="Output directory name",
        required=False,
        default="data/simulation/",
    )
    parser.add_argument(
        "-o", "--output_file_name",
        help="Output file name",
        required=True,
    )
    parser.add_argument(
        "-xmin", "--x_min",
        help="x-axis minimum",
        type=float,
    )
    parser.add_argument(
        "-xmax", "--x_max",
        help="x-axis maximum",
        type=float,
    )
    parser.add_argument(
        "-ymin", "--y_min",
        help="y-axis minimum",
        type=float,
    )
    parser.add_argument(
        "-ymax", "--y_max",
        help="y-axis maximum",
        type=float,
    )

    return parser.parse_args()


def __check_arguments(args):
    if ((args.xmin is not None and args.xmax is None)
        or (args.xmin is None and args.xmax is not None)
    ):
        print("Error: xmin and xmax must be both None or noth not None!")
        sys.exit(1)


def time_string_to_timestamp(time_string, time_format="%Y-%m-%d %H:%M:%S"):
    """Convert time in an arbitrary format to timestamp.

    Args:
        time_string (str): String representing time in an arbitrary format
        time_format (str): The format in which the time is provided

    Returns:
        float
    """

    return dt.datetime.strptime(time_string, time_format).timestamp()


def __read_data_from_simulation(simulation_file_name):
    simulation_file = ROOT.TFile.Open(simulation_file_name)

    keys_to_read = [
        "U_dep",
        "N_benef_anneal_g1",  #get_gA
        "N_revers_anneal_g1", #get_gY
        "N_nadefects_g1",     #get_gY
        "N_constdamage_g1",   #get_gC
        "fluence",
        "N_donor",
        "N_donor_stable_donorremoval",
    ]

    simulation_data = {
        k: simulation_file.Get(k) for k in keys_to_read
    }

    return simulation_data


def __read_depletion_voltage_data(depletion_voltage_data_file_name):
    return import_module_from(depletion_voltage_data_file_name).depletion_voltage_data


def __get_depletion_voltage_graph(depletion_voltage_data_file_name):
    depletion_voltage_data = __read_depletion_voltage_data(depletion_voltage_data_file_name)

    time_transformation = lambda x: time_string_to_timestamp(x, time_format="%Y-%m-%d %H:%M:%S")
    times = list(map(time_transformation, depletion_voltage_data.keys()))
    values = list(map(lambda x: x[0], depletion_voltage_data.values()))
    values_uncertainty = list(map(lambda x: x[1], depletion_voltage_data.values()))

    n_data = len(times)
    times = pyroot_array("d", times)
    times_uncertainty = pyroot_array("d", [3600 * 6] * n_data)  # 6 hours error
    values = pyroot_array("d", values)
    values_uncertainty = pyroot_array("d", values_uncertainty)

    graph = ROOT.TGraphErrors(
        n_data, times, values, times_uncertainty, values_uncertainty
    )

    graph.SetMarkerStyle(8)
    graph.SetMarkerSize(1)

    graph.GetXaxis().SetTimeDisplay(1)
    graph.GetXaxis().SetNdivisions(6, 2, 0)
    graph.GetXaxis().SetTimeFormat("%d/%m/%Y")
    graph.GetXaxis().SetTimeOffset(0, "gmt")

    graph.GetXaxis().SetTitle("Date")
    graph.GetYaxis().SetTitle("Depletion voltage [V]")
    graph.SetTitle("")

    return graph


def __plot_depletion_voltage_fit(
        output_directory,
        output_file_name,
        data,
        fits,
        legends,
        x_min=None,
        x_max=None,
        y_min=None,
        y_max=None,
    ):


    canvas = ROOT.TCanvas("", "", 800, 600)
    legend = ROOT.TLegend(0.14, 0.8, 0.86, 0.88)
    legend.SetNColumns(2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)

    if y_min is not None:
        data.SetMinimum(y_min)
    if y_max is not None:
        data.SetMaximum(y_max)

    if x_min is not None and x_max is not None:
        x_min_timestamp = time_string_to_timestamp(x_min, time_format="%Y-%m-%d")
        x_max_timestamp = time_string_to_timestamp(x_max, time_format="%Y-%m-%d")
        data.GetXaxis().SetLimits(x_min_timestamp, x_max_timestamp)
    if y_min is not None:
        data.SetMinimum(y_min)
    if y_max is not None:
        data.SetMaximum(y_max)

    data.Draw("AP SAME")
    legend.AddEntry(data, "Depletion Voltage Data", "lep")

    for fit, legend_label in zip(fits, legends):
        fit.Draw("SAME")
        legend.AddEntry(fit, legend_label, "l")
    
    legend.Draw("SAME")

    plot_name = f"{output_directory}/{output_file_name}"
    for extension in ["png", "pdf"]:
        canvas.SaveAs(f"{plot_name}.{extension}")



def main():
    args = __get_arguments()
    __check_arguments(args)

    Path(args.output_directory).mkdir(parents=True, exist_ok=True)

    depletion_voltage_graph = __get_depletion_voltage_graph(args.depletion_voltage_data_file_name)

    simulation_data = __read_data_from_simulation(args.input_simulation_file_name)

    v_depl = simulation_data["U_dep"].Clone()

    # hack for linear scaling
    #v_depl.Scale(0.44)
    #legends = ("Simulation (x 0.44)", )
    # Otherwise
    legends = ("Simulation", )

    v_depl.SetLineColor(ROOT.kBlue)
    v_depl.SetLineWidth(2)
    fits = (v_depl, )

    __plot_depletion_voltage_fit(
        args.output_directory,
        args.output_file_name,
        depletion_voltage_graph,
        fits,
        legends,
        x_min=args.x_min,
        x_max=args.x_max,
        y_min=args.y_min,
        y_max=args.y_max,
    )



if __name__ == "__main__":
    main()

