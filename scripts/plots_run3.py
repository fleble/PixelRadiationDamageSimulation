import ROOT as rt
import argparse

def make_plots(lumi,end_year,variable,histo_name_sim,histo_name_data,add_log,yaxis_title,tgr,title,disk,ring):
    files_arr = []
    for i in range(2,5):
        if variable == "vdepl":
            files_arr.append(rt.TFile.Open("Vdepl_sim_r%sd%srog1BpO_clean_run3_wp%sm_run21%s_l%s.root"%(ring,disk,i,end_year,lumi)))
        else:
            files_arr.append(rt.TFile.Open(
                "sim_r%sd%srog1BpO_ileak_clean_run3_wp%sm_run21%s_l%s.root" % (ring, disk, i, end_year, lumi)))
    files_arr_log = []
    
    if variable == "vdepl" and add_log:
        for i in range(2,5):
            files_arr_log.append(rt.TFile.Open(
                "Vdepl_sim_r%sd%srog1BpO_clean_run3_log_wp%sm_run21%s_l%s.root" % (ring, disk, i, end_year, lumi)))
    
    colors = [rt.kBlack,rt.kGreen+3,rt.kRed,rt.kBlue]
    histos_sim = []
    histos_sim_log = []
    for f in files_arr:
        histos_sim.append(f.Get(histo_name_sim))
    for f in files_arr_log:
        histos_sim_log.append(f.Get(histo_name_sim))
    if histo_name_data is not "":
        histo_data = files_arr[0].Get(histo_name_data)
        if tgr:
            histo_data.SetFillStyle(0)
    leg = rt.TLegend(0.1,0.6,0.4,0.9)
    # leg.SetFillStyle(0)
    # leg.SetHeader("Warm period duration:")
    if histo_name_data is not "": leg.AddEntry(histo_data,"Data, Run2")
    
    for i,h in enumerate(histos_sim):
        h.SetLineColor(colors[i])
        h.SetLineWidth(1)
        h.SetMarkerSize(0)
        h.SetMarkerColor(colors[i])
        h.SetTitle("%s, %s fb^{-1}"%(title,lumi))
        h.GetXaxis().SetTitle("Date")
        h.GetYaxis().SetTitle(yaxis_title)
        h.SetFillStyle(0)
        leg.AddEntry(h,"Warm period: %s month(s)"%(i+2))
    
    for i,h in enumerate(histos_sim_log):
        h.SetLineColor(colors[i])
        h.SetMarkerColor(colors[i])
        h.SetLineWidth(1)
        h.SetLineStyle(3)
        h.SetMarkerSize(0)
        h.SetTitle("Depletion voltage, %s fb^{-1}" % lumi)
        h.GetXaxis().SetTitle("Date")
        h.GetYaxis().SetTitle(yaxis_title)
        leg.AddEntry(h, "Warm period: %s month(s), log" % (i+2))
    if histo_name_data is not "":
        histo_data.SetMarkerStyle(20)
        histo_data.SetMarkerSize(1)
        histo_data.SetMarkerColor(rt.kViolet)

    if add_log:
        outfile = rt.TFile.Open("%s_%s_%s_r%sd%s_run3_log.root" % (
            variable, ring, disk, lumi, histo_name_sim), "RECREATE")
    else:
        outfile = rt.TFile.Open("%s_%s_%s_r%sd%s_run3.root" % (
            variable, ring, disk, lumi, histo_name_sim), "RECREATE")
    canv = rt.TCanvas(variable,variable,0,0,1600,800)
    canv.Draw()
    print histos_sim
    if histo_name_sim == "fluence":
        histos_sim[0].Draw("aly")
    else:
        for i in range(2,-1,-1):
            if histo_name_sim == "cp_Temperature": histos_sim[i].GetYaxis().SetRangeUser(250,330)
            if i == 2 and tgr: histos_sim[i].Draw("aly")
            else: histos_sim[i].Draw("same l y")
            if add_log: histos_sim_log[i].Draw("same l")
    if histo_name_data is not "": histo_data.Draw("same p")
    
    if histo_name_sim is not "fluence": leg.Draw("same")
    rt.gPad.SetTicky()
    canv.Write()
    if add_log: canv.Print("vdepl_run3_pdfs_r%sd%s/%s_%s_%s_run3_log.pdf" % (ring,disk,variable, lumi, histo_name_sim))
    else: canv.Print("vdepl_run3_pdfs_r%sd%s/%s_%s_%s_run3.pdf" % (ring,disk,variable, lumi, histo_name_sim))
    outfile.Close()

def main(args):
    disk = args.disk
    ring = args.ring
    yaxis_title = "Depletion voltage [V]"
    title = "Depletion voltage"
    # make_plots("200","23","vdepl","Func","Graph",False,yaxis_title,False,title,disk,ring)
    make_plots("250","23","vdepl","Func","Graph",False,yaxis_title,False,title,disk,ring)
    make_plots("350","24","vdepl","Func","Graph",False,yaxis_title,False,title,disk,ring)
    # make_plots("200","23","vdepl","Func","Graph",True,yaxis_title,False,title,disk,ring)
    make_plots("250","23","vdepl","Func","Graph",True,yaxis_title,False,title,disk,ring)
    make_plots("350","24","vdepl","Func","Graph",True,yaxis_title,False,title,disk,ring)
    
    yaxis_title = "I_{leak} (@ -7C) [mA], 1 ROG"
    title = "Leakage current"
    # make_plots("200","23","ileak","I_leak_per_module","I_leak_per_module_data",False,yaxis_title,True,title,disk,ring)
    make_plots("250","23","ileak","I_leak_per_module","I_leak_per_module_data",False,yaxis_title,True,title,disk,ring)
    make_plots("350","24","ileak","I_leak_per_module","I_leak_per_module_data",False,yaxis_title,True,title,disk,ring)

    yaxis_title = "Fluence [N/cm^{2}]"
    title = "Fluence"
    # make_plots("200","23","ileak","fluence","",False,yaxis_title,True,title,disk,ring)
    make_plots("250","23","ileak","fluence","",False,yaxis_title,True,title,disk,ring)
    make_plots("350","24","ileak","fluence","",False,yaxis_title,True,title,disk,ring)

    # yaxis_title = "Temperature [K]"
    # title = "Temperature"
    # make_plots("200","23","ileak","cp_Temperature","",False,yaxis_title,True,title,disk,ring)
    make_plots("250","23","ileak","cp_Temperature","",False,yaxis_title,True,title,disk,ring)
    make_plots("350","24","ileak","cp_Temperature","",False,yaxis_title,True,title,disk,ring)


# I_leak_per_module
# I_leak_per_module_data
# fluence
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--batch", dest='batch', action="store_true")
    parser.add_argument("-d","--disk", dest='disk', action="store", type=int)
    parser.add_argument("-r", "--ring", dest='ring',
                        action="store", type=int)
    args = parser.parse_args()
    main(args)
