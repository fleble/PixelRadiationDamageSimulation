#!/usr/bin/env python
import math
from array import array
import ROOT as rt
from datetime import datetime
from datetime import timedelta

##########################################################################
## Data points: Depletion voltage
##########################################################################

# Manually read-off kink from HV scans as depletion voltage data point (count days from 2017-05-23):
# TODO: Properly extract values with JL Agrams script (use curvature method for highly irradiated)
# 2017-08-14
# 2017-09-23
# 2017-11-09
# 2018-04-19 full scan (until 600V) // 331
# 2018-05-05 mini bias scan // 347 //L1: weird Ileak values
# 2018-05-12 mini bias scan // 354 //L1: weird Ileak values
# 2018-05-25 mini bias scan // 367 //L1: weird Ileak values
# 2018-06-09 mini bias scan // 382 //L1: weird Ileak values
# 2018-07-11 mini bias scan // 414 //L1: weird Ileak values
# 2018-07-31 full scan (until 400V) // 433
# 2018-08-17 mini scan (until 600V) // 451
# 2018-09-01 mini scan (until 600V) // 466
# 2018-09-07 mini scan (until 600V) new hvgrp // 472
# 2018-09-26 mini scan // 491

#READING OFF -- from cluster charge and cluster size

BEGIN_RUN  = "2017-05-23 12:00:00"
FORMAT_DATE = "%Y-%m-%d %H:%M:%S"
daydata_L1      = [ 123*24*60*60, 170*24*60*60, 331*24*60*60, 347*24*60*60, 354*24*60*60, 433*24*60*60, 451*24*60*60, 466*24*60*60, 472*24*60*60, 491*24*60*60]
lumdata_L1      = [ 24.83, 49.93, 51.67, 56.3, 60.27, 84.83, 94.37, 100.67, 103.85, 109.40 ]
CluChargeDataL1 = [370., 470., 240., 280., 330., 370., 390., 400., 500., 500.]
CluSizDataL1    = [ 260., 440., 200., 230., 260., 300., 390., 420., 450., 450. ]
daydata_L2      = [83*24*60*60, 123*24*60*60, 170*24*60*60, 331*24*60*60, 347*24*60*60, 354*24*60*60, 367*24*60*60, 382*24*60*60, 414*24*60*60, 433*24*60*60, 451*24*60*60, 466*24*60*60, 472*24*60*60, 491*24*60*60]
lumdata_L2      = [ 15.93, 24.83, 49.93, 51.67, 56.3, 60.27, 67.02, 74.9706, 82.92, 84.83, 94.37, 100.67, 103.85, 109.40]
CluChargeDataL2 = [87.,175.,230.,110.,130.,150.,180.,200.,200.,150.,195.,190.,220.,220.]
CluSizDataL2    = [85.,85.,200.,110.,130.,145.,170.,210.,175.,140.,210.,220.,210.,210.]
daydata_L3      =   [83*24*60*60, 123*24*60*60, 170*24*60*60, 331*24*60*60, 347*24*60*60, 354*24*60*60, 367*24*60*60, 382*24*60*60, 414*24*60*60, 433*24*60*60, 451*24*60*60, 466*24*60*60, 472*24*60*60, 491*24*60*60]
lumdata_L3      = [ 15.93, 24.83, 49.93, 51.67, 56.3, 60.27, 67.02, 74.9706, 82.92, 84.83, 94.37, 100.67, 103.85, 109.40 ]
CluChargeDataL3 = [ 72., 120., 150., 60., 75., 90., 115., 120., 110., 85., 130., 140., 140., 150. ]
CluSizDataL3    = [ 50., 50., 140., 70., 80., 90., 105., 125., 125., 100., 125., 130., 150., 150. ]
daydata_L4      = [83*24*60*60, 123*24*60*60, 170*24*60*60, 331*24*60*60, 347*24*60*60, 354*24*60*60, 367*24*60*60, 382*24*60*60, 414*24*60*60, 433*24*60*60, 451*24*60*60, 466*24*60*60, 472*24*60*60, 491*24*60*60]
lumdata_L4      = [15.93,24.83,49.93,51.67,56.3,60.27,67.02,74.9706,82.92,84.83,94.37,100.67,  103.85,109.40]
CluChargeDataL4 = [60.,110.,140.,50.,55.,60.,65.,85.,70.,60.,80.,90.,120.,120.]
CluSizDataL4    = [45.,45.,100.,50.,55.,60.,65.,85.,100.,70.,85.,100.,110.,110.]
dataentries     = [len(lumdata_L1), len(lumdata_L2), len(lumdata_L3), len(lumdata_L4)]


#CURVATURE FIT -- cluster charge
#L1 fit failed always
dayCurvdata_L2      = [347*24*60*60, 414*24*60*60, 433*24*60*60, 451*24*60*60, 466*24*60*60, 472*24*60*60]
CluChargeCurvDataL2 = [122.233, 163.814, 145.80, 179.65, 190.57, 198.00]
dayCurvdata_L3      = [331*24*60*60, 347*24*60*60, 414*24*60*60, 433*24*60*60, 451*24*60*60, 466*24*60*60, 472*24*60*60]
CluChargeCurvDataL3 = [61.39, 71.32, 99.05, 84.40, 103.64, 125.22, 137.03]
dayCurvdata_L4      = [347*24*60*60, 414*24*60*60, 433*24*60*60, 451*24*60*60, 466*24*60*60, 472*24*60*60]
CluChargeCurvDataL4 = [44.32, 58.64, 57.56, 74.02, 82.58, 92.57]
dataCurventries     = [len(dayCurvdata_L2), len(dayCurvdata_L3), len(dayCurvdata_L4)]

#CURVATURE FIT -- cluster size
dayCurvCSdata_L2  = [ 347*24*60*60, 433*24*60*60, 451*24*60*60, 466*24*60*60]
CluSizCurvDataL2  = [ 110.82,       159.88,       208.60,       229.90]
dayCurvCSdata_L3  = [ 331*24*60*60, 347*24*60*60, 414*24*60*60, 433*24*60*60, 451*24*60*60, 466*24*60*60, 472*24*60*60]
CluSizCurvDataL3  = [ 63.43,        74.02,        92.64,        85.36,       102.92,       115.26,       140.55]
dayCurvCSdata_L4  = [ 347*24*60*60, 414*24*60*60, 433*24*60*60, 451*24*60*60]
CluSizCurvDataL4  = [ 45.78,        52.51,        60.81,        52.32]
dataCurvCSentries = [len(dayCurvdata_L2), len(dayCurvdata_L3), len(dayCurvdata_L4)]

def createArray(myList):
    return array('d', map(lambda x: float(x), myList))

def create_tgraph(x,y):
    return rt.TGraph(len(x),x,y)

def load_tcanv(fileName,canvName):
    tfile = rt.TFile.Open(fileName)
    canv1 = tfile.Get(canvName)
    return canv1

def draw_all(canv1, tgraph, outputFile):
    outputFile = rt.TFile.Open(outputFile,"RECREATE")
    canv1.Draw()
    tgraph[0].SetMarkerStyle(20)
    tgraph[0].SetFillColor(0)
    tgraph[0].SetLineColor(0)
    tgraph[0].SetMarkerColor(rt.kBlue)
    tgraph[0].SetName("UdepCluChargeDataL1")
    tgraph[0].SetTitle("U_{dep}, from Cluster Charge L1")
    tgraph[0].Draw("P")
    tgraph[1].SetMarkerStyle(20)
    tgraph[1].SetFillColor(0)
    tgraph[1].SetLineColor(0)
    tgraph[1].SetMarkerColor(rt.kGreen+3)
    tgraph[1].SetName("UdepCluSizDataL1")
    tgraph[1].SetTitle("U_{dep}, from Cluster Size L1")
    tgraph[1].Draw("P")
    draw_legend(canv1,tgraph)
    canv1.Write()

def get_seconds(begin_date, seconds_after_start):
    bd = datetime.strptime(begin_date,FORMAT_DATE)
    d = bd + timedelta(0,seconds_after_start)
    return float(rt.TDatime(d.strftime(FORMAT_DATE)).Convert())

def date_array(begin_date,date_list):
    return [ get_seconds(begin_date,d) for d in date_list ]


def draw_legend(canv1,gr):
    leg = canv1.GetListOfPrimitives().FindObject("TPave")
    for i,g in enumerate(gr):
        leg.AddEntry(g,g.GetTitle())
    leg.Draw()

def main():
    # daydata_L1
    # CluChargeDataL1
    # CluSizDataL1
    fileName = "Vdepl.root"
    canvName = "Vdepl_LAYER1"
    xL1 = date_array(BEGIN_RUN,daydata_L1)
    xL1 = createArray(xL1)
    yL1 = createArray(CluChargeDataL1)
    yL2 = createArray(CluSizDataL1)
    gr = [create_tgraph(xL1,yL1),create_tgraph(xL1,yL2)]
    canv1 = load_tcanv(fileName,canvName)
    draw_all(canv1,gr,"Vdepl_plusfinn.root")

if __name__ == "__main__":
    main()
