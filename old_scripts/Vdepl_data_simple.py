import os
import numpy as np
import math
import optparse
import csv
from datetime import datetime
import ROOT as rt
import collections
from array import array

# TODO: This is file is duplicated 5 times, need to find the "correct" one...

def mysort(list1, list2): #sort 2 lists by values in list1
    sorted_list1 = sorted(list1)
    sorted_list2 = [x for _, x in sorted(zip(list1, list2))]
    return sorted_list1, sorted_list2

def initDict(arg1,arg2):
    myDict = {}
    if isinstance(arg1, list) and isinstance(arg2, list):
        myDict = dict(zip(arg1,arg2))
    elif isinstance(arg1, list) and not isinstance(arg2, list):
        myDict = dict(zip(arg1,[arg2 for k in arg1]))
    else:
        print "The first argument should be a list!"
        exit()
    return myDict

class Vdepl:
    def __init__(self, BEGIN_RUN):
        self.begin_run = BEGIN_RUN
        self.date_format = "%Y-%m-%d %H:%M:%S"

    # def get_start_run(self, BEGIN_RUN):
    #     self.begin_run = BEGIN_RUN

    def getFile(self, fileName):
        self.input_file = open(fileName, "r")

    def get_csv(self):
        self.csv_obj = csv.DictReader(self.input_file)

    def get_keys(self):
        keys = self.csv_obj.fieldnames
        self.keys = [i for i in keys if i and i != "Date"]

    def vdepl_data_init(self):
        dates_Vdepl = [{"Date": [], "Vdepl": []} for k in self.keys]
        self.vdepl = initDict(self.keys, dates_Vdepl)

    def initDataStruct(self,fileName):
        self.getFile(fileName)
        self.get_csv()
        self.get_keys()
        self.vdepl_data_init()

    def get_date(self, date_str):
        myDate = datetime.now()
        if date_str is not '':
            if date_str == 'Start':
                myDate = datetime.strptime(self.begin_run, self.date_format)
            else:
                dateList = date_str.split(".")
                dateList = map(lambda x: int(x), dateList)
                date_time = "%4d-%02d-%02d 12:00:00"%(dateList[2], dateList[1], dateList[0])
                myDate = datetime.strptime(date_time, self.date_format)
            return myDate
        else:
            return "You should add \"if date_str is not ''\" before using this function!"
            exit()

    def get_dates_vdepl(self):
        self.nPoints = initDict(self.keys, None)
        for d in self.csv_obj:
            date_str_from_csv = d["Date"]
            for k in self.keys:
                vdepl_str_from_csv = d[k]
                if date_str_from_csv is not '' and vdepl_str_from_csv is not '':
                    dataPoint_date = self.get_date(date_str_from_csv)
                    self.vdepl[k]["Date"].append(dataPoint_date)
                    self.vdepl[k]["Vdepl"].append(vdepl_str_from_csv.rstrip("?"))
        for k in self.keys:
            self.vdepl[k]["Date"], self.vdepl[k]["Vdepl"] = mysort(self.vdepl[k]["Date"], self.vdepl[k]["Vdepl"])
            self.vdepl[k]["Date"] = map(lambda x: float(rt.TDatime(datetime.strftime(x, self.date_format)).Convert()), self.vdepl[k]["Date"])
            self.vdepl[k]["Vdepl"] = map(lambda x: float(x), self.vdepl[k]["Vdepl"])
            self.vdepl[k]["Date"] = array('d',self.vdepl[k]["Date"])
            self.vdepl[k]["Vdepl"] = array('d',self.vdepl[k]["Vdepl"])
            self.nPoints[k] = len(self.vdepl[k]["Date"])

    def init_tgraphs(self, chooseKeys):
        self.keys = chooseKeys
        self.tg_dict = initDict(self.keys,rt.TGraph())

    def create_tgraphs(self):
        for k, tgr in self.tg_dict.iteritems():
            self.tg_dict[k] = rt.TGraph(self.nPoints[k], self.vdepl[k]["Date"], self.vdepl[k]["Vdepl"])

    def set_draw_options(self,drawOptionsFile):
        # should be smth else, taken from config file
        # creating self.draw_options struct
        return 0

    def apply_draw_options(self):
        # should be smth else using self.draw_options
        self.create_tgraphs()
        print self.tg_dict_sim
        self.colors = initDict(self.keys, [rt.kRed,rt.kBlue,rt.kMagenta,rt.kAzure,rt.kBlack,rt.kPink,rt.kCyan])
        print self.colors
        print self.tg_dict
        print self.tg_dict_sim
        for k in self.keys:
            self.tg_dict[k].SetMarkerStyle(20)
            self.tg_dict[k].SetMarkerColor(self.colors[k])
            self.tg_dict_sim[k].SetLineColor(self.colors[k])
            self.tg_dict_sim[k].SetLineWidth(2)

    def init_vdepl_sim(self, outputSimFile):
        self.simFile = rt.TFile.Open(outputSimFile)
        self.tg_dict_sim = initDict(self.keys,rt.TGraph())
        print self.tg_dict_sim

    def get_vdepl_sim(self):
        # future
        # for k,tgr in self.tg_dict_sim.iteritems():
        #     self.tg_dict_sim[k] = rt.Get("U_dep_%s"%k)
        #now
        k="LAYER 1"
        self.tg_dict_sim[k] = self.simFile.Get("U_dep")
        self.tg_dict_sim[k].SetTitle("U_{depletion} for %s"%k)

    def draw_graphs(self,canv,vdepl,vdepl_sim):
        canv.Draw()
        vdepl_sim.Draw("AL")
        vdepl.Draw("P")

    def save_vdepl(self):
        outputFile = rt.TFile.Open("Vdepl.root","RECREATE")
        # future
        # for value in variable:
            # pass
        # now
        k="LAYER 1"
        self.canv1 = rt.TCanvas("Vdepl_%s"%k.replace(" ", ""),"Vdepl_%s"%k.replace(" ", ""),1600,800)
        self.draw_graphs(self.canv1,self.tg_dict[k],self.tg_dict_sim[k])
        self.canv1.Write()
        outputFile.Close()


def main():
    BEGIN_RUN = "2017-05-23 14:32:22"
    fileName="Vdepl_data/CMS_PIXELS_VDEPL.csv"
    outputSimFile="simulation_results_test_new_datetime.root"
    chooseKeys=["LAYER 1"]

    myVdepl = Vdepl(BEGIN_RUN)

    myVdepl.initDataStruct(fileName)
    myVdepl.get_dates_vdepl()
    myVdepl.init_vdepl_sim(outputSimFile)
    myVdepl.get_vdepl_sim()
    myVdepl.init_tgraphs(chooseKeys)
    myVdepl.apply_draw_options()
    myVdepl.save_vdepl()

if __name__ == "__main__":
    main()
