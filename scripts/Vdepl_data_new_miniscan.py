import os
import numpy as np
import math
import optparse
import csv
from datetime import datetime
from datetime import timedelta
import ROOT as rt
import collections
from array import array
import numpy
from fluka_l1 import fl_pos_dict
from get_fluence import GetFlunce

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

def calculate_chi2(data_tgraphasymerr, sim_tf1, t0):
    chi2 = 0
    n_points = 0
    # print data_tgraphasymerr.GetN()
    x=rt.Double(0)
    y=rt.Double(0)
    for point in range(0,data_tgraphasymerr.GetN()):
        data_tgraphasymerr.GetPoint(point,x,y)
        if x > t0:
            n_points += 1
            chi2 += (y-sim_tf1.Eval(x))**2/(data_tgraphasymerr.GetErrorYhigh(point)**2)
            # print "%s :: %s: fit: %4.2f -> %4.2f :tgraph: unc: %4.2f"%(point, x,sim_tf1.Eval(x),y,data_tgraphasymerr.GetErrorYhigh(point))
        else:
            pass #print "Why not %s: %s, %s"%(point, x,y)
    # print "Chi2(root) = %s"%data_tgraphasymerr.Chisquare(sim_tf1)
    return [chi2, n_points]

# def vDeplFitFunc(x, pars):
    # return eval

class Vdepl:
    def __init__(self, BEGIN_RUN, disk, ring, is_log, is_all, plot_add):
        self.begin_run = BEGIN_RUN
        self.date_format = "%Y-%m-%d %H:%M:%S"
        self.isfit = True
        self.log_dependence = is_log
        self.disk = disk
        self.ring = ring
        self.is_all = is_all
        self.vdepl_add = {}
        self.plot_add = plot_add

    # def get_start_run(self, BEGIN_RUN):
    #     self.begin_run = BEGIN_RUN

    def getFile(self, fileName):
        self.input_file = open(fileName, "r")

    def get_csv(self):
        self.csv_obj = csv.DictReader(self.input_file)

    def get_keys(self,mykeys):
        # mykeys = 'l1,l2,l3,l4,d1,d2,d3', 'allbpix', 'allfpix'
        keys = self.csv_obj.fieldnames
        if mykeys == 'all':
            self.keys = [i for i in keys if i and i != "Date"]
        elif mykeys == 'allbpix':
            self.keys = [i for i in keys if "LAYER" in i]
        elif mykeys == 'allfpix':
            self.keys = [i for i in keys if "FDISK" in i]
        else:
            self.keys = [k1.replace("l","LAYER ").replace("d", "FDISK") for k1 in mykeys.split(",")]
        self.choose_key = self.keys[0]
        # print self.keys

    def vdepl_data_init(self):
        self.keys = ["LAYER1"]
        dates_Vdepl = [{"Date": [], "Vdepl": [], "eVdepl": []} for k in self.keys]
        self.vdepl = initDict(self.keys, dates_Vdepl)
    
    def vdepl_add_data_init(self, mods):
        # self.keys = ["LAYER1"]
        self.mods_to_plot = mods
        for k in self.keys:
            self.vdepl_add[k] = {}
            for modx in mods:
                self.vdepl_add[k][modx] = {}
                self.vdepl_add[k][modx]["Date"]   = list()
                self.vdepl_add[k][modx]["Vdepl"]  = list()
                self.vdepl_add[k][modx]["eVdepl"] = list()

    def initDataStructFromCSV(self,fileName,mykeys='all'):
        self.getFile(fileName)
        self.get_csv()
        self.get_keys(mykeys)
        self.vdepl_data_init()
    
    def get_vdepl_from_fit(self, vdepl_date_time, tag):
        for d in vdepl_date_time:
            self.vdepl[tag]["Date"].append(float(rt.TDatime(d).Convert()))
            self.vdepl[tag]["Vdepl"].append(vdepl_date_time[d][0])
            self.vdepl[tag]["eVdepl"].append(vdepl_date_time[d][1])
        # print self.vdepl[tag]["Date"]
        # print self.vdepl[tag]["Vdepl"]
        # print self.vdepl[tag]["eVdepl"]
        # map(lambda x: float(rt.TDatime(datetime.strftime(
        #     x, self.date_format)).Convert()), self.vdepl[k]["Date"])
        self.vdepl[tag]["Date"] = array('d', self.vdepl[tag]["Date"])
        self.vdepl[tag]["Vdepl"] = array('d', self.vdepl[tag]["Vdepl"])
        self.vdepl[tag]["eVdepl"] = array('d', self.vdepl[tag]["eVdepl"])
        # self.vdepl[k]["Date"], self.vdepl[k]["Vdepl"] = mysort(self.vdepl[k]["Date"], self.vdepl[k]["Vdepl"])
        # self.vdepl[k]["Date"], self.vdepl[k]["Vdepl"] = mysort(self.vdepl[k]["Date"], self.vdepl[k]["eVdepl"])
    
    def get_vdepl_additional(self, vdepl_date_time, tag, modx):
        for d in vdepl_date_time:
            # self.vdepl_add[tag][modx]["Date"].append(float(rt.TDatime(d).Convert()))
            # self.vdepl_add[tag][modx]["Vdepl"].append(vdepl_date_time[d][0])
            # self.vdepl_add[tag][modx]["eVdepl"].append(vdepl_date_time[d][1])
            self.vdepl_add[tag][modx]["Date"] = [(float(rt.TDatime(d).Convert()))]
            self.vdepl_add[tag][modx]["Vdepl"] = [(vdepl_date_time[d][0])]
            self.vdepl_add[tag][modx]["eVdepl"] = [(vdepl_date_time[d][1])]
        # print self.vdepl_add[tag]["Date"]
        # print self.vdepl_add[tag]["Vdepl"]
        # print self.vdepl_add[tag]["eVdepl"]
        # map(lambda x: float(rt.TDatime(datetime.strftime(
        #     x, self.date_format)).Convert()), self.vdepl[k]["Date"])
        self.vdepl_add[tag][modx]["Date"]   = array('d', self.vdepl_add[tag][modx]["Date"])
        self.vdepl_add[tag][modx]["Vdepl"]  = array('d', self.vdepl_add[tag][modx]["Vdepl"])
        self.vdepl_add[tag][modx]["eVdepl"] = array('d', self.vdepl_add[tag][modx]["eVdepl"])

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

    def init_tgraphs(self):
        self.tg_dict = initDict(self.keys,rt.TGraphErrors())
        if self.plot_add:
            self.tg_dict_add = {}
            for k in self.keys:
                self.tg_dict_add[k] = initDict(self.mods_to_plot, rt.TGraphErrors())


    def create_tgraphs(self):
        # syst_err_func = lambda x: 0.1*x #10% of Vdepl
        syst_err_func = lambda x: 10 #Volts
        syst_err_func_modx = lambda x: 10. #Volts
        stat_errors = {
            "LAYER1":{
                "eyup": self.vdepl['LAYER1']["eVdepl"], # array('d',[0,30,30,30,30,30,30,30,30,30,30,30,30,30]),
                "eylow": self.vdepl['LAYER1']["eVdepl"], # array('d',[0,25,25,25,25,25,25,25,20,20,20,20,20,20])
            },
            "FDISK1":{
                "eyup": array('d',[0,30,30,30,30,30,30,30,30,30,30,30,30,30]),
                "eylow": array('d',[0,25,25,25,25,25,25,25,20,20,20,20,20,20])
            },
            "FDISK2":{
                "eyup": array('d',[0,30,30,30,30,30,30,30,30,30]),
                "eylow": array('d',[0,20,20,20,20,20,20,20,20,20])
            },
            # "FDISK3":{
            #     "eyup": []
            #     "eylow":[]
            # },
        }
        stat_errors_add = {
            "LAYER1":{
                # "mod1":{
                #     "eyup":  self.vdepl_add['LAYER1']['mod1']["eVdepl"],
                #     "eylow": self.vdepl_add['LAYER1']['mod1']["eVdepl"],
                # },
                "mod2":{
                    "eyup":  self.vdepl_add['LAYER1']['mod2']["eVdepl"],
                    "eylow": self.vdepl_add['LAYER1']['mod2']["eVdepl"],
                },
                # "mod3":{
                #     "eyup":  self.vdepl_add['LAYER1']['mod3']["eVdepl"],
                #     "eylow": self.vdepl_add['LAYER1']['mod3']["eVdepl"],
                # },
                "mod4":{
                    "eyup":  self.vdepl_add['LAYER1']['mod4']["eVdepl"],
                    "eylow": self.vdepl_add['LAYER1']['mod4']["eVdepl"],
                }
            }
        }
        # nP = self.nPoints[k]
        nP = len(self.vdepl["LAYER1"]["Date"])
        for k in self.keys:
            # ex=array('d',[0]*self.nPoints[k])
            # ey=array('d',[15]*self.nPoints[k])
            # self.tg_dict[k] = rt.TGraphErrors(self.nPoints[k], self.vdepl[k]["Date"], self.vdepl[k]["Vdepl"], ex, ey)
            syst_errors = map(syst_err_func, self.vdepl[k]["Vdepl"])
            # print syst_errors
            errors = [math.sqrt(syst_errors[i]**2 + stat_errors[k]["eylow"][i]**2) for i,_ in enumerate(syst_errors)]
            errors = array('d', errors)
            self.tg_dict[k] = rt.TGraphAsymmErrors(
                nP, self.vdepl[k]["Date"], self.vdepl[k]["Vdepl"], array('d', [0]*nP), array('d', [0]*nP), errors, errors)
            # self.tg_dict[k] = rt.TGraph(self.nPoints[k], self.vdepl[k]["Date"], self.vdepl[k]["Vdepl"])
            for modx in self.mods_to_plot:
                syst_errors_add = map(syst_err_func_modx, self.vdepl_add[k][modx]["Vdepl"])
                errors_add = [math.sqrt(syst_errors_add[i]**2 + stat_errors_add[k][modx]["eylow"][i]**2) for i,_ in enumerate(syst_errors_add)]
                errors_add = array('d', errors_add)
                # print "self.vdepl_add[%s][%s][\"Date\"] -> %s"%(k,modx,self.vdepl_add[k][modx]["Date"])
                # print "self.vdepl_add[%s][%s][\"Vdepl\"] -> %s"%(k,modx,self.vdepl_add[k][modx]["Vdepl"])
                self.tg_dict_add[k][modx] = rt.TGraphAsymmErrors(1, self.vdepl_add[k][modx]["Date"], self.vdepl_add[k][modx]["Vdepl"], array('d', [0]), array('d', [0]), errors_add, errors_add)
                print(self.tg_dict_add[k][modx].GetN())
                print("what?")

    def set_draw_options(self,drawOptionsFile):
        # should be smth else, taken from config file
        # creating self.draw_options struct
        return 0

    def apply_draw_options(self):
        # should be smth else using self.draw_options
        self.create_tgraphs()
        isfit = True
        # print self.tg_dict_sim
        self.colors = initDict(self.keys, [
                               rt.kRed, rt.kBlue, rt.kMagenta, rt.kAzure, rt.kBlack, rt.kRed, rt.kRed])
        # print self.colors
        # print self.tg_dict
        # print self.tg_dict_sim
        for k in self.keys:
            self.tg_dict[k].SetMarkerStyle(20)
            # self.tg_dict[k].SetLineColor(0)
            self.tg_dict[k].SetMarkerColor(self.colors[k])
            self.tg_dict[k].GetXaxis().SetTimeDisplay(1)
            self.tg_dict[k].GetXaxis().SetNdivisions(6,2,0)
            self.tg_dict[k].GetXaxis().SetTimeFormat("%d/%m/%Y")
            self.tg_dict[k].GetXaxis().SetTimeOffset(0,"gmt")
            if not self.isfit:
                self.tg_dict_sim[k].SetLineColor(self.colors[k])
                self.tg_dict_sim[k].SetMarkerColor(self.colors[k])
                self.tg_dict_sim[k].SetLineWidth(2)
            else:
                self.tg_dict_sim[k].SetLineWidth(0)
                self.tg_dict_sim[k].SetMarkerColor(0)
            # remove black boxes from TLegend
            self.tg_dict[k].SetFillColor(0)
            self.tg_dict_sim[k].SetFillColor(0)

    def init_vdepl_sim(self, outputSimFile):
        self.simFile = rt.TFile.Open(outputSimFile)
        self.tg_dict_sim = initDict(self.keys,rt.TGraph())
        # print self.tg_dict_sim

    def get_vdepl_sim(self):
        # future
        # for k,tgr in self.tg_dict_sim.iteritems():
        #     self.tg_dict_sim[k] = rt.Get("U_dep_%s"%k)
        #now
        k = self.choose_key
        self.tg_dict_sim[k] = self.simFile.Get("U_dep")
        self.tg_dict_sim[k].SetTitle("Depletion voltage")

    def get_Vdepl_terms(self):
        # self.N_benef_anneal_g1 = []
        # self.N_revers_anneal_g1 = []
        # self.N_nadefects_g1 = []
        # self.N_constdamage_g1 = []
        # self.N_donor = []
        # for i in range(0,64):
        #     self.N_benef_anneal_g1.append(self.simFile.Get("N_benef_anneal_g1_"+str(i)))#get_gA
        #     self.N_revers_anneal_g1.append(self.simFile.Get("N_revers_anneal_g1_"+str(i)))#get_gY
        #     self.N_nadefects_g1.append(self.simFile.Get("N_nadefects_g1_"+str(i)))#get_gY
        #     self.N_constdamage_g1.append(self.simFile.Get("N_constdamage_g1_"+str(i)))#get_gC
        #     self.N_donor.append(self.simFile.Get("N_donor_"+str(i)))
        #     self.fluence = self.simFile.Get("fluence"+str(i))
        
        self.N_benef_anneal_g1 = self.simFile.Get("N_benef_anneal_g1")#get_gA
        self.N_revers_anneal_g1 = self.simFile.Get("N_revers_anneal_g1")#get_gY
        self.N_nadefects_g1 = self.simFile.Get("N_nadefects_g1")#get_gY
        self.N_constdamage_g1 = self.simFile.Get("N_constdamage_g1")#get_gC
        self.fluence = self.simFile.Get("fluence")
        self.N_donor = self.simFile.Get("N_donor")
        self.N_donor_stable_donorremoval = self.simFile.Get("N_donor_stable_donorremoval")

    def fsum(self, x, p0, p1, p2, N_donor, N_benef_anneal_g1, N_constdamage_g1, N_revers_anneal_g1):
        ret=0
        for i in range(0,64):
            ret +=  abs(-N_donor[i].Eval(x) +
                    p0*N_benef_anneal_g1[i].Eval(x) +
                    p1*N_constdamage_g1[i].Eval(x) +
                    p2*N_revers_anneal_g1[i].Eval(x))
        return 0.1/6.4*ret
    
    def fmiddle(self, x, p0, p1, p2, N_donor, N_benef_anneal_g1, N_constdamage_g1, N_revers_anneal_g1):
        index = 0
        ret =  abs(-N_donor[index].Eval(x) +
                p0*N_benef_anneal_g1[index].Eval(x) +
                p1*N_constdamage_g1[index].Eval(x) +
                p2*N_revers_anneal_g1[index].Eval(x))
        return ret

    def flog(self, x, p0, p1, p2, N_donor, N_benef_anneal_g1, N_constdamage_g1, N_revers_anneal_g1):
        index = 0
        ret = abs(-N_donor[index].Eval(x) +
                p0*N_benef_anneal_g1[index].Eval(x) +
                p1*math.log(abs(N_constdamage_g1[index].Eval(x))+1) +
                p2*N_revers_anneal_g1[index].Eval(x))
        return ret
    
    def flog_l1(self, x, p0, p1, p2, N_donor, N_benef_anneal_g1, fluence, N_revers_anneal_g1):
        ret = abs(-N_donor.Eval(x) +
            p0*N_benef_anneal_g1.Eval(x) +
            p1*math.log(fluence.Eval(x)+1) +
            p2*N_revers_anneal_g1.Eval(x))
        return ret
    
    def falphalog_l1(self, x, p0, p1, p2, p3, corr_phi, N_donor, N_benef_anneal_g1, fluence, N_revers_anneal_g1):
        ret = abs(-corr_phi*p3*N_donor.Eval(x) +
            corr_phi*p3*p0*N_benef_anneal_g1.Eval(x) +
            p1*math.log(corr_phi*p3*fluence.Eval(x)+1) +
            corr_phi*p3*p2*N_revers_anneal_g1.Eval(x))
        return ret
    
    def ndonor(self, x, p0, N_donor):
        ret = p0*N_donor.Eval(x)
        return ret
    def nbaneal(self, x, p0, N_benef_anneal_g1):
        ret = p0*N_benef_anneal_g1.Eval(x)
        return ret
    def nraneal(self, x, p0, N_revers_anneal_g1):
        ret = p0*N_revers_anneal_g1.Eval(x)
        return ret
    def nconst(self, x, p0, fluence):
        ret = p0*fluence.Eval(x)
        return ret
    def nconst_log(self, x, p0, fluence):
        ret = p0*math.log(fluence.Eval(x)+1)
        return ret
    
    def plot_constituents(self):
        thickness = 285.
        multiplier = 1.6021766208e-13/(2*11.68*8.854187817)*thickness*thickness
        multiplier = 1.
        self.__ndonor__ = lambda x, pars: multiplier*\
            self.ndonor(x[0], pars[0], self.N_donor)
        self.__nbaneal__ = lambda x, pars: multiplier*\
            self.nbaneal(x[0], pars[0], self.N_benef_anneal_g1)
        self.__nraneal__ = lambda x, pars: multiplier*\
            self.nraneal(x[0], pars[0], self.N_revers_anneal_g1)
        self.__nconst__ = lambda x, pars: multiplier*\
            self.nconst(x[0], pars[0], self.fluence)
        self.__nconst_log__ = lambda x, pars: multiplier*\
            self.nconst_log(x[0], pars[0], self.fluence)
        
        self.Vfunc_ndonor = rt.TF1("Vfunc_ndonor", self.__ndonor__, 1.49076e9,1.54862e9, 1)
        self.Vfunc_nbaneal = rt.TF1("Vfunc_nbaneal", self.__nbaneal__, 1.49076e9,1.54862e9, 1)
        self.Vfunc_nraneal = rt.TF1("Vfunc_nraneal", self.__nraneal__, 1.49076e9,1.54862e9, 1)
        self.Vfunc_nconst = rt.TF1("Vfunc_nconst", self.__nconst__, 1.49076e9,1.54862e9, 1)
        self.Vfunc_nconst_log = rt.TF1("Vfunc_nconst_log", self.__nconst_log__, 1.49076e9,1.54862e9, 1)

        canva1 = rt.TCanvas("Constituents","Neff",1600,800)
        canva1.Draw()
        self.Vfunc_ndonor.SetParameter(0,1.)
        self.Vfunc_nbaneal.SetParameter(0,self.Vfunc.GetParameter(0))
        self.Vfunc_nraneal.SetParameter(0,self.Vfunc.GetParameter(2))
        self.Vfunc_nconst.SetParameter(0,self.Vfunc.GetParameter(1))
        self.Vfunc_nconst_log.SetParameter(0,self.Vfunc_log.GetParameter(1))
        # print self.Vfunc.GetParameter(1)
        self.Vfunc_ndonor.SetLineColor(rt.kRed)
        self.Vfunc_nbaneal.SetLineColor(rt.kBlack)
        self.Vfunc_nraneal.SetLineColor(rt.kBlue)
        self.Vfunc_nconst.SetLineColor(rt.kGreen)
        self.Vfunc_nconst_log.SetLineColor(rt.kCyan)
        
        self.Vfunc_ndonor.GetYaxis().SetRangeUser(0,1e13)
        self.Vfunc_ndonor.Draw()
        self.Vfunc_nbaneal.Draw("same l")
        self.Vfunc_nraneal.Draw("same l")
        self.Vfunc_nconst.Draw("same l")
        self.Vfunc_nconst_log.Draw("same l")

        leg_pos = [0.135,0.5,.48,.8]
        leg = rt.TLegend(leg_pos[0],leg_pos[1],leg_pos[2],leg_pos[3])
        
        leg.AddEntry(self.Vfunc_ndonor,     "ndonor"    )
        leg.AddEntry(self.Vfunc_nbaneal,    "nbaneal"   )
        leg.AddEntry(self.Vfunc_nraneal,    "nraneal"   )
        leg.AddEntry(self.Vfunc_nconst,     "nconst"    )
        leg.AddEntry(self.Vfunc_nconst_log, "nconst_log")

        leg.Draw("same")

        canva1.SaveAs("constituents.pdf")

    def flog2(self, x, p0, p1, p2, N_donor, N_benef_anneal_g1, N_constdamage_g1, N_revers_anneal_g1):
        index = 0
        ret = abs(-N_donor[index].Eval(x) +
                p0*N_benef_anneal_g1[index].Eval(x) +
                p1*math.log(p2*fluence.Eval(x))*N_constdamage_g1[index].Eval(x) +
                p3*N_revers_anneal_g1[index].Eval(x))
        return ret
    
    def f(self, x, p0, p1, p2, N_donor, N_benef_anneal_g1, N_constdamage_g1, N_revers_anneal_g1):
        ret = abs(-N_donor.Eval(x) +
                  p0*N_benef_anneal_g1.Eval(x) +
                  p1*N_constdamage_g1.Eval(x) +
                  p2*N_revers_anneal_g1.Eval(x))
        return ret

    def falpha(self, x, p0, p1, p2, p3, corr_phi, N_donor, N_benef_anneal_g1, N_constdamage_g1, N_revers_anneal_g1):
        ret = corr_phi*p3*abs(-N_donor.Eval(x) +
                  p0*N_benef_anneal_g1.Eval(x) +
                  p1*N_constdamage_g1.Eval(x) +
                  p2*N_revers_anneal_g1.Eval(x))
        return ret
    
    def chi2_vs_alpha(self, func, data, start_date):
        start_alpha = 0.3
        end_alpha = 2.3
        diff_alpha = 0.05
        
        alpha_arr = list()
        alpha_i = start_alpha
        while alpha_i <= end_alpha: 
            alpha_arr.append(alpha_i)
            alpha_i += diff_alpha
        data_points = list()
        # print "alpha_arr = %s"%alpha_arr
        DataPoint = collections.namedtuple("DataPoints", 'alpha chi2')

        for i, alpha in enumerate(alpha_arr):
            func.SetParameter(3, alpha)
            chi2, ndof = calculate_chi2(data, func, start_date)
            # chi2ndof = chi2/float(ndof-1)
            chi2ndof = chi2
            dp = DataPoint(alpha, chi2ndof)
            data_points.append(dp)
        
        return data_points
    
    def fit_alpha_mods(self, modx, layer, dates_sorted, is_log, fileName, gC, P0,corr_phi):
        thickness = 285.
        multiplier = 1.6021766208e-13/(2*11.68*8.854187817)*thickness*thickness
        canva1 = rt.TCanvas("fit_vdepl_alpha_modx","fit_vdepl_alpha_modx",800,800)
        canva1.Draw()
        if not is_log:
            self.falphaf = lambda x, pars: multiplier * \
                self.falpha(x[0], pars[0], pars[1], pars[2], pars[3], pars[4], self.N_donor,
                            self.N_benef_anneal_g1, self.fluence, self.N_revers_anneal_g1)
        else:
            self.falphaf = lambda x, pars: multiplier * \
                self.falphalog_l1(x[0], pars[0], pars[1], pars[2], pars[3], pars[4], self.N_donor,
                            self.N_benef_anneal_g1, self.fluence, self.N_revers_anneal_g1)
        
        self.Vfunc_alpha = rt.TF1("VFuncAlpha", self.falphaf, 1.49076e9, 1.54862e9, 5)
        # self.Vfunc_alpha = rt.TF1("VFuncAlpha", self.falpha_log, 1.49076e9, 1.54862e9, 4)
        
        self.Vfunc_alpha.FixParameter(0,1.4e-2)
        self.Vfunc_alpha.FixParameter(1,gC)
        self.Vfunc_alpha.FixParameter(2,6.7e-2)
        self.Vfunc_alpha.SetParameter(3,1.)
        self.Vfunc_alpha.FixParameter(4,corr_phi)
        if not is_log: print("cor_phi = %s"%corr_phi)
        self.tg_dict_add[layer][modx].Fit("VFuncAlpha","QC","",
                dates_sorted[P0-1]-1,dates_sorted[-1]+1)
        
        canva1.SaveAs(fileName)

        return self.Vfunc_alpha.GetParameter(3)

    def h_chi2_vs_alpha(self, data_points):
        diff_alpha = data_points[1].alpha - data_points[0].alpha
        histo = rt.TH1D("chi2_vs_alpha", "chi2 vs alpha", len(data_points), 
                        data_points[0].alpha - diff_alpha/2., data_points[-1].alpha + diff_alpha/2.)
        for i,dp in enumerate(data_points):
            histo.SetBinContent(i, dp.chi2)
        return histo
        
    def draw_chi2_vs_alpha(self, histo, fileName):
        canva1 = rt.TCanvas("chi2","chi2",800,800)
        histo.SetMarkerStyle(20)
        histo.Draw("P")
        rt.gStyle.SetOptStat(0000000)
        canva1.SaveAs(fileName)

    def fit_gConstants(self, dir_out):
        thickness = 285.
        multiplier = 1.6021766208e-13/(2*11.68*8.854187817)*thickness*thickness
        if self.log_dependence:
            self.f1 = lambda x, pars: multiplier * \
                self.flog_l1(x[0], pars[0], pars[1], pars[2], self.N_donor,
                            self.N_benef_anneal_g1, self.fluence, self.N_revers_anneal_g1)
        else:
            self.f1 = lambda x, par: multiplier * \
                self.f(x[0], par[0], par[1], par[2], self.N_donor,
                        self.N_benef_anneal_g1, self.fluence, self.N_revers_anneal_g1)
        self.f2 = lambda x, pars: multiplier * \
                self.flog_l1(x[0], pars[0], pars[1], pars[2], self.N_donor,
                            self.N_benef_anneal_g1, self.fluence, self.N_revers_anneal_g1)
        
        print "Start fit gA,gC,gY!"
        
        
        self.falpha_lin = lambda x, pars: multiplier * \
            self.falpha(x[0], pars[0], pars[1], pars[2], pars[3], pars[4], self.N_donor,
                        self.N_benef_anneal_g1, self.fluence, self.N_revers_anneal_g1)
        
        self.falpha_log = lambda x, pars: multiplier * \
            self.falphalog_l1(x[0], pars[0], pars[1], pars[2], pars[3], pars[4], self.N_donor,
                        self.N_benef_anneal_g1, self.fluence, self.N_revers_anneal_g1)

        k = self.choose_key
        self.f1_models = lambda x, par: multiplier * \
            self.f(x[0], par[0], par[1], par[2], self.N_donor,
                         self.N_benef_anneal_g1, self.N_constdamage_g1, self.N_revers_anneal_g1)
        self.Vfunc = rt.TF1("VFuncFit", self.f1, 1.49076e9,
                            1.54862e9, 3)
        self.Vfunc_log = rt.TF1("VFuncFit_log", self.f2, 1.49076e9,
                            1.54862e9, 3)
        self.Vfunc_oldAtlas = rt.TF1("Vfunc_oldAtlas",self.f1_models,1.49076e9,1.54862e9,3)
        self.Vfunc_cb = rt.TF1("Vfunc_cb",self.f1_models,1.49076e9,1.54862e9,3)
        self.Vfunc_rd48 = rt.TF1("Vfunc_rd48",self.f1_models,1.49076e9,1.54862e9,3)
        self.Vfunc.SetParName(0,"gA")
        self.Vfunc.FixParameter(0,1.4e-2)
        # self.Vfunc.SetParameter(0,0.7e-2)
        # self.Vfunc.SetParLimits(0,0.5e-2,2.5e-2)
        self.Vfunc.SetParName(1,"gC")
        self.Vfunc.SetParameter(1,0.5e-2)
        self.Vfunc.SetParLimits(1,0.2e-2,1.5e-2)
        if self.log_dependence:
            self.Vfunc.SetParameter(1,1e8)
            self.Vfunc.SetParLimits(1,0.,1e11)
        self.Vfunc.SetParName(2,"gY")
        # self.Vfunc.SetParameter(2,6.0e-2)
        # self.Vfunc.SetParLimits(2,4.e-2,10.0e-2)
        self.Vfunc.FixParameter(2,6.7e-2)
        # self.Vfunc.FixParameter(2,7.e-2)
        # self.Vfunc.SetParName(3, "g")
        # self.Vfunc.SetParameter(3, 1.5)
        # self.Vfunc.SetParLimits(3, 1., 3.)
        self.Vfunc.SetNpx(10000)
        self.Vfunc.SetLineColor(rt.kMagenta)
        self.Vfunc.SetLineWidth(2)
        self.Vfunc_oldAtlas.SetLineWidth(1)
        self.Vfunc_cb.SetLineWidth(1)
        self.Vfunc_rd48.SetLineWidth(1)
        attempts=1
        P0=1
        nP = len(self.vdepl[k]["Date"])
        dates_sorted = []
        for i in range(0,nP): dates_sorted.append(self.vdepl[k]["Date"][i])
        dates_sorted=sorted(dates_sorted)
        # print dates_sorted
        if k=="FDISK2": P0=1
        fit_options=""
        # for i in range(0,attempts): self.tg_dict[k].Fit("VFuncFit","VNM","",self.vdepl[k]["Date"][P0],self.vdepl[k]["Date"][self.nPoints[k]-1])
        print "DATES SORTED: %s"%dates_sorted
        print "DATA: %s"%self.tg_dict[k]
        for i in range(0,attempts): self.tg_dict[k].Fit("VFuncFit",fit_options,"",
                dates_sorted[P0-1]-1,dates_sorted[-1]+1)
        self.Vfunc_log.SetParName(0,"gA")
        self.Vfunc_log.FixParameter(0, 1.4e-2)
        # self.Vfunc.SetParameter(0,0.7e-2)
        # self.Vfunc.SetParLimits(0,0.5e-2,2.5e-2)
        self.Vfunc_log.SetParName(1,"gC")
        # self.Vfunc.SetParameter(1,0.5e-2)
        # self.Vfunc.SetParLimits(1,0.2e-2,1.0e-2)
        self.Vfunc_log.SetParameter(1,1e8)
        self.Vfunc_log.SetParLimits(1,0.,1e11)
        self.Vfunc_log.SetParName(2,"gY")
        # self.Vfunc.SetParameter(2,6.0e-2)
        # self.Vfunc.SetParLimits(2,4.e-2,10.0e-2)
        self.Vfunc_log.FixParameter(2,6.7e-2)
        # self.Vfunc_log.FixParameter(2,7.e-2)
        # self.Vfunc.SetParName(3, "g")
        # self.Vfunc.SetParameter(3, 1.5)
        # self.Vfunc.SetParLimits(3, 1., 3.)
        self.Vfunc_log.SetNpx(10000)
        self.Vfunc_log.SetLineColor(rt.kBlue)
        self.Vfunc_log.SetLineWidth(2)
        attempts=1
        
        if self.log_dependence or self.is_all:
            for i in range(0, attempts):
                self.tg_dict[k].Fit("VFuncFit_log", fit_options, "", dates_sorted[P0]-1,dates_sorted[-1]+1)
        # self.Vfunc.Draw("same")
        # print("self.vdepl[k][\"Date\"][2] = %s"%(self.vdepl[k]["Date"][2]))
        d0 = datetime.strptime(self.begin_run,self.date_format) - datetime(1970,1,1,0,0,0)
        self.check_values(d0.total_seconds(), d0.total_seconds()+86400*60,86400, self.fluence,
                            self.Vfunc, self.Vfunc_log, self.N_donor)
        self.chi2_lin, self.ndof_lin = calculate_chi2(self.tg_dict[k], self.Vfunc,dates_sorted[P0-1]-1)
        self.chi2_log, self.ndof_log = calculate_chi2(self.tg_dict[k], self.Vfunc_log,dates_sorted[P0-1]-1)
        self.ndof_lin -= 1
        self.ndof_log -= 1
        # self.ndof=nP-P0-1
        
        print "ndof_lin = %s"%self.ndof_lin
        print "ndof_log = %s"%self.ndof_log
        print "Chi2(lin) = %s"%self.chi2_lin
        print "Chi2(log) = %s"%self.chi2_log
        print "Chi2/ndof(lin): %4.2f / %s = %s"%(self.chi2_lin,self.ndof_lin,self.chi2_lin/float(self.ndof_lin))
        print "Chi2/ndof(log): %4.2f / %s = %s"%(self.chi2_log,self.ndof_log,self.chi2_log/float(self.ndof_log))
        
        aver_phi = 0.0859955711871
        corr_phi = 0
        
        self.Vfunc_alpha = rt.TF1("VFuncAlpha", self.falpha_lin, 1.49076e9, 1.54862e9, 5)
        self.Vfunc_alpha.SetParameter(0, self.Vfunc.GetParameter(0))
        self.Vfunc_alpha.SetParameter(1, self.Vfunc.GetParameter(1))
        self.Vfunc_alpha.SetParameter(2, self.Vfunc.GetParameter(2))
        self.Vfunc_alpha.SetParameter(3, 0.3)
        self.Vfunc_alpha.FixParameter(4, corr_phi)

        self.Vfunc_alpha_log = rt.TF1("VFuncAlphaLog", self.falpha_log, 1.49076e9, 1.54862e9, 5)
        self.Vfunc_alpha_log.SetParameter(0, self.Vfunc_log.GetParameter(0))
        self.Vfunc_alpha_log.SetParameter(1, self.Vfunc_log.GetParameter(1))
        self.Vfunc_alpha_log.SetParameter(2, self.Vfunc_log.GetParameter(2))
        self.Vfunc_alpha_log.SetParameter(3, 0.3)
        self.Vfunc_alpha_log.FixParameter(4, corr_phi)
        
        # dir_out = "r_vs_z_2linesfit_aver_zabs_log_minibias"
        of_alpha = open("%s/output.txt"%dir_out, "w+")
        
        ldrs = [5,5,6]
        z_arr = list()
        z_arr_err = list()
        r_data_aversim_arr = list()
        z_arr_sim = [ -i/10. for i in range(0,350) ]
        my_Fl = GetFlunce("../PixelMonitoring/FLUKA/fluence_field.root")
        fl_sim_f = lambda x: my_Fl.fluence(x[0],x[1])
        zr_arr_sim = [[2.92,z] for i,z in enumerate(z_arr_sim)]
        z_arr_sim_err = [0.05 for i,z in enumerate(z_arr_sim)]
        z_arr_sim = array('d', z_arr_sim)
        # zr_arr_sim = array('d', r_arr_sim)
        z_arr_sim_err = array('d', z_arr_sim_err)
        r_data_sim_arr = list()
        r_data_sim_arr_err = list()
        r_data_aversim_arr_err = list()
        r_data_norm_arr = list()
        r_data_norm_arr_err = list()
        fl_aver_mods = {
            "mod1": 0.0875378424909, 
            "mod2": 0.0940240145989, 
            "mod3": 0.0807983035748, 
            "mod4": 0.0633744011094,
        }
        # fluence = [
        #     None,
        #     0.0940240145989,
        #     None,
        #     0.0633744011094,
        #     ]
        # dfluence = [
        #     None,
        #     0.00683477366132,
        #     None,
        #     0.00691847444178,
        # ]
        # fl_aver_mods = {"mod1": 0.08980608172725, 
        #                 "mod2": 0.08623534926899565, 
        #                 "mod3": 0.08238739196374782, 
        #                 "mod4": 0.08085061123838258}

        fl_aver_mods_err = {
            "mod1": 0.00136665214972,
            "mod2": 0.00683477366132,
            "mod3": 0.00153435147896,
            "mod4": 0.00691847444178,
        }
        z_arr_dict = {
            "mod1":3.35,
            "mod2":10.05,
            "mod3":16.75,
            "mod4":23.45
        }
        # fl_aver_mods = {"mod1": 0.08980608172725, 
        #                 "mod2": 0.08623534926899565, 
        #                 "mod3": 0.08238739196374782, 
        #                 "mod4": 0.08085061123838258}
        # z_arr = {
        #             "mod1":3.35,
        #             "mod2":10.05,
        #             "mod3":16.75,
        #             "mod4":23.45}
        
        for i,modx in enumerate(["mod2", "mod4"]):
            # z_arr.append(abs(fl_pos_dict["BPix_BmO_SEC7_LYR1_LDR%sF_MOD%s"%(ldrs[i],modx[-1])].z))
            # z_arr.append(z_arr_arr[modx])
            z_arr_err.append(6.4/2.)
            # corr_phi = fl_pos_dict["BPix_BmO_SEC2_LYR1_LDR2F_MOD3"].fluence/aver_phi
            # corr_phi = fl_pos_dict["BPix_BmO_SEC2_LYR1_LDR2F_MOD4"].fluence/aver_phi
            corr_phi = fl_aver_mods[modx]/aver_phi
            self.Vfunc_alpha.FixParameter(4, corr_phi)
            self.Vfunc_alpha_log.FixParameter(4, corr_phi)
            my_data_points = self.chi2_vs_alpha(self.Vfunc_alpha, self.tg_dict_add[k][modx], dates_sorted[P0-1]-1)
            my_histo = self.h_chi2_vs_alpha(my_data_points)
            self.draw_chi2_vs_alpha(my_histo, "chi2_vs_alpha_%s.pdf"%modx)
            my_data_points_log = self.chi2_vs_alpha(self.Vfunc_alpha_log, self.tg_dict_add[k][modx], dates_sorted[P0-1]-1)
            my_histo_log = self.h_chi2_vs_alpha(my_data_points_log)
            self.draw_chi2_vs_alpha(my_histo_log, "chi2_vs_alpha_log_%s.pdf"%modx)
            
            # corr_phi = fl_pos_dict["BPix_BmO_SEC7_LYR1_LDR%sF_MOD%s"%(ldrs[i],modx[-1])].fluence/aver_phi
            corr_phi = fl_aver_mods[modx]/aver_phi
            corr_phi_def = corr_phi
            alpha = self.fit_alpha_mods(modx,k,dates_sorted,is_log = True, fileName = "Vdepl_alpha/%s_alpha_log.pdf"%modx, gC = self.Vfunc_log.GetParameter(1), P0=P0, corr_phi=corr_phi)
            # alpha = self.fit_alpha_mods(modx,k,dates_sorted,is_log = True, fileName = "Vdepl_alpha/%s_alpha_log.pdf"%modx, gC = 86.2e9, P0=P0, corr_phi=corr_phi)
            alpha_real_log = alpha
            # position_modx = fl_pos_dict["BPix_BmO_SEC7_LYR1_LDR%sF_MOD%s"%(ldrs[i],modx[-1])]
            of_alpha.write("%s (%s, %s, %s)\n"%(modx, 2.93, z_arr_dict[modx], fl_aver_mods[modx]))
            # of_alpha.write("real\n")
            # of_alpha.write("log: ")
            # of_alpha.write("alpha = %s\n"%alpha)
                
            corr_phi = 1.
            alpha = self.fit_alpha_mods(modx,k,dates_sorted,is_log = True, fileName = "Vdepl_alpha/%s_alpha_log.pdf"%modx, gC = self.Vfunc_log.GetParameter(1), P0=P0, corr_phi=corr_phi)
            alpha_aver_log = alpha

            corr_phi = fl_aver_mods[modx]/aver_phi
            alpha = self.fit_alpha_mods(modx,k,dates_sorted,is_log = not True, fileName="Vdepl_alpha/%s_alpha_lin.pdf"%modx, gC = self.Vfunc.GetParameter(1), P0=P0, corr_phi=corr_phi)
            # alpha = self.fit_alpha_mods(modx,k,dates_sorted,is_log = not True, fileName="Vdepl_alpha/%s_alpha_lin.pdf"%modx, gC = 0.66e-2, P0=P0, corr_phi=corr_phi_def)
            alpha_real_lin = alpha
            evdepl = self.vdepl_add[k][modx]["eVdepl"][0]
            vdepl = self.vdepl_add[k][modx]["Vdepl"][0]
            dalpha = alpha*math.sqrt((evdepl/vdepl)*
                                     (evdepl/vdepl)+
                                     (fl_aver_mods_err[modx]/fl_aver_mods[modx])*
                                     fl_aver_mods_err[modx]/fl_aver_mods[modx])
            of_alpha.write("lin: ")
            of_alpha.write("alpha = %s +- %s\n"%(alpha, dalpha))
            
            # of_alpha.write("aver\n")
            # of_alpha.write("log: ")
            # of_alpha.write("alpha = %s\n"%alpha_aver_log)
            alpha = self.fit_alpha_mods(modx,k,dates_sorted,is_log = not True, fileName="Vdepl_alpha/%s_alpha_lin.pdf"%modx, gC = self.Vfunc.GetParameter(1), P0=P0, corr_phi=1.)
            alpha_aver_lin = alpha
            # of_alpha.write("lin: ")
            # of_alpha.write("alpha = %s\n"%alpha)

            r_data_sim_arr.append(alpha_real_log)
            # r_data_sim_arr.append(alpha_real_lin)
            # r_data_sim_arr.append(position_modx.fluence)
            r_data_sim_arr_err.append(0)
            r_data_aversim_arr.append(alpha_aver_log)
            # r_data_aversim_arr.append(alpha_aver_lin)
            r_data_aversim_arr_err.append(0)
            r_data_norm_arr.append(r_data_sim_arr[i]/r_data_sim_arr[0])
            r_data_norm_arr_err.append(0)

        z_arr = array('d', [3.35, 10.05, 16.75, 23.45])
        z_arr_err = array('d', z_arr_err)
        r_data_aversim_arr = array('d', r_data_aversim_arr)
        r_data_aversim_arr_err = array('d', r_data_aversim_arr_err)
        r_data_sim_arr = array('d', r_data_sim_arr)
        r_data_sim_arr_err = array('d', r_data_sim_arr_err)
        r_data_norm_arr = array('d', r_data_norm_arr)
        r_data_norm_arr_err = array('d', r_data_norm_arr_err)

        r_data_aversim_gr       = rt.TGraphErrors(3, z_arr, r_data_aversim_arr, z_arr_err, r_data_aversim_arr_err)
        # r_data_aversim_gr_atlas = rt.TGraphErrors("r_data_aversim_atlas", "r_data_aversim_atlas", 3, z_arr_atlas, r_data_aversim_arr_atlas, z_arr_err_atlas, r_data_aversim_arr_err_atlas)
        r_data_sim_gr           = rt.TGraphErrors(3, z_arr, r_data_sim_arr, z_arr_err, r_data_sim_arr_err)
        # r_data_sim_gr_atlas     = rt.TGraphErrors("r_data_aver_atlas", "r_data_aver_atlas", 3, z_arr_atlas, r_data_sim_arr_atlas, z_arr_err_atlas, r_data_sim_arr_err_atlas)
        r_data_norm             = rt.TGraphErrors(3, z_arr, r_data_norm_arr, z_arr_err, r_data_norm_arr_err)
        # r_data_norm_atlas       = rt.TGraphErrors("r_data_aver_atlas", "r_data_aver_atlas", 3, z_arr_atlas, r_data_aver_arr_atlas, z_arr_err_atlas, r_data_aver_arr_err_atlas)
        
        canva2 = rt.TCanvas("c2","c2",0,0,800,800)
        canva2.Draw()
        r_data_aversim_gr.SetTitle("r vs z")
        r_data_aversim_gr.GetYaxis().SetTitle("r value")
        r_data_aversim_gr.GetXaxis().SetTitle("|z| [cm]")
        r_data_aversim_gr.SetMarkerColor(rt.kBlack)
        r_data_aversim_gr.SetMarkerStyle(20)
        r_data_aversim_gr.Draw("ap")
        canva2.SaveAs("%s/aver.pdf"%dir_out)
        
        canva3 = rt.TCanvas("c2","c2",0,0,800,800)
        canva3.Draw()
        r_data_sim_gr.SetTitle("r vs z")
        r_data_sim_gr.GetYaxis().SetTitle("r value")
        r_data_sim_gr.GetXaxis().SetTitle("|z| [cm]")
        r_data_sim_gr.SetMarkerColor(rt.kBlack)
        r_data_sim_gr.SetMarkerStyle(20)
        r_data_sim_gr.Draw("ap")
        canva3.SaveAs("%s/sim.pdf"%dir_out)
        
        canva4 = rt.TCanvas("c2","c2",0,0,800,800)
        canva4.Draw()
        r_data_norm.SetTitle("r vs z")
        r_data_norm.GetYaxis().SetTitle("r value")
        r_data_norm.GetXaxis().SetTitle("|z| [cm]")
        r_data_norm.SetMarkerColor(rt.kBlack)
        r_data_norm.SetMarkerStyle(20)
        r_data_norm.Draw("ap")
        canva4.SaveAs("%s/norm.pdf"%dir_out)
       
       # r_sim_aversim_gr        = rt.TGraph("r_sim_aversim", "r_sim_aversim", 3, z_arr, r_data_aversim_arr)
        # r_sim_aversim_gr_atlas  = rt.TGraph("r_sim_aversim_atlas", "r_sim_aversim_atlas", 3, z_arr_atlas, r_data_aversim_arr_atlas)
        # r_sim_sim_gr            = rt.TGraph("r_sim_aver", "r_sim_aver", 3, z_arr, r_data_sim_arr)
        # r_sim_sim_gr_atlas      = rt.TGraph("r_sim_aver_atlas", "r_sim_aver_atlas", 3, z_arr_atlas, r_data_sim_arr_atlas)
        # r_sim_norm              = rt.TGraph("r_sim_aver", "r_sim_aver", 3, z_arr, r_data_aver_arr)
        # r_sim_norm_atlas        = rt.TGraph("r_sim_aver_atlas", "r_sim_aver_atlas", 3, z_arr_atlas, r_data_aver_arr_atlas)

        of_alpha.close()
        
        self.Vfunc_oldAtlas.SetParameter(0,0.7e-2)
        self.Vfunc_oldAtlas.SetParameter(2,6.0e-2)
        self.Vfunc_oldAtlas.SetParameter(1,0.7e-2)
        # self.Vfunc_oldAtlas.SetParameter(3,0.7e-10)
        # self.Vfunc_oldAtlas.SetParameter(0,0)
        # self.Vfunc_oldAtlas.SetParameter(2,6.0e-2)
        # self.Vfunc_oldAtlas.SetParameter(1,0)
        # self.Vfunc_oldAtlas.SetLineColor(rt.kBlack)
        # self.Vfunc_oldAtlas.SetNpx(10000)
        # self.Vfunc_oldAtlas.Draw("same")
        self.gC_oxy = 0
        if self.ring == "1":
            self.gC_oxy_cb = 1.39e-2
            self.gC_oxy_rd = 1.03e-2
            self.gY_oxy_cb = 6.49e-2
            self.gY_oxy_rd = 5.12e-2
        elif self.ring == "11":
            #layer 1
            self.gC_oxy_cb = 0.78e-2
            self.gC_oxy_rd = 0.88e-2
            self.gY_oxy_cb = 5.0e-2
            self.gY_oxy_rd = 6.7e-2
        else:
            self.gC_oxy_cb = 1.03e-2
            self.gC_oxy_rd = 1.03e-2
            self.gY_oxy_cb = 1.03e-2
            self.gY_oxy_rd = 1.03e-2
        self.gC_oxy_cb_str = "%s#upoint10^{-2}" % round(self.gC_oxy_cb*1e2, 2)
        self.gC_oxy_rd_str = "%s#upoint10^{-2}" % round(self.gC_oxy_rd*1e2, 2)
        self.Vfunc_cb.SetParameter(0,1.4e-2)
        self.Vfunc_cb.SetParameter(2, self.gY_oxy_cb)
        self.Vfunc_cb.SetParameter(1, self.gC_oxy_cb)
        # self.Vfunc_cb.SetParameter(0, 0.7e-2)
        # self.Vfunc_cb.SetParameter(2, 0)
        # self.Vfunc_cb.SetParameter(1, 0)
        self.Vfunc_cb.SetLineColor(rt.kBlue)
        self.Vfunc_cb.SetNpx(10000)
        self.Vfunc_cb.SetLineWidth(2)
        # self.Vfunc_cb.Draw("same")
        self.Vfunc_rd48.SetParameter(0,1.4e-2)
        self.Vfunc_rd48.SetParameter(2, self.gY_oxy_rd)
        self.Vfunc_rd48.SetParameter(1, self.gC_oxy_rd)
        # self.Vfunc_rd48.SetParameter(0, 0)
        # self.Vfunc_rd48.SetParameter(2, 0)
        # self.Vfunc_rd48.SetParameter(1, 0.7e-2)
        self.Vfunc_rd48.SetLineColor(rt.kGreen+1)
        self.Vfunc_rd48.SetNpx(10000)
        # self.Vfunc_rd48.Draw("same")
        self.Vfunc_rd48.SetLineWidth(2)
        self.gA = self.Vfunc.GetParameter(0)
        self.gC = self.Vfunc.GetParameter(1)
        self.gC_logarithmic = self.Vfunc_log.GetParameter(1)
        self.gY = self.Vfunc.GetParameter(2)
        # self.isfit = True
        # xmin = self.N_benef_anneal_g1[0].GetX()[0]
        # xmax = self.N_benef_anneal_g1[0].GetX()[self.N_benef_anneal_g1[0].GetN()-1]
        xmin = self.N_benef_anneal_g1.GetX()[0]
        xmax = self.N_benef_anneal_g1.GetX()[self.N_benef_anneal_g1.GetN()-1]
        self.Vfunc_ext = rt.TF1("VFuncFit_ext", self.f1, xmin,
                                xmax, 3)
        self.Vfunc_ext.SetNpx(10000)
        self.Vfunc_ext.SetParameter(0, self.gA)
        self.Vfunc_ext.SetParameter(2, self.gY)
        self.Vfunc_ext.SetParameter(1, self.gC)
        # self.Vfunc.Draw("same")
        print "Fit is finished: gA = %s; gC = %s; gY = %s"%(self.gA,self.gC,self.gY)
        self.params = rt.TH1F("params", "", 3, 0, 3)
        self.params.GetXaxis().SetBinLabel(1,"gA")
        self.params.GetXaxis().SetBinLabel(2,"gC")
        self.params.GetXaxis().SetBinLabel(3,"gY")
        self.params.Fill(0.5, self.gA)
        self.params.Fill(1.5, self.gC)
        self.params.Fill(2.5, self.gY)

    def draw_graphs(self,canv,vdepl,vdepl_sim):
        canv.Draw()
        vdepl_sim.Draw("AL")
        vdepl.Draw("P")

    def draw_data(self,vdepl,opt=""):
        vdepl.SetMarkerColor(rt.kBlack)
        vdepl.Draw(opt)
    
    def check_values(self,xmin,xmax,dx,phi_eq,flin,flog,ndonor):
        fout = open("check_vdepl_vals.txt", "w+")
        x = xmin
        d0 = datetime(1970,1,1,0,0,0)
        fout.write("     Date and time     |        Phieq        |      ln(Phieq)      |        vlin         |        vlog         |        ndonor\n")
        fout.write("-----------------------+---------------------+---------------------+---------------------+---------------------+---------------------|\n")
        while x < xmax:
            date = d0 + timedelta(seconds = x)
            fout.write("%s    |%20.13e |%20.13e |%20.13e |%20.13e |%20.13e\n"%(date.strftime(self.date_format),
            phi_eq.Eval(x),86.2e9*math.log(phi_eq.Eval(x)),flin.Eval(x),flog.Eval(x),ndonor.Eval(x)))
            x += dx

    
    def add_legend(self, leg_pos, disk, ring):
        leg=rt.TLegend(leg_pos[0],leg_pos[1],leg_pos[2],leg_pos[3])
        for k in self.keys:
            leg.AddEntry(self.tg_dict[k], "Data: From Cluster Charge")
            if not self.isfit:
                leg.AddEntry(self.tg_dict_sim[k],"V_{dep} simulation, %s"%k)
        self.legend = leg
        self.legend.SetFillStyle(0)
        self.legend.SetBorderSize(0)
        # self.legend.SetLineStyle(0)

    def save_vdepl(self,outputFile_name,dir_out):
        outputFile = rt.TFile.Open(outputFile_name, "RECREATE")
        # future
        # for value in variable:
            # pass
        # now
        # k = "FDISK1"
        self.choose_key = "LAYER1"
        k = self.choose_key
        self.canv1 = rt.TCanvas("Vdepl_%s"%k.replace(" ", ""),"Vdepl_%s"%k.replace(" ", ""),1600,800)
        # if not self.isfit:
        # self.draw_graphs(self.canv1,self.tg_dict[k],self.tg_dict_sim[k])
        if self.isfit: self.fit_gConstants(dir_out)
        self.canv1.Draw()
        h0 = rt.TH1D("h0def", "h0def", 10000, 1.49076e9,1543762080.0)
        # h0 = self.Vfunc.GetHistogram()
        h0.SetTitle("Phase-1 Barrel Pixel - Full depletion voltage vs day")
        h0.GetYaxis().SetTitle("V_{depl} [V]")
        h0.GetXaxis().SetTitle("Date")
        h0.GetYaxis().SetRangeUser(0,1000)
        h0.GetXaxis().SetNdivisions(509)
        h0.GetXaxis().SetTimeFormat("%d/%m/%Y%F1970-01-01 00:00:00s0 GMT")
        if self.plot_add:
            xmin_ = min(self.vdepl_add[k]['mod2']["Date"]) - 3600*75
            xmax_ = max(self.vdepl_add[k]['mod2']["Date"]) + 3600*75
        h0.GetXaxis().SetRangeUser(xmin_, xmax_)
        h0.GetYaxis().SetRangeUser(100,500)
        rt.gStyle.SetOptStat(0000)
        h0.Draw()
        # directory = "../PixelMonitoring/BPix_studies"
        # self.canv1.SaveAs("%s/%s.pdf" %(directory,outputFile_name.split('.')[0]))
        # exit()
        self.draw_data(self.tg_dict[k], opt="p same")
        # print self.tg_dict[k].GetXaxis().GetXmax()
        hist_vfunc_draw_lin = self.Vfunc.GetHistogram()
        first_bin = hist_vfunc_draw_lin.FindBin((datetime.strptime(self.begin_run,self.date_format)-datetime(1970,1,1,0,0,0)).total_seconds())
        for i in range(1,first_bin+1):
            hist_vfunc_draw_lin.SetBinContent(i,110.)
        hist_vfunc_draw_lin.Draw("same")
        
        hist_vfunc_draw_log = self.Vfunc_log.GetHistogram()
        first_bin_log = hist_vfunc_draw_log.FindBin((datetime.strptime(self.begin_run,self.date_format)-datetime(1970,1,1,0,0,0)).total_seconds())
        # for i in range(1,first_bin_log+1):
        #     hist_vfunc_draw_log.SetBinContent(i,110.)
        hist_vfunc_draw_log.Draw("same")
        
        # self.Vfunc.Draw("same")
        gC = self.gC
        gC_str = ""
        
        cmsText = "CMS"
        cmsTextFont = 61
        writeExtraText = True
        extraText = "Preliminary"
        extraText2 = "2020"
        extraText3 = "CMS FLUKA study v3.23.1.0"
        extraTextFont = 52
        lumiTextSize = 0.5
        lumiTextOffset = 0.15
        cmsTextSize = 0.5
        cmsTextOffset = 0.1
        # only used in outOfFrame version
        relPosX = 0.045
        relPosY = 0.035
        relExtraDY = 1.2
        # ratio of "CMS" and extra text size
        extraOverCmsTextSize = 0.65
        # lumi_13TeV = "20.1 fb^{-1}"
        # lumi_8TeV = "19.7 fb^{-1}"
        # lumi_7TeV = "5.1 fb^{-1}"
        # lumiText
        # lumiText += lumi_8TeV
        # lumiText += " (13 TeV)"
        # lumiText = "#sqrt{s} = 13 TeV "
        t = self.canv1.GetTopMargin()
        b = self.canv1.GetBottomMargin()
        r = self.canv1.GetRightMargin()
        l = self.canv1.GetLeftMargin()
        latex = rt.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(rt.kBlack)
        extraTextSize = extraOverCmsTextSize*cmsTextSize*1.1
        latex.SetTextFont(42)
        latex.SetTextAlign(31)
        latex.SetTextSize(lumiTextSize*t)
        # latex.DrawLatex(1-r, 1-t+lumiTextOffset*t, lumiText)
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11)
        latex.SetTextSize(cmsTextSize*t)
        latex.DrawLatex(l+0.05, 1-t+lumiTextOffset*t-0.09, cmsText)
        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize*t)
        # latex.DrawLatex(l+0.05, 1-t+lumiTextOffset*t-0.09-0.06, extraText)
        latex.DrawLatex(l+0.05+0.06, 1-t+lumiTextOffset*t-0.09, extraText)
        latex.SetTextFont(extraTextFont-10)
        # latex.DrawLatex(l+0.14, 1-t+lumiTextOffset*t-0.09-0.06, extraText2)
        # latex.DrawLatex(l+0.14+0.06, 1-t+lumiTextOffset*t-0.09, extraText2)
        latex.SetTextFont(51)
        latex.SetTextSize(extraTextSize*t*0.85)
        latex.DrawLatex(l+0.6, 1-t+lumiTextOffset*t-0.05, extraText3)
        # fpix_logo = "Forward Pixel Ring %s Disk %s" % (self.ring, self.disk)
        fpix_logo = "Barrel Pixel Layer 1"
        TrackSelctionText = fpix_logo
        # latex.SetTextFont(61)
        latex.SetTextSize(extraTextSize*t)
        latex.SetTextFont(extraTextFont+10)
        latex.DrawLatex(l+0.4, 1-t+lumiTextOffset*t-0.15, TrackSelctionText)
        TrackSelctionText2 = ""
        if self.disk == "1":
            TrackSelctionText2 = "z = 32 cm"
        elif self.disk == "2":
            TrackSelctionText2 = "z = 40 cm"
        elif self.disk == "3":
            TrackSelctionText2 = "z = 50 cm"
        # TrackSelctionText2 = "z = 0 cm"
        TrackSelctionText2 = ""
        latex.SetTextFont(61)
        latex.SetTextSize(extraTextSize*t)
        latex.DrawLatex(l+0.4, 1-t+lumiTextOffset*t -
                        0.15-0.04, TrackSelctionText2)
        
        chi2_lin_text = "#chi_{lin}^{2} / ndof = %4.1f / %s = %4.1f"%(self.chi2_lin,self.ndof_lin,self.chi2_lin/float(self.ndof_lin))
        chi2_log_text = "#chi_{log}^{2} / ndof = %4.1f / %s = %4.1f"%(self.chi2_log,self.ndof_log,self.chi2_log/float(self.ndof_log))

        latex.DrawLatex(l+0.4, 1-t+lumiTextOffset*t -
                        0.15-0.04, TrackSelctionText2)
        
        # latex.SetTextSize(51)
        latex.SetTextSize((extraTextSize-0.05)*t)
        latex.SetTextFont(extraTextFont-10)
        latex.DrawLatex(0.4, 0.65, chi2_lin_text)
        latex.DrawLatex(0.4, 0.55, chi2_log_text)

        if not self.log_dependence:
            gC_str = "%s#upoint10^{-2}" % round(gC*1e2, 2)
        else:
            gC_str = "%s#upoint10^{9}" % round(gC/1e9, 1)
        model = "linear"
        model_fit = model
        unit = "cm^{-1}"
        if self.log_dependence:
            model_fit = "logarithm"
            unit = "cm^{-3}"
        self.legend.AddEntry(
            self.Vfunc, "Fit, %s, g_{C} = %s %s" % (model_fit, gC_str, unit))
        if self.is_all:
            self.legend.AddEntry(self.Vfunc_log,"Fit, %s, g_{C}^{log} = %s cm^{-3}"%("logarithm","%s#upoint10^{9}" % round(self.gC_logarithmic/1e9, 1)))
            # self.Vfunc_log.Draw("same")    
        # self.Vfunc_oldAtlas.Draw("same")
        # self.legend.AddEntry(self.Vfunc_oldAtlas,"Atlas 2017")
        if not self.is_all:
            self.Vfunc_cb.Draw("same")
            self.legend.AddEntry(self.Vfunc_cb,"CB-oxy, %s, g_{C} = %s cm^{-1}"%(model, self.gC_oxy_cb_str))
        if not self.is_all:
            self.Vfunc_rd48.Draw("same")
            self.legend.AddEntry(self.Vfunc_rd48, "RD48-oxy, %s, g_{C} = %s cm^{-1}" % (model, self.gC_oxy_rd_str))
        # self.Vfunc.Draw("same")
        
        self.tg_dict[k].Draw("same p")
        
        # self.tg_dict_add[k]["mod1"].SetMarkerStyle(3)
        # self.tg_dict_add[k]["mod1"].SetMarkerSize(1.7)
        # self.tg_dict_add[k]["mod1"].SetMarkerColor(rt.kRed)
        # self.tg_dict_add[k]["mod1"].SetLineColor(rt.kRed)
        # self.tg_dict_add[k]["mod1"].SetFillColor(0)
        # self.legend.AddEntry(self.tg_dict_add[k]["mod1"],"Module1")
        # self.tg_dict_add[k]["mod1"].Draw("same p")
        
        self.tg_dict_add[k]["mod2"].SetMarkerStyle(4)
        self.tg_dict_add[k]["mod2"].SetMarkerSize(1.7)
        self.tg_dict_add[k]["mod2"].SetMarkerColor(rt.kRed)
        self.tg_dict_add[k]["mod2"].SetLineColor(rt.kRed)
        self.tg_dict_add[k]["mod2"].SetFillColor(0)
        self.legend.AddEntry(self.tg_dict_add[k]["mod2"],"Module2")
        self.tg_dict_add[k]["mod2"].Draw("same p")

        # self.tg_dict_add[k]["mod3"].SetMarkerStyle(5)
        # self.tg_dict_add[k]["mod3"].SetMarkerSize(1.7)
        # self.tg_dict_add[k]["mod3"].SetMarkerColor(rt.kRed)
        # self.tg_dict_add[k]["mod3"].SetLineColor(rt.kRed)
        # self.tg_dict_add[k]["mod3"].SetFillColor(0)
        # self.legend.AddEntry(self.tg_dict_add[k]["mod3"],"Module3")
        # self.tg_dict_add[k]["mod3"].Draw("same p")

        self.tg_dict_add[k]["mod4"].SetMarkerStyle(7)
        self.tg_dict_add[k]["mod4"].SetMarkerSize(1.7)
        self.tg_dict_add[k]["mod4"].SetMarkerColor(rt.kRed)
        self.tg_dict_add[k]["mod4"].SetLineColor(rt.kRed)
        self.tg_dict_add[k]["mod4"].SetFillColor(0)
        self.legend.AddEntry(self.tg_dict_add[k]["mod4"],"Module4")
        self.tg_dict_add[k]["mod4"].Draw("same p")

            # raw_input("...")
        # else:
        #     self.draw_data(self.canv1,self.tg_dict[k])
        #     self.fit_gConstants()
        self.legend.Draw("same")
        self.canv1.Write()
        directory = "../PixelMonitoring/BPix_studies"
        self.canv1.SaveAs("%s/%s.pdf" %(directory,outputFile_name.split('.')[0]))
        self.canv1.SaveAs("%s/%s.png" %(directory,outputFile_name.split('.')[0]))
        # self.params.Write()
        # self.N_constdamage_g1[0].Write()
        # self.canv2 = rt.TCanvas("Vdepl_%s_sim" % k.replace(" ", ""), "Vdepl_%s_sim" % k.replace(" ", ""), 1600, 800)
        # self.canv2.Draw()
        # self.tg_dict[k].GetXaxis().SetTimeFormat(
        #     "%d/%m/%Y%F1970-01-01 00:00:00s0 GMT")
        # xmin = self.N_benef_anneal_g1[0].GetX()[0]
        # xmax = self.N_benef_anneal_g1[0].GetX()[self.N_benef_anneal_g1[0].GetN()-1]
        # xmin = self.N_benef_anneal_g1.GetX()[0]
        # print "xmin = %s"%xmin
        # xmin = 1.49077e9
        # xmax = self.N_benef_anneal_g1.GetX()[self.N_benef_anneal_g1.GetN()-1]
        # self.tg_dict[k].GetXaxis().SetRangeUser(xmin, xmax)
        # self.Vfunc_ext.GetXaxis().SetTimeFormat("%d/%m/%Y%F1970-01-01 00:00:00s0 GMT")
        # h0 = rt.TH1F("h0","h0",10000,xmin,xmax)
        # h0=self.Vfunc_ext.GetHistogram()
        # h0.GetXaxis().SetRangeUser(xmin,xmax)
        # h0.GetXaxis().SetNdivisions(505)
        # h0.Draw()
        # h0.GetXaxis().SetTimeFormat("%d/%m/%Y%F1970-01-01 00:00:00s0 GMT")
        # self.Vfunc_ext.SetLineColor(rt.kBlack)
        # self.Vfunc_ext.SetLineWidth(2)
        # self.Vfunc_ext.Draw("same l")
        
        # self.draw_data(self.canv1, self.tg_dict[k],"same p")
        # self.Vfunc_ext.GetXaxis().SetTimeFormat('%d/%m/%Y')
        # self.Vfunc_ext.GetXaxis().SetTimeOffset(0, "gmt")
        # self.canv2.Write()
        # self.canv2.Print('vdepl_pdf/'+outputFile_name.split('.')[0]+'.pdf')
        # self.tg_dict[k].Write()
        # h0.Write()
        outputFile.Close()


def main(options, args):
    BEGIN_RUN = "2017-05-23 14:32:22"
    # fileName="Vdepl_data/CMS_PIXELS_VDEPL.csv"
    # outputSimFile="simulation_results_gA_1p4.root"
    # chooseKeys=["FDISK1"]
    outFile       = options.output
    fileName      = options.data
    outputSimFile = options.input
    chooseKeys    = options.layerTag
    # print chooseKeys
    mykeys = {"l1": "LAYER1"}
    leg_pos = [0.135,0.5,.48,.8]
    is_log = options.log
    is_all = options.all
    plot_add = options.plot_add
    if is_all:
        leg_pos = [0.135,0.5,.37,.8]
    myVdepl = Vdepl(BEGIN_RUN, options.disk, options.ring, is_log, is_all, plot_add)
    myVdepl.vdepl_data_init()
    myVdepl.vdepl_add_data_init(["mod2","mod4"])
    # myVdepl.initDataStructFromCSV(fileName, mykeys=chooseKeys)
    # myVdepl.initDataStructFromCSV(fileName)
    # myVdepl.get_dates_vdepl()
    date_vdepl_dict = { 
        # '2017-05-23 14:32:22':  [110.00, 0.0],
        '2017-08-14 12:00:00':	[279.8770436083903, 2.861688650771493],
        '2017-09-23 12:00:00':	[287.2182156683212, 50.45926421614772],
        '2017-10-04 12:00:00':	[427.9244731211447, 86.39914685054431],
        '2017-10-27 12:00:00':	[412.2285242994644, 10.039207643040575],
        '2017-11-09 12:00:00':	[444.6075593345522, 13.329370633395312],
        
        '2018-05-05 12:00:00':	[219.70760583947464, 21.381952721439305],
        '2018-05-12 12:00:00':	[290.49303204464945, 18.013619051200934],
        '2018-05-24 12:00:00':	[332.8671113232275, 25.670488513703003],
        '2018-07-11 12:00:00':	[328.2472067937181, 43.353323844328884],
        '2018-07-30 12:00:00':	[306.18106342603124, 13.755913498527221],
        '2018-08-17 12:00:00':	[400.62356587524613, 23.104338842885902],
        '2018-09-01 12:00:00':	[425.3308148998366, 35.029054223180644],
        '2018-09-07 12:00:00':	[482.0076955636585, 63.93361480649718],
        '2018-09-26 12:00:00':	[444.56970431780763, 30.068930894297235],
        '2018-10-20 12:00:00':	[496.04567857680826, 107.6843645649626],

        # '2018-08-17 12:00:00':	[400.517398396, 0.0],
        # '2018-09-01 12:00:00':	[418.740359493, 0.0],
        # '2018-09-07 12:00:00':	[475.481560204, 0.0],
        # '2018-09-26 12:00:00':	[450.967069295, 0.0],
        # '2018-10-20 12:00:00':	[484.932440031, 0.0],
    }

    myVdepl.get_vdepl_from_fit(date_vdepl_dict, tag=mykeys[chooseKeys])
    vdepl_vals = [260.89, 249.69, 226.32, 210.66]
    vdepl_input_file = dict()
    mod_names = ["mod2", "mod4"]
    parent_input_dir = "../HVBiasScans/1f_mods_mini_scan/"
    input_data = dict()
    my_dates = [
        "2018-08-17 12:00:00",
        "2018-09-01 12:00:00",
        "2018-09-07 12:00:00",
        "2018-09-26 12:00:00",
        "2018-10-20 12:00:00",
    ]
    for my_date in my_dates:
        input_data = dict()
        dir_out = "r_vs_z_minibias/onelinefit_%s"%(my_date.split(' ')[0].replace('-','_'))
        for modx in mod_names:
            input_data[modx] = dict()
            vdepl_input_file[modx] = open("%s/pdfs_%s/vdepl_values_%s.txt"%(parent_input_dir,modx,modx))
            input_lines = vdepl_input_file[modx].readlines()
            for i,l in enumerate(input_lines):
                if i == 0:
                    continue
                else:
                    l = l.split('\t')
                    if my_date in l[0]:
                        input_data[modx][l[0]] = [float(l[1]),float(l[2])]
    
        print(str(input_data))
        myVdepl.get_vdepl_additional(input_data["mod2"], tag=mykeys[chooseKeys], modx = "mod2")
        myVdepl.get_vdepl_additional(input_data["mod4"], tag=mykeys[chooseKeys], modx = "mod4")
        myVdepl.init_vdepl_sim(outputSimFile)
        myVdepl.init_tgraphs()
        myVdepl.get_Vdepl_terms()
        myVdepl.apply_draw_options()
        myVdepl.add_legend(leg_pos, options.disk,options.ring)
        myVdepl.save_vdepl(options.output, dir_out)
        myVdepl.plot_constituents()


usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-i", "--input-file", action="store", type="string", dest="input")
parser.add_option("-d", "--data-file", action="store", type="string", dest="data")
parser.add_option("-b", "--batch", 	action="store_true", dest="batch")
parser.add_option("-o", "--output-file", 	action="store", type="string", dest="output")
parser.add_option("-l", "--layer-tag", 	action="store", type="string", dest="layerTag")
parser.add_option("", "--log-dep", 	action="store_true", dest="log")
parser.add_option("", "--all-dep", 	action="store_true", dest="all")
parser.add_option("", "--plot-add", action="store_true", dest="plot_add")
parser.add_option("", "--disk", action="store", dest="disk")
parser.add_option("", "--ring", action="store", dest="ring")
# parser.add_option("--histo", action="store", type="string", dest="histoName", default="dijet_mass")
(options, args) = parser.parse_args()

if __name__ == "__main__":
    main(options, args)
