import numpy
import datetime
np_arr = numpy.load("../pixel-annealing/workdir/LS2_fpix_bmi.npy")
d0 = datetime.datetime.strptime('2018-12-02 16:10:07.222363', "%Y-%m-%d %H:%M:%S.%f")
dd = []
for npv in np_arr:
 dd.append(npv[0])

dd=numpy.array(dd)

dt = [datetime.datetime(2018, 12, 2, 16, 11, 6, 606000)-d0]
diffd=numpy.diff(dd)
for d in diffd:
 dt.append(d)

for i,t in enumerate(dt):
 dt[i] = t.total_seconds()

fout = open("../pixel-annealing/workdir/LS2_bmi_1.txt", "w+")
T=0
t0=0
for i,t in enumerate(dt):
 T+=t
 if abs(np_arr[i][1]+273.15-t0) > 1. or T > 86400:
  if int(round(T)) > 0: fout.write('%s\t273.15\t0\t%s\t-99.00\t-99.00\n' % (int(round(T)), np_arr[i][1]+273.15))
  else: pass
  T=0
 t0=np_arr[i][1]+273.15

fout.close()
