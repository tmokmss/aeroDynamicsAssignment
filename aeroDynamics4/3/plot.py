import csv
import pylab as pl
import numpy as np
import sys

filename = 'result.csv'
if (len(sys.argv)>1):
  filename = sys.argv[1]

reader = csv.reader((open(filename,'r')))
x = []
rho = []
u = []
p = []
e = []
rowname = reader.next()
print rowname
for row in reader:
  x.append(float(row[0]))
  rho.append(float(row[1]))
  u.append(float(row[2]))
  p.append(float(row[3]))
  e.append(float(row[4]))

def read_t(csvname):
  treader = csv.reader((open(csvname,'r')))
  vt = []
  tt = []
  for row in treader:
    vt.append([])
    tt.append(float(row[0]))
    for col in row[1:]:
      vt[-1].append(float(col))
  return vt, tt

def plot_contour(Z, fname, title):
  xm, tm = np.meshgrid(x, time)
  pl.contourf(xm, tm, Z,100)
  pl.title(title, fontsize=20)
  pl.colorbar()
  pl.xlabel('x')
  pl.ylabel('t [sec]')
  pl.savefig(fname+'.png')
  pl.clf()

ut, time = read_t('ut38.csv')
rhot, time = read_t('rhot38.csv')
pt, time = read_t('pt38.csv')
et, time = read_t('et38.csv')

plot_contour(ut,'ut','Time and Velocity')
plot_contour(rhot, 'rhot', 'Time and Density')
plot_contour(et, 'et', 'Time and Temperature')
plot_contour(pt, 'pt', 'Time and Pressure')

idx = int(len(x)*0.4)
pl.plot(time, [rhoi[idx] for rhoi in rhot], label='dens.')
pl.plot(time, [pi[idx] for pi in pt], label='pres.')
pl.title('Changes along time at x = ' + str(x[idx]))
pl.xlabel('time [sec]')
pl.ylabel('value')
pl.xlim(time[0],time[-1])
pl.legend()
pl.savefig('rhoptime.eps')
pl.show()

pl.plot(time, [ui[idx] for ui in ut], label='veloc.')
pl.plot(time, [ei[idx] for ei in et], label='temp.')
pl.title('Changes along time at x = ' + str(x[idx]))
pl.xlabel('time [sec]')
pl.ylabel('value')
pl.xlim(time[0],time[-1])
pl.legend()
pl.savefig('uttime.eps')
pl.show()

timeidx = int(len(time)*0.6)
pl.plot(x, rhot[timeidx], label='dens.')
pl.plot(x, pt[timeidx], label='pres.')
pl.legend()
pl.xlim([x[0],x[-1]])
pl.xlabel('x')
pl.ylabel('value')
pl.title('Changes along x at time = ' + str(time[timeidx])+' sec')
pl.savefig('rhopx.eps')
pl.show()
pl.plot(x, ut[timeidx], label='veloc.')
pl.plot(x, et[timeidx], label='temp.')
pl.legend()
pl.xlim([x[0],x[-1]])
pl.xlabel('x')
pl.ylabel('value')
pl.title('Changes along x at time = ' + str(time[timeidx])+' sec')
pl.savefig('utx.eps')
pl.show()

