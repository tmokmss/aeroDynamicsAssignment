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
  pl.contour(xm, tm, Z)
  pl.title(title, fontsize=20)
  pl.colorbar()
  pl.show()

ut, time = read_t('ut38.csv')
rhot, time = read_t('rhot38.csv')
pt, time = read_t('pt38.csv')
et, time = read_t('et38.csv')

plot_contour(ut,'ut.png','Time and Velocity')
plot_contour(rhot, 'rhot.png', 'Time and Density')
plot_contour(et, 'et,png', 'Time and Temperature')

pl.plot(x, rho, label='density')
pl.plot(x, u, label='u')
pl.plot(x, p, label='p')
pl.plot(x, e, label='T')
pl.legend()
pl.xlim([x[0],x[-1]])
pl.show()
