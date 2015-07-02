import csv
import pylab as pl
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

pl.plot(x, rho, label='density')
pl.plot(x, u, label='u')
pl.plot(x, p, label='p')
pl.plot(x, e, label='T')
pl.legend()
pl.xlim([x[0],x[-1]])
pl.show()
