import csv
import pylab as pl
import sys

filename = 'result.csv'
if (len(sys.argv)>1):
  filename = sys.argv[1]

reader = csv.reader((open(filename,'r')))
j = []
exact = []
numerical = []
rowname = reader.next()
numeNum = 2
print rowname 
for i in xrange(numeNum):
  numerical.append([])
  exact.append([])
for row in reader:
  j.append(float(row[0]))
  exact[0].append(float(row[3]))
  numerical[0].append(float(row[4]))
  exact[1].append(float(row[5]))
  numerical[1].append(float(row[6]))


pl.plot(j, exact[0], label='Exact u')
pl.plot(j, exact[1], label='Exact v')
pl.plot(j, numerical[0], label='Numerical u')
pl.plot(j, numerical[1], label='Numerical v')
pl.legend()
pl.xlim([j[0],j[-1]])
pl.show()
