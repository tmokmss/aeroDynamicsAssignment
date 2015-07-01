import csv
import pylab as pl

filename = 'result.csv'
reader = csv.reader((open(filename,'r')))
j = []
exact = []
numerical = []
reader.next()
for row in reader:
  j.append(float(row[0]))
  exact.append(float(row[2]))
  numerical.append(float(row[3]))

pl.plot(j, exact, label='Exact')
pl.plot(j, numerical, label='Numerical')
pl.legend()
pl.show()
