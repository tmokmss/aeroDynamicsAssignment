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
numeNum = len(rowname)-3
print rowname 
for i in xrange(numeNum):
  numerical.append([])
for row in reader:
  j.append(float(row[0]))
  exact.append(float(row[2]))
  for i in xrange(numeNum):
    numerical[i].append(float(row[3+i]))

pl.plot(j, exact, label='Exact')
for nume in numerical:
  pl.plot(j, nume, label='Numerical')
pl.legend()
pl.xlim([j[0],j[-1]])
pl.show()
