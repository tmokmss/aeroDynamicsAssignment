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

#  calc errors
def calc_error(exact, nume):
  errorsum = 0
  for ei, ni in zip(exact, nume):
    error = (ei-ni)
    if (abs(ei) > 0.0000001):
      error /= ei
    errorsum += abs(error)
  return errorsum/len(exact)

print 'calculated error:'
for nume, row in zip(numerical, rowname[3:]):  
  print row, calc_error(exact, nume)

pl.plot(j, exact, label=' Exact')
for nume, row in zip(numerical, rowname[3:]):
  pl.plot(j, nume, label=row)
pl.legend()
pl.xlim([j[0],j[-1]])
pl.xlabel('x', fontsize=20)
pl.ylabel('u', fontsize=20)
pl.savefig(filename+'.png')
pl.savefig(filename+'.eps')
pl.show()


