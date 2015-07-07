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
init = []
rowname = reader.next()
numeNum = 2
print rowname 
for i in xrange(numeNum):
  numerical.append([])
  exact.append([])
  init.append([])
for row in reader:
  j.append(float(row[0]))
  init[0].append(float(row[1]))
  init[1].append(float(row[2]))
  exact[0].append(float(row[3]))
  numerical[0].append(float(row[5]))
  exact[1].append(float(row[4]))
  numerical[1].append(float(row[6]))

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
for nume, row, ex in zip(numerical, rowname[5:], exact):  
  print row, calc_error(ex, nume)


pl.plot(j, init[0], '--', label='Init u')
pl.plot(j, exact[0], label='Exact u')
pl.plot(j, numerical[0], label='Numerical u')
pl.legend()
pl.xlabel('x', fontsize=16)
pl.ylabel('variable', fontsize=16)
pl.xlim([j[0],j[-1]])
pl.savefig('result27_u.eps')
pl.savefig('result27_u.png')
pl.show()

pl.plot(j, init[1], '--',  label='Init v')
pl.plot(j, exact[1], label='Exact v')
pl.plot(j, numerical[1], label='Numerical v')
pl.legend()
pl.xlim([j[0],j[-1]])
pl.xlabel('x', fontsize=16)
pl.ylabel('variable', fontsize=16)
pl.savefig('result27_v.eps')
pl.savefig('result27_v.png')
pl.show()
