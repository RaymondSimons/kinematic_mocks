from joblib import Parallel, delayed
import time
import astropy

def test(x = 5):
    for i in range(100000000): pass
    print x
    f = open('./temp/'+str(x), 'w+')
    f.close()

def test2(x = 10):
    for i in range(100000000): pass
    print x
    f = open('./temp/new_'+str(x), 'w+')
    f.close()

a = time.time()
Parallel(n_jobs = -1)(delayed(test)(i) for i in range(32))
b = time.time()
print 'Finished parallel in:', b - a

'''
e = time.time()
Parallel(n_jobs = 1)(delayed(test2)(i) for i in range(16))
f = time.time()
print 'Finished series in:', f - e
'''