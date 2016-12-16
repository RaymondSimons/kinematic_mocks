from joblib import Parallel, delayed
import time
import astropy

def test(x = 5):
    f = open('./temp/'+str(x), 'w+')
    f.close()

def test2(x = 10):
    f = open('./temp/new_'+str(x), 'w+')
    f.close()

a = time.time()
Parallel(n_jobs = -1)(delayed(test)(i) for i in range(10))
b = time.time()
print 'Finished in:', b - a


a = time.time()
Parallel(n_jobs = 1)(delayed(test2)(i) for i in range(10))
b = time.time()
print 'Finished in:', b - a