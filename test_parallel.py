from joblib import Parallel, delayed
import time


def test(x = 5):
    for i in range(30000000): pass
    print x


a = time.time()
Parallel(n_jobs = -1)(delayed(test)(i) for i in range(10))
b = time.time()
print 'Finished in:', b - a


a = time.time()
Parallel(n_jobs = 1)(delayed(test)(i) for i in range(10))
b = time.time()
print 'Finished in:', b - a