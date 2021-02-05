from multiprocessing import Process
import os
import time
# git remote set-url origin https://mgrecu35@github.com/mgrecu35/cmbv7.git

def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())
    
def fsh(fname_p):
    fname=fname_p[0]
    p=fname_p[1]
    lines=open(fname,"r").readlines()
    print(lines[-3])
    logName="log."+lines[-5][9:]
    cmd='bpsh %i ./combAlg.exe junk %s >&out/out.%s'%(247-p,fname,logName)
    print(cmd)
    os.system(cmd)
    #time.sleep(1)
import glob
t11=time.time()
if __name__ == '__main__':
    fs=glob.glob("paramF/*")
    fs=sorted(fs)
    for i in range(int(len(fs)/10)-1):
        jobs=[]
        t1=time.time()
        for k in range(10):
            p = Process(target=fsh, args=((fs[10*i+k],k),))
            jobs.append(p)
            p.start()
        for j in jobs:
            j.join()
        print('all done')
        print(time.time()-t1)
    print(time.time()-t11)
