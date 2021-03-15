import os,numpy as np,sys
from setuptools import setup


python_version=int(sys.version.split(' ')[0].split('.')[0])
print ('Python version :'+str(sys.version.split(' ')[0]))
if python_version!=3:
	print ('Python version is less than 3. aNKflag can only run with python>3\n')	
	os._exit(0)

cwd=os.getcwd()
LD_LIBRARY_PATH=cwd+'/gsl/lib'
INCLUDE_PATH=cwd+'/gsl/include/'

os.chdir(cwd+'/aNKflag')
makefil=open('Makefile','r')
lines=makefil.readlines()
for i in range(len(lines)):
	if 'GSL_INCLUDE_DIR=' in lines[i]:
		lines[i]='GSL_INCLUDE_DIR='+INCLUDE_PATH+'\n'

for i in range(len(lines)):
	if 'GSL_LIBRARIES=' in lines[i]:
		lines[i]='GSL_LIBRARIES=-L'+LD_LIBRARY_PATH+' -Wl,\"-R '+LD_LIBRARY_PATH+'\"\n'

makefil.close()

lines=lines[:18]
with open("Makefile", "w") as output:
    for line in lines:
        output.write(line)
output.close()


if os.path.isfile('ankflag')==True:
	os.system('make clean')
os.system('make')
np.save('LDPATH',LD_LIBRARY_PATH)
os.chdir(cwd)
try:
	import casatools
	print ('casatools is already installed\n')
except:
	os.system('python3 -m pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools --user')

try:
	import casatasks
	print ('casatasks is already installed\n')
except:
	os.system('python3 -m pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatasks --user')

os.system('rm -rf casa*log')

setup(
    name='aNKflag',
    version='1.0',
    packages=['aNKflag'],
	package_data={'aNKflag':['aNKflag/*.c', 'aNKflag/*.npy', 'aNKflag/ankflag', 'aNKflag/*.h','aNKflag/*.dat']},
    author='Apurba Bera, Python wrapper by Devojyoti Kansabanik',
    description='Flagger',
    install_requires=["numpy", "astropy", "matplotlib", "scipy>=0.15.1"],
    )
final_library_path='~/.local/lib/python3.6/site-packages/aNKflag-1.0-py3.6.egg/aNKflag'
os.system('cp -r '+cwd+'/aNKflag/*.c '+cwd+'/aNKflag/*.h '+cwd+'/aNKflag/ankflag '+cwd+'/aNKflag/*.dat '+cwd+'/aNKflag/*.npy '+final_library_path)

