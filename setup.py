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
make=open('Makefile','r+')
lines=make.readlines()
for i in range(len(lines)):
	if 'GSL_INCLUDE_DIR=' in lines[i]:
		lines[i]='GSL_INCLUDE_DIR='+INCLUDE_PATH+'\n'

for i in range(len(lines)):
	if 'GSL_LIBRARIES=' in lines[i]:
		lines[i]='GSL_LIBRARIES=-L'+LD_LIBRARY_PATH+' -Wl,\"-R '+LD_LIBRARY_PATH+'\"\n'

make.seek(0)
make.writelines(lines)
make.close()
if os.path.isfile('ankflag')==True:
	os.system('make clean')
os.system('make')

np.save('LDPATH',LD_LIBRARY_PATH)
os.chdir(cwd)
try:
	import casatools
	print ('casatools==6.0.0.27 is already installed\n')
except:
	os.system('python3 -m pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools==6.0.0.27 --user')

try:
	import casatasks
	print ('casatasks==6.0.0.27 is already installed\n')
except:
	os.system('python3 -m pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatasks==6.0.0.27 --user')

os.system('rm -rf casa*log')

setup(
    name='aNKflag',
    version='1.0',
    packages=['aNKflag'],
    author='Apurba Bera, Python wrapper by Devojyoti Kansabanik',
    description='Flagger',
    install_requires=["numpy", "astropy", "skyfield", "matplotlib", "scipy>=0.15.1"],
    extras_require={'skymap':["ephem", "Pillow"]}   # Needed only to generate sky maps in mwa_pb/skymap.py
)


