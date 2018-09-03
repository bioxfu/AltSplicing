conda create -y -n spanki pip parallel samtools=0.1.19 python=2.7
source activate spanki
pip install numpy

cd Spanki_0.5.0_xfu
python setup.py install
