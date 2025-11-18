# 1st, in powerShell, install Ubuntu:
# wsl --install -d Ubuntu
echo "-------------------------------------------------------------------------------"
echo "This script will install the Fortran, MPI and NetCDF packages and the"
echo "systemwide Python modules necessary for ARCA (and other typical Python scripts)"
echo ""
yesno=""

while getopts 'y' flag; do
  case "${flag}" in
    y) yesno='--yes'
       echo "will use apt --assume_yes" ;;
    *) echo "only one option, -y if you want to skip the y/n questions of apt (using apt --assume_yes)"
       exit 1 ;;
  esac
done

read -p "Continue? (y/n)? " answer
case ${answer:0:1} in
    y|Y )
        echo  "OK! This needs your sudo pass."
    ;;
    * )
        exit
    ;;
esac

cd ~
sudo apt-get --yes update
if [ "$EUID" -ne 0 ]
  then echo "root password is needed to install packages"
  exit
fi
sudo apt-get --yes install xdg-utils
sudo apt-get --yes install desktop-file-utils
sudo apt-get --yes install make
sudo apt-get --yes install gfortran
sudo apt --yes install zlib1g-dev
sudo apt-get --yes install zip
sudo apt-get --yes install g++
sudo apt-get --yes install hdf5-tools
sudo apt-get --yes install libnetcdff-dev
sudo apt-get --yes install libnetcdf-mpi-dev
sudo apt-get --yes install netcdf-bin
sudo apt-get --yes install python3-matplotlib
sudo apt-get --yes install python3-pyqt5
sudo apt-get --yes install python3-pyqtgraph
sudo apt-get --yes install python3-netcdf4
git clone --depth 1 https://github.com/pclusius/arca-box.git
cd arca-box/install/
python3 setup.py
cd ..
sh run_arca.sh
