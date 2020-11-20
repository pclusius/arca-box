#!/bin/bash

here=$(pwd)
echo "#!/bin/bash

cd $here

python3 ModelLib/gui/run.py
" > run_arca.sh

chmod a+x run_arca.sh

#cp $here/ModelLib/gui/thebox_ico.png ~/.icons/arca.png

echo "[Desktop Entry]
Name=ARCA Box Model 0.9
Comment=Atmospherically Relevant Chemistry and Aerosol Box model
Exec=$here/run_arca.sh
Icon=$here/ModelLib/gui/thebox_ico.png
Terminal=true
Type=Application
Categories=GNOME;GTK;Core;
StartupWMClass=ARCAbox utility" > arca.desktop
