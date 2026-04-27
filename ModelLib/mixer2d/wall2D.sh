#/usr/bash/bin

[ "$1" == 'init' ] && CONT=1
[ "$1" == '' ] && CONT=0
[ "$1" == 'continue' ] && CONT=0
[ "$1" == 'post' ] && CONT="post"

THIS=$(pwd)

[ -z ${CONT+x} ] && echo Missing CONT && exit 1
[ -z ${locationPlumeData+x} ] && echo Missing locationPlumeData && exit 1
[ -z ${initialisationSettings+x} ] && echo Missing initialisationSettings && exit 1
[ -z ${plumeInitialisationSettings+x} ] && echo Missing plumeInitialisationSettings && exit 1
[ -z ${locInitialBinary+x} ] && echo Missing locInitialBinary && exit 1
[ -z ${locPlumeBinary+x} ] && echo Missing locPlumeBinary && exit 1
[ -z ${backgroundContSettings+x} ] && echo Missing backgroundContSettings && exit 1
[ -z ${bottomContSettings+x} ] && echo Missing bottomContSettings && exit 1
[ -z ${dx+x} ] && echo Missing dx && exit 1
[ -z ${ddz+x} ] && echo Missing ddz && exit 1
[ -z ${nVerticalLayers+x} ] && echo Missing nVerticalLayers && exit 1
[ -z ${nHorisontalColumns+x} ] && echo Missing nHorisontalColumns && exit 1
[ -z ${zIndexPlumeStr+x} ] && echo Missing zIndexPlumeStr && exit 1
[ -z ${yIndexPlumeStr+x} ] && echo Missing yIndexPlumeStr && exit 1
[ -z ${tLapseRate+x} ] && echo Missing tLapseRate && exit 1
[ -z ${pLapseRate+x} ] && echo Missing pLapseRate && exit 1
[ -z ${mixingTimestep+x} ] && echo Missing mixingTimestep && exit 1
[ -z ${rounds+x} ] && echo Missing rounds && exit 1
[ -z ${boxRuntime+x} ] && echo Missing boxRuntime && exit 1
[ -z ${scaleConcWithAltitude+x} ] && echo Missing scaleConcWithAltitude && exit 1
[ -z ${namelist+x} ] && echo Missing namelist && exit 1
[ -z ${CaseCells+x} ] && echo Missing CaseCells && exit 1

initialisationSettings=${locationPlumeData}/${initialisationSettings}
plumeInitialisationSettings=${locationPlumeData}/${plumeInitialisationSettings}
backgroundContSettings=${locationPlumeData}/${backgroundContSettings}
bottomContSettings=${locationPlumeData}/${bottomContSettings}
locInitialBinary=${locationPlumeData}/${locInitialBinary}
locPlumeBinary=${locationPlumeData}/${locPlumeBinary}
# figure out NCHEM
tmp=$(grep NAMESDAT $initialisationSettings)
namesdatPath=$(sed s/\'//g <<< $tmp |awk '{print $3}')
SPC_NAMES_Path=$( dirname ${arca}/$namesdatPath)/SPC_NAMES.txt
NCHEM=$(cat $SPC_NAMES_Path | wc -l )
tmp=$(grep N_BINS_PAR $initialisationSettings)
nBins=$(sed s/\'//g <<< $tmp |awk '{print $3}')
N_dilute=$(( $NCHEM + $nBins ))

echo
echo "Settings for this simulation:"
echo "------------------------------------------------------------------------------------------------------------------"
echo "Initialize?                : $CONT"
echo "locationPlumeData          : $locationPlumeData"
echo "initialisationSettings     : $initialisationSettings"
echo "plumeInitialisationSettings: $plumeInitialisationSettings"
echo "locInitialBinary           : $locInitialBinary"
echo "locPlumeBinary             : $locPlumeBinary"
echo "backgroundContSettings     : $backgroundContSettings"
echo "bottomContSettings         : $bottomContSettings"
echo "dx                         : $dx"
echo "ddz                        : $ddz"
echo "nVerticalLayers            : $nVerticalLayers"
echo "nHorisontalColumns         : $nHorisontalColumns"
echo "zIndexPlumeStr             : $zIndexPlumeStr"
echo "yIndexPlumeStr             : $yIndexPlumeStr"
echo "tLapseRate                 : $tLapseRate"
echo "pLapseRate                 : $pLapseRate"
echo "mixingTimestep             : $mixingTimestep"
echo "rounds                     : $rounds"
echo "boxRuntime                 : $boxRuntime"
echo "Total simulation time:     : $(( rounds * boxRuntime ))"
echo "KyKz                       : $KyKz"
echo "randomize_k                : $randomize_k"
echo "friction_vel               : $friction_vel"
echo "pblh                       : $pblh"
echo "scaleConcWithAltitude      : $scaleConcWithAltitude"
echo "namelist                   : $namelist"
echo "CaseCells                  : $CaseCells"
echo "openMP                     : $openMP"
IFS=', ' read -r -a yIndexPlume <<< "$yIndexPlumeStr"
IFS=', ' read -r -a zIndexPlume <<< "$zIndexPlumeStr"


MIXERSRC=${arca}/ModelLib/mixer2d
if [ $openMP == '1' ] ; then
  export OMP_NUM_THREADS=$(nproc --all)
fi


bginit=$initialisationSettings
plumeinit=$plumeInitialisationSettings
columnSettingDir=${bginit}/settings/column/
cells=$locationPlumeData/${CaseCells}
r16bin=CFINAL_0000.r16
bginitDir=$(dirname $bginit)
plumeinitDir=$(dirname $plumeinit)
mkdir -p ${cells}
if [ $CONT == 'post' ] ; then
  for d in $(seq ${nVerticalLayers}); do
  for c in $(seq ${nHorisontalColumns}); do
    s=${cells}/$( printf "%03d" $d )/$( printf "%03d" $c )
    cat $s/CFINAL_* > $s/c_ts.r16 ;
  done ;
  done

  python ${arca}/ModelLib/required/readWriteBin.py 1 || exit
  exit
fi

if [ "$CONT" == '1' ] ; then

rm -f ${locInitialBinary}/CFINAL*
rm -f ${locPlumeBinary}/CFINAL*
cd ${arca}
./arcabox.exe ${bginit}
./arcabox.exe ${plumeinit}
cd $THIS

NSPEC=$(<${locPlumeBinary}/${r16bin} wc -c)
NSPEC=$(( $NSPEC / 8 ))
tmpPyFile=F209fn0932f_OJ9OIDWNm94321b.py

echo "------------------------------------------------------------------------------------------------------------------"
dim=$(( $nVerticalLayers * $nHorisontalColumns ))
echo "Number of scalars to mix   : $NSPEC"
echo "Number of cells to mix     : $dim"

echo "Mixing plume with background"
echo """
from numpy import fromfile, array, float64
import sys
def read_bin(file):
    f = open(file, 'rb')
    c = fromfile(f, dtype=float64)
    f.close()
    return c
def write_bin(file,c):
    f = open(file, 'wb')
    f.write(c);f.close(); return
ndil = ${N_dilute}
infiles  = [f for f in sys.argv[2:]]
c = read_bin(infiles[0])
for f in infiles[1:]:
    c[:ndil] = c[:ndil] + read_bin(f)[:ndil]
write_bin(sys.argv[1],c)
""">${tmpPyFile}
python ${tmpPyFile} \
  ${locPlumeBinary}/${r16bin}   \
  ${locInitialBinary}/${r16bin} \
  ${locPlumeBinary}/${r16bin}
echo "Done mixing "

# echo mkdir -p $cells # DONE already at this point

Tinit=$(grep "MODS.*TEMPK" ${bottomContSettings} |awk '{print $6}')
Pinit=$(grep "MODS.*PRESSURE" ${bottomContSettings} |awk '{print $6}')
T0=${Tinit/d/e}
P0=${Pinit/d/e}

echo """
from numpy import linspace,savetxt,exp,append,cumsum
dz=linspace(0,(${nVerticalLayers}-2)*${ddz},${nVerticalLayers}-1)
h = append(0,cumsum(dz+${dx}))/1000
print(h)
savetxt('${locationPlumeData}/settings/t.txt',( ${T0} - h*${tLapseRate} ).T,fmt='%8.5e')
""" >$tmpPyFile
if [ $pLapseRate == 'b' ] ; then
  echo """
pnz=${P0}*exp(-h/8.431)
savetxt('${locationPlumeData}/settings/p.txt',pnz.T,fmt='%8.5e')""" >>$tmpPyFile
else
echo """
savetxt('${locationPlumeData}/settings/p.txt',( ${P0}-h*${pLapseRate}   ).T,fmt='%8.5e')
""" >>$tmpPyFile
fi

python $tmpPyFile

if [ "$scaleConcWithAltitude" == '1' ] ; then
echo """
from numpy import fromfile, array, float64
import sys
def read_bin(file):
    f = open(file, 'rb')
    return fromfile(f, dtype=float64)
def write_bin(file,c):
    f = open(file, 'wb')
    f.write(c);f.close()
t0 = float(sys.argv[1])
p0 = float(sys.argv[2])
t1 = float(sys.argv[3])
p1 = float(sys.argv[4])
f  = sys.argv[5]
kelvin = 0 if t0>150 else 273.15
factor = (p1*(t0+kelvin))/((t1+kelvin)*p0)
#if kelvin>0:print('Assuming C°', factor,f)
c = read_bin(f)
c[:$N_dilute] = c[:$N_dilute] * factor
write_bin(f,c)
""" > $tmpPyFile
fi
echo "Processing output directories and data..."
rm -rf ${cells}/*

mkdir -p ${cells}/k-values
rm -f ${cells}/k-values/*

echo > onecycle.sh
for d in $(seq ${nVerticalLayers}); do
  Temp=$(sed -n ${d}p ${locationPlumeData}/settings/t.txt )
  Pres=$(sed -n ${d}p ${locationPlumeData}/settings/p.txt )
  dir=${cells}/$( printf "%03d" $d )/001
  mkdir -p $dir
  rm -f $dir/*
  binary=$dir/${r16bin}
  cp ${locInitialBinary}/${r16bin} $binary
  if [ $scaleConcWithAltitude == '1' ] ;then
    python $tmpPyFile $T0 $P0 $Temp $Pres $binary &
  fi
done
wait

for d in $(seq ${nVerticalLayers}); do
  Temp=$(sed -n ${d}p ${locationPlumeData}/settings/t.txt )
  Pres=$(sed -n ${d}p ${locationPlumeData}/settings/p.txt )
  for c in $(seq ${nHorisontalColumns}); do
    s=${cells}/$( printf "%03d" $d )/$( printf "%03d" $c )
    mkdir -p $s
    # Clear directory, start copying stuff
    if [ "$c" != '1' ] ; then
      rm -f $s/*
      cp ${cells}/$( printf "%03d" $d )/001/${r16bin} $s/. &
    fi
    if printf '%s\0' "${zIndexPlume[@]}" | grep -xqnF -- $d && \
    printf '%s\0' "${yIndexPlume[@]}" | grep -xqnF -- $c ; then
      wait
      echo
      echo copying plume file from ${locPlumeBinary} to $s
      cp ${locPlumeBinary}/${r16bin} $s/.
      [ $scaleConcWithAltitude == '1' ] && python $tmpPyFile $T0 $P0 $Temp $Pres $s/${r16bin}
    fi

    # if [ $d == $zIndexPlume ] && [ $c == $yIndexPlume ] ; then
    #   wait
    #   echo
    #   echo copying plume file from ${locPlumeBinary} to $s
    #   cp ${locPlumeBinary}/${r16bin} $s/.
    #   [ $scaleConcWithAltitude == '1' ] && python $tmpPyFile $T0 $P0 $Temp $Pres $s/${r16bin}
    # fi

    if [ $d == '2' ] ; then
      initfile=${bottomContSettings}
    else
      initfile=${backgroundContSettings}
    fi
    cp ${initfile} $s/cont.init
    sed -i "s|<XINITY>|${s}|g" ${s}/cont.init
    sed -i "s|MODS(1)   = -1  -1  1.00000d+00.*${Tinit}|MODS(1)   = -1  -1  1.00000d+00  ${Temp}|g" ${s}/cont.init
    sed -i "s|MODS(2)   = -1  -1  1.00000d+00.*${Pinit}|MODS(2)   = -1  -1  1.00000d+00  ${Pres}|g" ${s}/cont.init
    echo ./arcabox.exe ${s}/cont.init \& >> onecycle.sh
    # echo ${s}/cont.init \& >> onecycle.sh
  done
  printf '%s' "$d.."
  echo wait >> onecycle.sh
done
echo
echo done

rm $tmpPyFile


fi
[ ! -f ${cells}/001/001/${r16bin} ] && echo "!!! No simulation to continue !!!" && exit

NSPEC=$(<${cells}/001/001/${r16bin} wc -c)
NSPEC=$(( $NSPEC / 8 ))

update_nml() {
echo """&NML_meteo
 KyKz         = ${1}
 randomize_k  = ${2}
 friction_vel = ${3}
 pblh         = ${4}
/
""" > ${locationPlumeData}/meteo.namelist
}

if [ $namelist == '1' ] ; then
  use_nml="use_namelist=1"
  update_nml ${KyKz} ${randomize_k} ${friction_vel} ${pblh}
  meteonamelistfile=$(pwd)/meteo.namelist
else
  use_nml="use_KyKz=${KyKz} use_randomize_k=${randomize_k} use_friction_vel=${friction_vel} use_pblh=${pblh}"
fi
[ $openMP == '1' ] && use_nml+=" use_omp=1"

cd $MIXERSRC
make fdx=$dx fddz=$ddz nVerticalLayers=$nVerticalLayers nHorisontalColumns=$nHorisontalColumns \
  mixingTimestep=$mixingTimestep totalRuntime=$boxRuntime nSpec=$NSPEC Path=${cells} \
  yIndexPlume=1 zIndexPlume=1 $use_nml || exit
mv mixer2d.exe $THIS
cd $THIS

if [ "$CONT" == '1' ]; then
  time ./mixer2d.exe 0 0 $meteonamelistfile
fi

cd ${arca}

for i in $(seq $rounds )
do
  tot=$(ls ${cells}/001/001/CFINAL* |wc -l)
  echo Time: $(( tot * boxRuntime )) seconds
  echo  ${i} \(${tot}\): Running $dim ARCA-boxes
  time bash $THIS/onecycle.sh > /dev/null
  # time cat $THIS/onecycle.sh |xargs -n 1 -P $OMP_NUM_THREADS ./arcabox.exe >/dev/null
  [ -f ${locationPlumeData}/settings/blh.txt ] && update_nml ${KyKz} ${randomize_k} ${friction_vel} $(sed -n ${tot}p ${locationPlumeData}/settings/blh.txt )
  echo  ${i} \($tot\): mixing concentrations
  time ${locationPlumeData}/mixer2d.exe 0 0 $meteonamelistfile
[ -d 'stop' ] && rmdir stop && break
done

cd $THIS

for d in $(seq ${nVerticalLayers}); do
for c in $(seq ${nHorisontalColumns}); do
  s=${cells}/$( printf "%03d" $d )/$( printf "%03d" $c )
  cat $s/CFINAL_* > $s/c_ts.r16 ;
done ;
done

python ${arca}/ModelLib/required/readWriteBin.py 1
#



exit
