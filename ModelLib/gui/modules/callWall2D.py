import numpy as np
import os,re,shutil

def start2D(**kwds):
    # os.system(f'mkdir -p {kwds['locationPlumeData']}/cells')
    os.makedirs(f'{kwds['locationPlumeData']}/cells', exist_ok=True)
    os.makedirs(f'{kwds['locationPlumeData']}/settings', exist_ok=True)
    filesTo  = kwds['locationPlumeData']
    initial  = kwds['initialisationSettings']
    plumeIni = kwds['plumeInitialisationSettings']
    BGcont   = kwds['backgroundContSettings']
    BtmCont  = kwds['bottomContSettings']

    """Copy files to 'filesTo' """
    collect = []
    for i,source, target in zip(range(4),[initial,plumeIni,BGcont,BtmCont],['Initial','PlumeIni','BgCont','BtmCont']):
        shutil.copyfile(source, filesTo+f'/settings/{target}.init')
        collect.append(filesTo+f'/settings/{target}.init')

    """PreProcess the new files: Runtime, INIT_W_BIN """

    initial = collect[0] #  = kwds['initialisationSettings']
    plumeIni = collect[1] # = kwds['plumeInitialisationSettings']
    BGcont = collect[2] #   = kwds['backgroundContSettings']
    BtmCont = collect[3] #  = kwds['bottomContSettings']
    locInitialBinary = f'{kwds['locationPlumeData']}/initial/background/'
    locPlumeBinary = f'{kwds['locationPlumeData']}/initial/plume/'
    os.makedirs(locInitialBinary, exist_ok=True)
    os.makedirs(locPlumeBinary, exist_ok=True)
    for source, target in zip([BGcont,BtmCont],['BgCont','BtmCont']):
        with open(source, 'r') as file:
            text = file.read()
            """Make sure model timestep does not exceed new runtime"""
            dtf = re.findall(r'DT\s*?=\s*?(\d*(\.\d+)?)\s*?\n',text)
            if len(dtf)>0:
                dt = float(dtf[0][0])
                if dt>float(kwds['chemTimestep']):
                    text = re.sub(
                        r'DT\s*?=\s*?(\d*(\.\d+)?)\s*?\n',
                        f'DT = {float(kwds['chemTimestep']):0.1f}\n',
                        text,count=1,flags=re.IGNORECASE)
            """insert new runtime"""
            text = re.sub(r'RUNTIME\s*?=\s*?(\d*(\.\d+)?)\s*?\n',f'RUNTIME = {float(kwds['chemTimestep']):0.1f}\n',text,count=1,flags=re.IGNORECASE)
            """Make sure these are seconds"""
            text = re.sub(r'RUNTIME_TIME_UNIT\s*?=\s*?\'(min|sec|hrs|day)\'\s*?\n',f'RUNTIME_TIME_UNIT = \'sec\'\n',text,count=1,flags=re.IGNORECASE)
            """replace binary output path with a placecholder, or write it if not present"""
            text,r = re.subn(r'INIT_W_BIN\s*?=\s*?\'[0-9a-zA-Z<>/\-_]*?\'\s*?\n',f'INIT_W_BIN = \'<XINITY>/\'\n',text,count=1,flags=re.IGNORECASE)
            if r==0:
                text = re.sub(r'&NML_CUSTOM\s*?\n',f'&NML_CUSTOM\n INIT_W_BIN = \'<XINITY>/\'\n',text,count=1,flags=re.IGNORECASE)
            """Don't write extrafiles, add option if not there """
            text,r = re.subn(r'EXTRAFILES\s*?=\s*?\.true\.\s*?\n',f'EXTRAFILES = .false.\n',text,count=1,flags=re.IGNORECASE)
            if r==0:
                text = re.sub(r'&NML_CUSTOM\s*?\n',f'&NML_CUSTOM\n EXTRAFILES = .false.\n',text,count=1,flags=re.IGNORECASE)
            """Don't write NetCDF, add option if not there """
            text,r = re.subn(r'NETCDF_OUT\s*?=\s*?\.true\.\s*?\n',f'NETCDF_OUT = .false.\n',text,count=1,flags=re.IGNORECASE)
            if r==0:
                text = re.sub(r'&NML_CUSTOM\s*?\n',f'&NML_CUSTOM\n NETCDF_OUT = .false.\n',text,count=1,flags=re.IGNORECASE)
            """Turn off verbose output  VERBOSE = .TRUE."""
            text,r = re.subn(r'VERBOSE\s*?=\s*?\.TRUE\.\s*?\n',f'VERBOSE = .FALSE.\n',text,count=1,flags=re.IGNORECASE)
        with open(source, 'w') as file:
            file.write(text)

    for source, target, bindir in zip([initial,plumeIni],['Initial','PlumeIni'],[locInitialBinary,locPlumeBinary]):
        with open(source, 'r') as file:
            text = file.read()

            """replace binary output path with a placecholder, or write it if not present"""
            text,r = re.subn(r'INIT_W_BIN\s*?=\s*?\'[0-9a-zA-Z<>/\-_]*?\'\s*?\n',f'INIT_W_BIN = \'{bindir}\'\n',text,count=1,flags=re.IGNORECASE)
            if r==0:
                text = re.sub(r'&NML_CUSTOM\s*?\n',f'&NML_CUSTOM\n INIT_W_BIN = \'{bindir}\'\n',text,count=1,flags=re.IGNORECASE)


        with open(source, 'w') as file:
            file.write(text)
    def relpath(p):
        return os.path.relpath(p, filesTo)
    paths_relative_to_locationPlumeData = True
    rounds = int(float(kwds['totalRuntime'])/float(kwds['chemTimestep']))
    file = open(f'{kwds['locationPlumeData']}/wall2dSettings.sh','w+')
    file.write('#/usr/bin/bash\n')
    file.write(f'export locationPlumeData={kwds['locationPlumeData']}\n')
    if paths_relative_to_locationPlumeData:
        file.write(f'export mixer=mixer2D.exe\n')
        file.write(f'export arca={relpath(kwds['arca'])}\n')
        file.write(f'export initialisationSettings={relpath(initial)}\n')
        file.write(f'export plumeInitialisationSettings={relpath(plumeIni)}\n')
        file.write(f'export backgroundContSettings={relpath(BGcont)}\n')
        file.write(f'export bottomContSettings={relpath(BtmCont)}\n')
        file.write(f'export locInitialBinary={relpath(locInitialBinary)}\n')
        file.write(f'export locPlumeBinary={relpath(locPlumeBinary)}\n')
    else:
        file.write(f'export mixer={kwds['locationPlumeData']}/mixer2D.exe\n')
        file.write(f'export arca={kwds['arca']}\n')
        file.write(f'export initialisationSettings={initial}\n')
        file.write(f'export plumeInitialisationSettings={plumeIni}\n')
        file.write(f'export backgroundContSettings={BGcont}\n')
        file.write(f'export bottomContSettings={BtmCont}\n')
        file.write(f'export locInitialBinary={locInitialBinary}\n')
        file.write(f'export locPlumeBinary={locPlumeBinary}\n')
    #
    file.write(f'export dx={kwds['dx']}\n')
    file.write(f'export nVerticalLayers={kwds['nVerticalLayers']}\n')
    file.write(f'export nHorisontalColumns={kwds['nHorisontalColumns']}\n')
    file.write(f'export zIndexPlume={kwds['zIndexPlume']}\n')
    file.write(f'export yIndexPlume={kwds['yIndexPlume']}\n')
    file.write(f'export tLapseRate={kwds['tLapseRate']}\n')
    file.write(f'export pLapseRate={kwds['pLapseRate']}\n')
    file.write(f'export mixingTimestep={kwds['mixingTimestep']}\n')
    file.write(f'export boxRuntime={kwds['chemTimestep']}\n')
    file.write(f'export KyKz={kwds['KyKz']}\n')
    file.write(f'export rounds={rounds:d}\n')
    file.write(f'export randomize_k={kwds['randomize_k']}\n')
    file.write(f'export scaleConcWithAltitude={kwds['scaleConcWithAltitude']}\n')
    file.write(f'export friction_vel={kwds['friction_vel']}\n')
    file.write(f'export pblh={kwds['pblh']}\n')
    file.write(f'export namelist={kwds['namelist']}\n')
    file.write('\n')
    file.write('[ "$1" == "init" ] || [ "$1" == "continue" ] || [ "$1" == "post" ] && bash ${arca}/ModelLib/mixer2d/wall2D.sh $@')
    file.close()

if __name__ == "__main__":
    print('Call start2D() with necessary kw arguments to write a 2D simulation setup and caller script')







#
