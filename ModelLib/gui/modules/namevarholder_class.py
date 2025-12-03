class namevarholder():
    def __init__(self,path_to_names,path_to_xtras):
        self.NAMES = []
        self.namesPyInds = {}
        self.namesFoInds = {}
        self.divider_i = 0
        self.divider_xtr_i = 0
        self.path_to_names = path_to_names
        self.path_to_xtras = path_to_xtras
        self.load_names(path_to_names,path_to_xtras)

    def load_names(self,path_to_names,path_to_xtras):
        self.NAMES = []
        namesPyInds = {}
        namesFoInds = {}
        self.path_to_names = path_to_names
        self.path_to_xtras = path_to_xtras
        ## Create lists and dictionaries related to NAMES -------------
        tempbfr = []
        tempemi = []
        i = 0
        with open(path_to_names) as f:
            for line in f:
                if line[:3].upper() == 'EMI':
                    tempemi.append(line[:-1])
                else:
                    tempbfr.append(line[:-1])
                if '#' in line:
                    self.divider_i=i
                i += 1
        tempbfr[self.divider_i+1:] = sorted(tempbfr[self.divider_i+1:], key=str.lower)
        tempbfr += sorted(tempemi, key=str.lower)
        i = 0
        # j = 0
        for line in tempbfr:
            # j += 1
            name = line
            if name == '':
                continue
            # if i+1!=j:
            #     print('\n\nEmpty line in the middle of NAMES.DAT, it will not work!\n\n')
            #     exit('ARCA will not start, check out '+path_to_names+' and remove empty lines.')
            if '#' in line:
                self.divider_i=i
                name = 'Components in current chemistry'
            self.NAMES.append(name)
            namesPyInds[name] = i
            namesFoInds[name] = i+1
            i += 1

        self.divider_xtr_i = i
        with open(path_to_xtras) as f:
            for line in f:
                # j += 1
                name = line.strip()
                if name == '':
                    print('\n\nEmpty line in the middle of AEMS.dat, it will not work!\n\n')
                    exit('ARCA will not start, check out '+path_to_xtras+' and remove empty lines.')
                self.NAMES.append(name)
                namesPyInds[name] = i
                namesFoInds[name] = i+1
                i += 1

        ## -----------------------------------------------------------
        self.namesPyInds = namesPyInds
        self.namesFoInds = namesFoInds
