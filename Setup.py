import tkinter as tk
from tkinter import ttk
import tkinter.filedialog as tkf
#from tkinter import filedialog
import tkinter.messagebox as tkm

import platform
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt


################################################################################
# Functions
################################################################################
class Comp:
    def __init__(self):
        self.index  = 0
        self.mode  = 0
        self.col  = -1
        self.multi = 1e0    # Multiplication factor in MODE0
        self.shift = 0e0    # Constant to be addded in MODE0
        self.min = 1e1      # Minimum value for the parametrized concentration OR constant value if max <= min
        self.max = 1e5      # Peak value
        self.sig = 1e0      # Standard deviation for the Gaussian=sig of the bell curve
        self.mju = 12e0     # Time of peak value
        self.fv  = 0e0      # Angular frequency [hours] of modifying sine function
        self.ph  = 0e0      # Angular frequency [hours] of modifying sine function
        self.am  = 1e0      # Amplitude of modificaion
        self.name = 'NONAME'# Human readable name for modified variable
        self.Find = 1

class indix:
    def __init__(self, F, P):
        self.Find = F
        self.Nind = ''
        self.savestatus = 0

def ind(name, lista):
    for i in range(len(lista)):
        if lista[i] == name:
            break
    return i

def tuhoa(name):
    i = ind(name,names) # indeksi freimien listassa
    print('Removing ', mods[i].name, ' from variables.')
    frames[i].destroy() # tuhotaan ko. freimi
    frames.pop(i)
    names.pop(i)
    mods.pop(i)
    ModsInUse.pop(i)
    ModsEntry.delete(i,None)
    SelectedEntry1.delete(i,None)
    for i,f in enumerate(frames):
        f.grid(row=i)

def add_variable_to_selected(x='x'):
    j=list_all_vars.curselection()[0]
    if j<lastenv:
        name=VARS[j]
    elif j>=lastenv:
        name=VARS[j+1]
    if name in names:
        print('- Note: ',name, 'is already included.')
    else:
        frames.append(tk.Frame(ScFrame))
        names.append(name)
        mods.append(Comp())
        ModsInUse.append(tk.IntVar())

        i = len(frames)-1
        frames[i].grid(row=i)
        mods[i].col = tk.IntVar()
        mods[i].col.set(-1)
        mods[i].multi = tk.IntVar()
        mods[i].multi.set(1)
        mods[i].shift = tk.IntVar()
        mods[i].shift.set(0)
        mods[i].name = name
        mods[i].Find = ind(name, VARS)+1
        print('Adding ', mods[i].name, ' to variables.')

        frames[i].grid(column=1, sticky=tk.W)
        tk.Label(frames[i], text='%s(%d)'%(name,mods[i].Find), width=20).grid(row=1, column=0, sticky=tk.W)
        tk.Checkbutton(frames[i], variable=ModsInUse[i], width=5).grid(row=1,column=1, sticky=tk.W)
        tk.Entry(frames[i], textvariable=mods[i].col, selectbackground=dark, width=8).grid(row=1, column=2,sticky=tk.W)
        tk.Entry(frames[i], textvariable=mods[i].multi, selectbackground=dark, width=8).grid(row=1, column=3,sticky=tk.W)
        tk.Entry(frames[i], textvariable=mods[i].shift, selectbackground=dark, width=8).grid(row=1, column=4,sticky=tk.W)
        tk.Button(frames[i], text='remove', command=lambda : tuhoa(name)).grid(row=1, column=5, sticky=tk.E)
        # adding variable also to tab1 and tab3
        ModsEntry.insert(tk.END, mods[i].name)
        SelectedEntry1.insert(tk.END, '%s(%s)'%(mods[i].name, mods[i].Find))

def onFrameConfigure(canvas):
    try:
        _onFrameConfigure(canvas)
    except:
        pass
def _onFrameConfigure(canvas):
    canvas.configure(scrollregion=canvas.bbox("all"))

def canv_width(event, canvas_frame):
    canvas_width = event.width
    canvas.itemconfig(canvas_frame, width = canvas_width)

def mouse_scroll(event, canvas):
    os_name = platform.system()
    if event.delta:
        if os_name == 'Darwin':
            canvas.yview_scroll(int(-1*(event.delta)), 'units') #MacOSX
        else:
            canvas.yview_scroll(int(-1*(event.delta/120)), 'units') #Windows, Linux
    else:
        if event.num == 5:
            move = 1
        else:
            move = -1
        canvas.yview_scroll(move, 'units')

def loadMod():
    compound = ModsEntry.curselection()
    short_i = compound[0]
    ax.set_title(names[short_i])
    if mods[short_i].mode>0:
        sigSc.set(mods[short_i].sig)
        peakSc.set(mods[short_i].mju)
        yscale.set(mods[short_i].mode)
        mini.set(mods[short_i].min)
        maxi.set(mods[short_i].max)
        AngWSc.set(mods[short_i].fv)
        PhaseSc.set(mods[short_i].ph)
        AmpSc.set(mods[short_i].am)
    tracker.Find = mods[short_i].Find
    tracker.Nind = names[short_i]

def saveModPara():
    if (tracker.Nind in names):
        i = ind(tracker.Nind, names)
        mods[i].mode  = yscale.get()
        mods[i].min   = mini.get()
        mods[i].max   = maxi.get()
        mods[i].sig   = sigSc.get()
        mods[i].mju   = peakSc.get()
        mods[i].fv    = AngWSc.get()
        mods[i].ph    = PhaseSc.get()
        mods[i].am    = AmpSc.get()
        ModsInUse[i].set(1)
    elif tracker.Nind != '':
        print('Variable %s was removed from list of variables.'%(tracker.Nind))
    else:
        pass

def updatenorm(val='x'):
    if len(names)>0:
        compound = ModsEntry.curselection()
        try:
            short_i = compound[0]
            ax.set_title(names[short_i])
        except IndexError:
            ax.set_title('Double click variables to select them for addition.')
        sigma = sigSc.get()#ssigma.val
        if abs(sigma)<0.01: sigma = 0.01
        peak = peakSc.get()#speak.val
        try:
            mini = float(minEntry.get())
        except:
            mini = 0
        try:
            maxi = float(maxEntry.get())
        except:
            maxi=0
        freq = AngWSc.get()
        phase = PhaseSc.get()
        amp = AmpSc.get()
        f = 1/np.sqrt(2*np.pi*sigma**2)
        D = peak + np.sin((x-peak)*freq)*amp + phase
        norm = 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-D)**2/(2*sigma**2))

        if yscale.get() == 1:
            f = (maxi-mini)/f
            norm = norm*f + mini
        else:
            f = (np.log10(maxi-mini+1))/f
            norm = 10**(norm*f)-1 + mini
        ax.grid(True)
        l.set_ydata(norm)
        if yscale.get() == 1:
            ax.set_ylim(max(mini*0.95, 0.01), maxi*1.05)
            ax.set_yscale('linear')
        if yscale.get() == 2:
            ax.set_ylim(max(mini*0.95, 0.01), maxi*1.5)
            ax.set_yscale('log')
        fig.canvas.draw_idle()
    else:
        pass

def check_changes():
    if tracker.Nind in names:
        i = ind(tracker.Nind, names)
        if mods[i].mode    != yscale.get() and mods[i].mode!= 0: print('mode');return 1
        elif mods[i].min   != mini.get(): print('min');return 1
        elif mods[i].max   != maxi.get(): print('max',mods[i].max, maxi.get());return 1
        elif mods[i].sig   != sigSc.get(): print('sig', mods[i].sig,sigSc.get());return 1
        elif mods[i].mju   != peakSc.get(): print('peak');return 1
        elif mods[i].fv    != AngWSc.get(): print('fv');return 1
        elif mods[i].ph    != PhaseSc.get(): print('ph');return 1
        elif mods[i].am    != AmpSc.get(): print('am');return 1
        else: return 0
    else: return 0

def focus_on_variable(x='x'):
    if check_changes() == 0:
        proceed = 'yes'
    else:
        proceed = tkm.askquestion("Save Changes?","Settings for %s have changed, save current settings?"%(tracker.Nind), icon='question')
    if proceed == 'yes':
        saveModPara()
    loadMod()
    updatenorm()

def printInit():
    for i,m in enumerate(mods):
        if m.col.get() == 1:
            tkm.showinfo("Stop right there","'Column' (used in %s) can not be 1, it is reserved for time."%(m.name))
            return
        print(m.col.get(),m.mode)
        if ModsInUse[i].get() == 0 and m.shift.get() == 0 and m.col.get() == -1 and m.mode>0:
            tkm.showinfo("Note","You are explicitly sending no input for %s."%(m.name))
            return
    for i,m in enumerate(mods):
        multistr = '%12.5e'%(m.multi.get())
        multistr = multistr.replace('e', 'd', 1)
        shiftstr = '%12.5e'%(m.shift.get())
        shiftstr = shiftstr.replace('e', 'd', 1)
        minstr = '%12.5e'%(m.min)
        minstr = minstr.replace('e', 'd', 1)
        maxstr = '%12.5e'%(m.max)
        maxstr = maxstr.replace('e', 'd', 1)
        if ModsInUse[i].get() == 1: mode = m.mode
        else: mode = 0
        strr = 'MODS(%d) = %d %d %s %s %s %s %fd0 %0fd0 %fd0 %fd0 %fd0 1d0 0d0 ! %s'%(
        m.Find,m.col.get(), mode, multistr,shiftstr,minstr, maxstr, m.sig,m.mju, m.fv,m.ph,m.am, m.name)
        print(strr)

def prnt():
    print('# Currently not implimented')

def SaveDef(name):
    pass

def saveInit():
    print('# Currently not quite implimented')
    # This prompts for file name and locaction
    name=tkf.asksaveasfile(initialdir = "../input/", mode='w',defaultextension="")

    if name != None:
        # Idea of this if is that if user chooses CANCEL, the function won't be upset about it
        # since name.name will be empty
        SaveDef(name.name)
        name.close()
    else: return


def Stop():
    plt.close()
    canvas.destroy()
    root.destroy()

###############################################################################
# End Functions
###############################################################################
# Find the system type (Windows, Linux, Darwin (MacOS)) for system specific configuration
os_name = platform.system()

rt = input('Please give runtime in hours (empty: 24 hours): ')

if rt == float():
    rt = float(rt)
else:
    rt = 24.0

# try:
#     rt = float(rt)
# except:
#     rt = 24.0

print('Model runtime: ', rt, ' hours')


root = tk.Tk()
root.geometry("1024x850")
root.title('Superbox configurator')

dark = '#666666'

# Read NAMES.DAT to find out the available variables
VARS = np.genfromtxt('NAMES.dat', dtype=str, usecols=0, comments='!')
for i,v in enumerate(VARS):
    if v == '#':
        lastenv = i
        break
env = VARS[0:lastenv]
mcm = VARS[lastenv+1:]

# Create frame for common buttons
toppen = tk.Frame(root)
toppen.pack(side=tk.TOP)
titleLabel = tk.Label(toppen, text='GAP Setup app',font=("Helvetica", 16),anchor=tk.W)
titleLabel.grid(row=i, padx=1, pady=1,ipadx=2, ipady=2, column=1, columnspan=4, sticky=tk.E+tk.W)


base = tk.Frame(root)
base.pack(side=tk.BOTTOM)
tk.Button(base, text='Print MODS only', command=printInit).pack(side=tk.LEFT,fill=tk.X)
tk.Button(base, text="Print INITFILE", command=prnt).pack(side=tk.LEFT,fill=tk.X)
tk.Button(base, text="Save INITFILE", command=saveInit).pack(side=tk.LEFT,fill=tk.X)
ttk.Separator(base, orient=tk.VERTICAL).pack(side=tk.LEFT,ipadx=10)
tk.Button(base, text='Exit without save', command=Stop).pack(side=tk.LEFT,fill=tk.X,padx=40)

# Define tabs
tabControl = ttk.Notebook(root)         # Create Tab Control
tab1 = ttk.Frame(tabControl)            # Create a tab
tab2 = ttk.Frame(tabControl)            # Create a tab
tab3 = ttk.Frame(tabControl)            # Create a tab
tab1 = ttk.Frame(tabControl)            # Create a tab
tab2 = ttk.Frame(tabControl)            # Create a tab
tab3 = ttk.Frame(tabControl)            # Create a tab
tabControl.add(tab1, text='General')  # Add the tab
tabControl.add(tab2, text='Input')    # Add the tab
tabControl.add(tab3, text='Parametrization')    # Add the tab
tabControl.pack(expand=1, fill="both")  # Pack to make visible

#################################################################################
# Tab 1
#################################################################################

# Variables primarily used in Tab 1
runtime  =    tk.StringVar();
PrInt  =    tk.StringVar();
FSInt  =    tk.StringVar();
case  =    tk.StringVar();
run   =    tk.DoubleVar();
ENVPATH  =    tk.StringVar();
MCMPATH  =    tk.StringVar();
PARPATH  =    tk.StringVar();
TempUnit  = tk.IntVar();
Aerosol_flag    = tk.IntVar();
Chemistry_flag  = tk.IntVar()
NUCLEATION      = tk.IntVar()
ACDC            = tk.IntVar()
Aerosol_flag    = tk.IntVar()
Particle_flag   = tk.IntVar()
Extra_data      = tk.IntVar()
Current_case    = tk.IntVar()

TempUnit.set(0)
Aerosol_flag.set(1)
Chemistry_flag.set(1)
Particle_flag.set(1)
Extra_data.set(1)
Current_case.set(1)
run.set(str(rt))

Chemistry_flag.set(1)
NUCLEATION.set(1)
ACDC.set(1)
Aerosol_flag.set(1)


# Layout the GUI elements
top=tk.Frame(tab1, height=12,width=850)
top.grid(row=0)
f1=tk.Frame(tab1, height=838-80,width=850)
f1.grid(row=1)
f1.grid_propagate(0)


frameL1 = tk.LabelFrame(f1, text=" General options ", relief=tk.RIDGE)
frameL1.grid(row=1, column=1, sticky=tk.E + tk.W + tk.N + tk.S, padx=30, pady=10,ipadx=10, ipady=10)

i = 0
tk.Label(frameL1, text='Model run time in hours: ').grid(row=i,column=1, sticky=tk.EW)
runEntry = tk.Entry(frameL1, selectbackground=dark, textvariable=run)
runEntry.grid(row=i, column=2)

i += 1
tk.Label(frameL1, text='YEAR yyyy: ').grid(row=i, column=1, sticky=tk.EW)
runtimeEntry = tk.Entry(frameL1, selectbackground=dark, textvariable=runtime)
runtimeEntry.grid(row=i, column=2)

i += 1
tk.Label(frameL1, text='PrInt mmdd: ').grid(row=i, column=1, sticky=tk.EW)
PrIntEntry = tk.Entry(frameL1, selectbackground=dark, textvariable=PrInt)
PrIntEntry.grid(row=i, column=2)

i += 1
tk.Label(frameL1, text='CASE [XX]: ').grid(row=i, column=1, sticky=tk.EW)
caseEntry = tk.Entry(frameL1, selectbackground=dark, textvariable=case)
caseEntry.grid(row=i, column=2)

i += 1
tk.Label(frameL1,  text="Temperature in Celsius: ").grid(row=i, column=1, sticky=tk.EW)
StRadBox1 = tk.Checkbutton(frameL1, variable=TempUnit).grid(row=i,column=2, sticky=tk.W)


frameL2 = tk.LabelFrame(f1, text=" Modules in use ", relief=tk.RIDGE)
frameL2.grid(row=3, column=1, sticky=tk.E + tk.W + tk.N + tk.S, padx=30, pady=10,ipadx=10, ipady=10)

i =0
tk.Label(frameL2,  text="Use chemistry: ").grid(row=i, column=1, sticky=tk.EW)
StRadBox2 = tk.Checkbutton(frameL2, variable=Chemistry_flag).grid(row=i,column=2, sticky=tk.W)

i += 1
tk.Label(frameL2,  text="Use nucleation: ").grid(row=i, column=1, sticky=tk.EW)
StRadBox2    = tk.Checkbutton(frameL2, variable=NUCLEATION).grid(row=i,column=2, sticky=tk.W)

i += 1
tk.Label(frameL2,  text="Use ACDC: ").grid(row=i, column=1, sticky=tk.EW)
StRadBox2 = tk.Checkbutton(frameL2, variable=ACDC).grid(row=i,column=2, sticky=tk.W)

i += 1
tk.Label(frameL2,  text="Use aerosol: ").grid(row=i, column=1, sticky=tk.EW)
StRadBox2 = tk.Checkbutton(frameL2, variable=Aerosol_flag).grid(row=i,column=2, sticky=tk.W)

frameR2 = tk.LabelFrame(f1, text=" Variables in use ", relief=tk.RIDGE)
frameR2.grid(row=1, column=2, sticky=tk.E + tk.W + tk.N + tk.S,rowspan=30, padx=30, pady=10,ipadx=10, ipady=10)

SelectedEntry1 = tk.Listbox(frameR2, selectbackground=dark, selectmode=tk.BROWSE, height=30, bd=0, relief=tk.FLAT)
SelectedEntry1.grid(row=1, column=0, padx=10, sticky=tk.E+tk.W)#;i+=1


frameL3 = tk.LabelFrame(f1, text=" Paths and files ", relief=tk.RIDGE, width=450, height=200)
frameL3.grid(row=4, column=1, sticky=tk.E + tk.W + tk.N + tk.S, padx=30, pady=10,ipadx=10, ipady=10)
frameL3.grid_propagate(0)
frameL3.grid_columnconfigure(1,weight=3)

def browse(f):
    filename = tkf.askdirectory()
    if filename != '':
        f.set(filename)

tk.Label(frameL3,text='ENV FILE', width=8).grid(row=0, column=0)
tk.Entry(frameL3,textvariable=ENVPATH, selectbackground=dark).grid(row=0, column=1, sticky=tk.EW)
tk.Button(frameL3,text="Browse for ENV", command=lambda f=ENVPATH:browse(f), width=13).grid(row=0, column=2)

tk.Label(frameL3,text='MCM FILE', width=8).grid(row=1, column=0)
tk.Entry(frameL3,textvariable=MCMPATH, selectbackground=dark).grid(row=1, column=1, sticky=tk.EW)
tk.Button(frameL3,text="Browse for MCM", command=lambda f=MCMPATH:browse(f), width=13).grid(row=1, column=2)

tk.Label(frameL3,text='PAR FILE', width=8).grid(row=2, column=0)
tk.Entry(frameL3,textvariable=PARPATH, selectbackground=dark).grid(row=2, column=1, sticky=tk.EW)
tk.Button(frameL3,text="Browse for PAR", command=lambda f=PARPATH:browse(f), width=13).grid(row=2, column=2)

###############################################################################################################
# tab 2
###############################################################################################################

frames = []
names  = []
mods  = []
ModsInUse  = []

tk.Frame(tab2, height=12,width=850).grid(row=0)
f2=tk.Frame(tab2, height=838-80,width=850)
f2.grid(row=1)
f2.grid_propagate(0)

tk.Label(f2, text='Key + Click the variables to pick them for the model:',font=("Helvetica", 14)).grid(row=1, column=1,  columnspan=3,sticky=tk.W,padx=20,pady=20)
ttk.Separator(f2, orient=tk.HORIZONTAL).grid(column=1, row=2, columnspan=3)

header = tk.Frame(f2)
header.grid(row=2, column=2, sticky=tk.N+tk.W)
tk.Label(header, text='Var. name (IND)',width=14).grid(row=1, column=0, sticky=tk.W,padx=20)
tk.Label(header, text='PM\n in use',width=6).grid(row=1, column=1, sticky=tk.W)
tk.Label(header, text='Column',width=8).grid(row=1, column=2, sticky=tk.W)
tk.Label(header, text='Multiply',width=9).grid(row=1, column=3, sticky=tk.W)
tk.Label(header, text='Shift by',width=9).grid(row=1, column=4, sticky=tk.W)

frameLeft = tk.Frame(f2)
frameLeft.grid(row=2,column=1, padx=20,sticky=tk.N, rowspan=30)
list_all_vars = tk.Listbox(frameLeft, selectmode=tk.BROWSE,selectbackground=dark, height=25)
list_all_vars.grid(row=1, column=1, sticky=tk.N+tk.S)

for c in VARS:
    if c != '#':
        list_all_vars.insert(tk.END, c)

if os_name == 'Darwin':
    # Button Release is important, prevents clicking related issues
    # "Command" is exclusive for Mac keyboards.
    list_all_vars.bind("<Command-ButtonRelease-1>", add_variable_to_selected)
else:
    # Windows/Linux use Shift instead of Command, it doesn't matter as long as it is set to ButtonRelease
    list_all_vars.bind("<Shift-ButtonRelease-1>", add_variable_to_selected)

frameRight = tk.Frame(f2)
frameRight.grid(row=3, column=2, sticky=tk.N, rowspan=25)

######
canvas = tk.Canvas(frameRight, width=600, height=450)
ScFrame = tk.Frame(canvas)
vertscroll = tk.Scrollbar(canvas, orient='vertical', command=canvas.yview)
canvas.configure(yscrollcommand=vertscroll.set)

txtvar = 'Variable(IND) = Variable in the model. IND refers to input index, handy if you manually modify INITFILE.'
txtpm = 'PM in use = Checked when the parametrization (PM) is first saved. Uncheck to override by column, check to use PM.'
txtcol = 'Column = Column number of whichever file the variable will be read in. CANNOT BE 1, it is reserved for time.'
txtmul = 'Multipl = Multiply values read from Column by this. If no input, will have no effect.'
txtshift = 'Shift = Shift values read from Column by this. Unless PM is used, will always add (or subtract) a constant.'

framebt = tk.Frame(f2)
framebt.grid(row=38,column=1, padx=3,sticky=tk.N, columnspan=3)
tk.Label(framebt, text=txtvar).pack(pady=3)
tk.Label(framebt, text=txtpm).pack(pady=3)
tk.Label(framebt, text=txtcol).pack(pady=3)
tk.Label(framebt, text=txtmul).pack(pady=3)
tk.Label(framebt, text=txtshift).pack(pady=3)

if os_name == 'Darwin':
    # on MacOSX root.bind must use onFrameConfigure() otherwise it creates TclError in the background
    root.bind('<Configure>', lambda event, canvas=canvas: onFrameConfigure(canvas))
    root.bind_all('<MouseWheel>', lambda event, canvas=canvas: mouse_scroll(event, canvas))
else:
    # on MacOSX these are not needed since you can scroll with mousewheel/trackpad
    # this part is needed for Linux/Windows
    root.bind_all('<Button-4>', lambda event, canvas=canvas: mouse_scroll(event, canvas))
    root.bind_all('<Button-5>', lambda event, canvas=canvas: mouse_scroll(event, canvas))
    # on MacOSX this part hangs tab2, adds vertical scroll bar for Windows/Linux
    # comment it if you don't want a fancy vertical scroll bar
    #vertscroll.pack(side=tk.RIGHT, fill=tk.Y)
    canvas.bind('<Configure>', lambda event, canvas_frame=canvas_frame: canv_width(event, canvas_frame))

canvas.grid(row=4, column=3, sticky=tk.N)
canvas_frame = canvas.create_window((4, 4), window=ScFrame, anchor="sw")

###############################################################################################################
# tab 3
###############################################################################################################
tk.Frame(tab3, height=12,width=1024).grid(row=0)
f3=tk.Frame(tab3, height=850,width=1024)
f3.grid(row=1)
f3.grid_propagate(0)

plt.style.use('ggplot')

tracker = indix(0,-1)
mod_defaults = Comp()

i=1
mini = tk.DoubleVar();
maxi = tk.DoubleVar();
mini.set(mod_defaults.min)
maxi.set(mod_defaults.max)
mini.trace("w", lambda name, index, mode, sv=mini: updatenorm(mini))
maxi.trace("w", lambda name, index, mode, sv=maxi: updatenorm(maxi))

yscale = tk.IntVar()
yscale.set(1)
scaleLabel = tk.Label(f3, text='Y-axis scale')
scaleLabel.grid(row=i,  sticky=tk.EW,padx=1, pady=1,ipadx=2, ipady=2, column=1)#;i+=1
loglin = tk.Radiobutton(f3, text="lin", variable=yscale, value=1, command=updatenorm).grid(row=i,  sticky=tk.EW,padx=1, pady=1,ipadx=2, ipady=2, column=2)
loglin = tk.Radiobutton(f3, text="log", variable=yscale, value=2, command=updatenorm).grid(row=i,  sticky=tk.EW,padx=1, pady=1,ipadx=2, ipady=2, column=3)

i+=1

fig = plt.figure(figsize=(3, 4.2), dpi=100)
ax = fig.add_subplot(111)
ax.set_title('First select input variables from "Input"-tab')
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95)

canvas3 = FigureCanvasTkAgg(fig, master=f3)
canvas3.get_tk_widget().grid(row=i, column=1, sticky=tk.EW, padx=1, pady=1,ipadx=2, ipady=2, columnspan=7);i+=1

tk.Label(f3, text='Width').grid(row=i, sticky=tk.E+tk.W+tk.N+tk.S, padx=1, pady=1,ipadx=2, ipady=2, column=1)
sigSc = tk.Scale(f3, command=updatenorm,from_=0.01, to=rt/2,resolution=0.01, orient=tk.HORIZONTAL)
sigSc.grid(row=i, column=2, sticky=tk.EW, padx=1, pady=1,ipadx=2, ipady=2, columnspan=6);i+=1
sigSc.set(mod_defaults.sig)

tk.Label(f3, text='Peaktime').grid(row=i, sticky=tk.E+tk.W+tk.N+tk.S, padx=1, pady=1,ipadx=2, ipady=2, column=1)
peakSc = tk.Scale(f3, command=updatenorm, from_=0, to=rt*1.2, resolution=0.01,orient=tk.HORIZONTAL)
peakSc.grid(row=i, column=2, sticky=tk.EW, padx=1, pady=1,ipadx=2, ipady=2, columnspan=6);i+=1
peakSc.set(mod_defaults.mju)

tk.Label(f3, text='Ang.Freq').grid(row=i, sticky=tk.E+tk.W+tk.N+tk.S, padx=1, pady=1,ipadx=2, ipady=2, column=1)
AngWSc = tk.Scale(f3, command=updatenorm, from_=-2, to=20, resolution=0.01, orient=tk.HORIZONTAL)
AngWSc.grid(row=i, column=2, sticky=tk.EW, padx=1, pady=1,ipadx=2, ipady=2, columnspan=6);i+=1
AngWSc.set(mod_defaults.fv)

tk.Label(f3, text='Phase').grid(row=i, sticky=tk.E+tk.W+tk.N+tk.S, padx=1, pady=1,ipadx=2, ipady=2, column=1)
PhaseSc = tk.Scale(f3, command=updatenorm, from_=-rt/0.2, to=rt/0.2,resolution=0.01, orient=tk.HORIZONTAL)
PhaseSc.grid(row=i, column=2, sticky=tk.E+tk.W, padx=1, pady=1,ipadx=2, ipady=2, columnspan=6);i+=1
PhaseSc.set(mod_defaults.ph)

tk.Label(f3, text='Amplitude').grid(row=i, sticky=tk.E+tk.W+tk.N+tk.S, padx=1, pady=1,ipadx=2, ipady=2, column=1)
AmpSc = tk.Scale(f3, command=updatenorm, from_=0, to=10,resolution=0.01, orient=tk.HORIZONTAL)
AmpSc.grid(row=i, column=2, sticky=tk.EW, padx=1, pady=1,ipadx=2, ipady=2, columnspan=6);i+=1
AmpSc.set(mod_defaults.am)

# create empty plot
x=np.linspace(0,rt,200)
l, = ax.plot(maxi.get(), mini.get(), lw=2, color='red')
l.set_xdata(x)
ax.axis([0, rt, 0,maxi.get()*1.05])#, range.minimum*0.95, range.maximum*1.05])
ax.set_yscale('linear')
ax.grid(True)
ax.set_xlabel('Time (h)', fontsize=11)
ax.set_ylabel('Numbers', fontsize=11)

scrollbaru = tk.Scrollbar(f3)
scrollbaru.grid(row=2, padx=1, pady=1,ipadx=2, ipady=2, column=9, sticky=tk.N+tk.S)

tk.Label(f3, text='Select variable for setup \nby double clicking below').grid(row=1, sticky=tk.E+tk.W+tk.N+tk.S, padx=1, pady=1,ipadx=2, ipady=2, column=8)#;i+=1
ModsEntry = tk.Listbox(f3,selectbackground=dark, selectmode=tk.BROWSE, yscrollcommand=scrollbaru.set,)
ModsEntry.grid(row=2, padx=1, pady=1,ipadx=2, ipady=2, column=8, rowspan=1, sticky=tk.N+tk.S+tk.E+tk.W)
scrollbaru.config(command=ModsEntry.yview)

SelectedBut = tk.Button(f3, text="Save values for\nselected variable", command=saveModPara)
SelectedBut.grid(row=3, sticky=tk.E+tk.W+tk.N+tk.S, padx=1, pady=1,ipadx=2, ipady=2, column=8)

tk.Label(f3, text='Minimum').grid(row=1, padx=1, pady=1,ipadx=2, ipady=2, column=4)#;i+=1
minEntry = tk.Entry(f3,selectbackground=dark, textvariable=mini, validate="focusout",  validatecommand=updatenorm)
minEntry.grid(row=1, padx=1, pady=1,ipadx=2, ipady=2, column=5)

tk.Label(f3, text='Maximum').grid(row=1, padx=1, pady=1,ipadx=2, ipady=2, column=6)#;i+=1
maxEntry = tk.Entry(f3,selectbackground=dark, textvariable=maxi, validate="focusout", validatecommand=updatenorm)
maxEntry.grid(row=1, padx=1, pady=1,ipadx=2, ipady=2, column=7)#;i+=1

if os_name == 'Darwin':
    # Button Release is important, prevents clicking related issues
    # "Command" is exclusive for Mac keyboards.
    ModsEntry.bind("<Command-ButtonRelease-1>", focus_on_variable)
else:
    # Windows/Linux use Shift instead of Command, it doesn't matter as long as it is set to ButtonRelease
    ModsEntry.bind("<Shift-ButtonRelease-1>", focus_on_variable)

tk.mainloop()
