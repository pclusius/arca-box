import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import string

# only relevant value you might occasionally have to change in the code
duration_in_hours = 24
#----------------------------------------------------------------------
#----------------------------------------------------------------------

indices = [
"CH3OH      = 20",
"C2H5OH     = 21",
"NPROPOL    = 22",
"IPROPOL    = 23",
"NBUTOL     = 24",
"BUT2OL     = 25",
"IBUTOL     = 26",
"TBUTOL     = 27",
"PECOH      = 28",
"IPEAOH     = 29",
"ME3BUOL    = 30",
"IPECOH     = 31",
"IPEBOH     = 32",
"CYHEXOL    = 33",
"MIBKAOH    = 34",
"ETHGLY     = 35",
"PROPGLY    = 36",
"MBO        = 37",
"HCHO       = 38",
"CH3CHO     = 39",
"C2H5CHO    = 40",
"C3H7CHO    = 41",
"IPRCHO     = 42",
"C4H9CHO    = 43",
"ACR        = 44",
"MACR       = 45",
"C4ALDB     = 46",
"CH4        = 47",
"C2H6       = 48",
"C3H8       = 49",
"NC4H10     = 50",
"IC4H10     = 51",
"NC5H12     = 52",
"IC5H12     = 53",
"NEOP       = 54",
"NC6H14     = 55",
"M2PE       = 56",
"M3PE       = 57",
"M22C4      = 58",
"M23C4      = 59",
"NC7H16     = 60",
"M2HEX      = 61",
"M3HEX      = 62",
"NC8H18     = 63",
"NC9H20     = 64",
"NC10H22    = 65",
"NC11H24    = 66",
"NC12H26    = 67",
"CHEX       = 68",
"C2H4       = 69",
"C3H6       = 70",
"BUT1ENE    = 71",
"CBUT2ENE   = 72",
"TBUT2ENE   = 73",
"MEPROPENE  = 74",
"PENT1ENE   = 75",
"CPENT2ENE  = 76",
"TPENT2ENE  = 77",
"ME2BUT1ENE = 78",
"ME3BUT1ENE = 79",
"ME2BUT2ENE = 80",
"HEX1ENE    = 81",
"CHEX2ENE   = 82",
"THEX2ENE   = 83",
"DM23BU2ENE = 84",
"C2H2       = 85",
"BENZENE    = 86",
"TOLUENE    = 87",
"OXYL       = 88",
"MXYL       = 89",
"PXYL       = 90",
"EBENZ      = 91",
"PBENZ      = 92",
"IPBENZ     = 93",
"TM123B     = 94",
"TM124B     = 95",
"TM135B     = 96",
"OETHTOL    = 97",
"METHTOL    = 98",
"PETHTOL    = 99",
"DIME35EB   = 100",
"DIET35TOL  = 101",
"STYRENE    = 102",
"BENZAL     = 103",
"CH3CL      = 104",
"CH2CL2     = 105",
"CHCL3      = 106",
"CH3CCL3    = 107",
"TCE        = 108",
"TRICLETH   = 109",
"CDICLETH   = 110",
"TDICLETH   = 111",
"CH2CLCH2CL = 112",
"CCL2CH2    = 113",
"CL12PROP   = 114",
"CHCL2CH3   = 115",
"CH3CH2CL   = 116",
"CHCL2CHCL2 = 117",
"CH2CLCHCL2 = 118",
"VINCL      = 119",
"C4H6       = 120",
"C5H8       = 121",
"CH3OCHO    = 122",
"METHACET   = 123",
"ETHACET    = 124",
"NPROACET   = 125",
"IPROACET   = 126",
"NBUTACET   = 127",
"SBUTACET   = 128",
"TBUACET    = 129",
"CH3OCH3    = 130",
"DIETETHER  = 131",
"MTBE       = 132",
"DIIPRETHER = 133",
"ETBE       = 134",
"MO2EOL     = 135",
"EOX2EOL    = 136",
"PR2OHMOX   = 137",
"BUOX2ETOH  = 138",
"BOX2PROL   = 139",
"CH3BR      = 140",
"DIBRET     = 141",
"CH3COCH3   = 142",
"MEK        = 143",
"MPRK       = 144",
"DIEK       = 145",
"MIPK       = 146",
"HEX2ONE    = 147",
"HEX3ONE    = 148",
"MIBK       = 149",
"MTBK       = 150",
"CYHEXONE   = 151",
"APINENE    = 152",
"BPINENE    = 153",
"LIMONENE   = 154",
"BCARY      = 155",
"HCOOH      = 156",
"CH3CO2H    = 157",
"PROPACID   = 158",
"DMM        = 159",
"DMC        = 160",
"DMS        = 161",
"ETHOX      = 162"]


x = np.linspace(0.0, duration_in_hours, 1000)
minimum_level = 1e1
maximum_level = 1e5
a0 = 1
f0 = 0.0
p0 = duration_in_hours/2
s0 = duration_in_hours/20
fs = 0

class u:
    def __init__(self):
        self.minimum = minimum_level
        self.maximum = maximum_level
range = u()
fig, ax = plt.subplots(figsize=(10,7), num='ParameterTweaker')
plt.subplots_adjust(left=0.15, bottom=0.4, right=0.7)

plt.title('Superbox function parameter tweaker')
ax.set_xlabel('time (hours)')
ax.set_ylabel('numbers')

D = p0 + np.sin((x-p0)*f0)*a0 + fs
f = 1/np.sqrt(2*np.pi*s0**2)
f = (range.maximum-range.minimum)/f
s = 1/np.sqrt(2*np.pi*s0**2)*np.exp(-(x-D)**2/(2*s0**2))
s = s*f + range.minimum

l, = plt.plot(x, s, lw=2, color='red')

plt.axis([0, duration_in_hours, range.minimum*0.95, range.maximum*1.05])
ax.set_yscale('linear')
ax.grid(True)

axcolor = 'lightgoldenrodyellow'
axsigma = plt.axes([0.15, 0.27, 0.55, 0.02], facecolor=axcolor)
axpeak = plt.axes([0.15, 0.24, 0.55, 0.02], facecolor=axcolor)
axfreq = plt.axes([0.15, 0.21, 0.55, 0.02], facecolor=axcolor)
axphase = plt.axes([0.15, 0.18, 0.55, 0.02], facecolor=axcolor)
axamp = plt.axes([0.15, 0.15, 0.55, 0.02], facecolor=axcolor)

ssigma = Slider(axsigma, 'Width', 0.01, duration_in_hours/2, valinit=s0, color='orange')
speak  = Slider(axpeak, 'Peaktime', 0.00, duration_in_hours, valinit=p0, color='orange')
sfreq  = Slider(axfreq, 'Ang.Freq', -2, 2, valinit=f0)
sphase = Slider(axphase, 'Phase', -duration_in_hours/2, duration_in_hours/2, valinit=fs)
samp   = Slider(axamp, 'Amplitude', 0.0, 10.0, valinit=a0)

def update(val):
    sigma = ssigma.val
    peak = speak.val
    amp = samp.val
    freq = sfreq.val
    phase = sphase.val
    f = 1/np.sqrt(2*np.pi*sigma**2)
    D = peak + np.sin((x-peak)*freq)*amp + phase
    norm = 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-D)**2/(2*sigma**2))
    if radio.value_selected == 'linear':
        f = (range.maximum-range.minimum)/f
        norm = norm*f + range.minimum
    else:
        f = (np.log10(range.maximum-range.minimum+1))/f
        norm = 10**(norm*f)-1 + range.minimum
    ax.grid(True)
    l.set_ydata(norm)
    fig.canvas.draw_idle()

sfreq.on_changed(update)
samp.on_changed(update)
speak.on_changed(update)
ssigma.on_changed(update)
sphase.on_changed(update)

resetax = plt.axes([0.2, 0.025, 0.1, 0.03])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    sfreq.reset()
    samp.reset()
    ssigma.reset()
    speak.reset()
    sphase.reset()
button.on_clicked(reset)

rax = plt.axes([0.15, 0.78, 0.1, 0.1], facecolor='w', alpha=0.1)
radio = RadioButtons(rax, ('linear','log'), active=0)

def scalefunc(scale):
    ax.set_yscale(scale)

radio.on_clicked(scalefunc)
radio.on_clicked(update)

non_numeric_chars = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-/:;<=>?@[\\]^_`{|}~ \t\n\r\x0b\x0c'
def pint(u):
    u=u.translate({ord(c): None for c in non_numeric_chars})
    range.minimum=float(u)
    ax.set_ylim(max(range.minimum*0.95, 0.01), range.maximum*1.05)
def pint2(u):
    u=u.translate({ord(c): None for c in non_numeric_chars})
    range.maximum=float(u)
    ax.set_ylim(max(range.minimum*0.95, 0.01), range.maximum*1.05)

tbax = plt.axes([0.15, 0.075, 0.15, 0.03], facecolor='none')
minb = TextBox(tbax, 'minimum', initial=str(range.minimum), color='.95', hovercolor='1', label_pad=0.01)
minb.on_text_change(pint)
minb.on_text_change(update)
tbax.annotate('Change minimum and maximum. Using TAB will cause freeze', wrap=True, xy=(-.8,0.1),
    xycoords=('axes fraction', 'figure fraction'), xytext=(0, 10), textcoords='offset points')

tb2ax = plt.axes([0.40, 0.075, 0.15, 0.03], facecolor='none')
maxb = TextBox(tb2ax, 'maximum', initial=str(range.maximum), color='.95', hovercolor='1', label_pad=0.01)
maxb.on_text_change(pint2)
maxb.on_text_change(update)

codeax = plt.axes([0.38, 0.025, 0.34, 0.03])
code = Button(codeax, 'Show fortran code on terminal', color=axcolor, hovercolor='0.975')

def ccode(event):
#    strr ='type(input_mod) :: [variable-name] = '
    strr = 'input_mod(min=%fd%d,max=%fd%d,sig=%fd0,mju=%0fd0,fv=%fd0,ph=%fd0,am=%fd0, MODE='%(
    range.minimum/(10**np.floor(np.log10(range.minimum))),np.floor(np.log10(range.minimum)),
    range.maximum/(10**np.floor(np.log10(range.maximum))),np.floor(np.log10(range.maximum)),
    ssigma.val,speak.val,
    sfreq.val,sphase.val,samp.val)
    if radio.value_selected == 'linear':
        strr = strr+'1)'
    else:
        strr = strr+'2)'
    print()
    print(strr)

def codeout(event,index,j):
    if radio.value_selected == 'linear':
        mode = 1
    else:
        mode = 2
    name = j
    strr = 'MODS(%d) = %d %fd%d %fd%d %fd0 %0fd0 %fd0 %fd0 %fd0 1d0 0d0 "%s"'%(
    index,mode,
    range.minimum/(10**np.floor(np.log10(range.minimum))),np.floor(np.log10(range.minimum)),
    range.maximum/(10**np.floor(np.log10(range.maximum))),np.floor(np.log10(range.maximum)),
    ssigma.val,speak.val,
    sfreq.val,sphase.val,samp.val, name)
    print()
    print(strr)

def c_temperature(event):
    codeout(event,1,'temperature')
def c_pressure(event):
    codeout(event,2,'pressure')
def c_relative_humid(event):
    codeout(event,3,'relative_humid')
def c_condens_sink(event):
    codeout(event,4,'condens_sink')
def c_shortwave_rad(event):
    codeout(event,5,'shortwave_rad')
def c_Ion_Prod_Rate(event):
    codeout(event,6,'Ion_Prod_Rate')
def c_H2SO4(event):
    codeout(event,7,'H2SO4')
def c_NH3(event):
    codeout(event,8,'NH3')
def c_DMA(event):
    codeout(event,9,'DMA')
def c_SO2(event):
    codeout(event,10,'SO2')
def c_NO(event):
    codeout(event,11,'NO')
def c_NO2(event):
    codeout(event,12,'NO2')
def c_CO(event):
    codeout(event,13,'CO')
def c_H2(event):
    codeout(event,14,'H2')
def c_O3(event):
    codeout(event,15,'O3')
def LIST(event):
    for c in indices:
        print(c)


compoundaxes =[]
j=16
tbax.annotate('Print input for init file:', wrap=True, xy=(0.8,0.025+0.03*(j+1)),
xycoords=('figure fraction', 'figure fraction'), xytext=(0, 10), textcoords='offset points')
while j>0:
    compoundaxes.append(plt.axes([0.78, 0.025+0.03*j, 0.20, 0.03]))
    j=j-1

i=-1
i=i+1;code1  = Button(compoundaxes[i], 'temperature', color=axcolor, hovercolor='0.975')
i=i+1;code2  = Button(compoundaxes[i], 'pressure', color=axcolor, hovercolor='0.975')
i=i+1;code3  = Button(compoundaxes[i], 'relative_humid', color=axcolor, hovercolor='0.975')
i=i+1;code4  = Button(compoundaxes[i], 'condens_sink', color=axcolor, hovercolor='0.975')
i=i+1;code5  = Button(compoundaxes[i], 'shortwave_rad', color=axcolor, hovercolor='0.975')
i=i+1;code6 = Button(compoundaxes[i], 'Ion_Prod_Rate', color=axcolor, hovercolor='0.975')
i=i+1;code7  = Button(compoundaxes[i], 'H2SO4', color=axcolor, hovercolor='0.975')
i=i+1;code8  = Button(compoundaxes[i], 'NH3', color=axcolor, hovercolor='0.975')
i=i+1;code9  = Button(compoundaxes[i], 'DMA', color=axcolor, hovercolor='0.975')
i=i+1;code10  = Button(compoundaxes[i], 'SO2', color=axcolor, hovercolor='0.975')
i=i+1;code11 = Button(compoundaxes[i], 'NO', color=axcolor, hovercolor='0.975')
i=i+1;code12 = Button(compoundaxes[i], 'NO2', color=axcolor, hovercolor='0.975')
i=i+1;code13 = Button(compoundaxes[i], 'CO', color=axcolor, hovercolor='0.975')
i=i+1;code14 = Button(compoundaxes[i], 'H2', color=axcolor, hovercolor='0.975')
i=i+1;code15 = Button(compoundaxes[i], 'O3', color=axcolor, hovercolor='0.975')
i=i+1;code16 = Button(compoundaxes[i], 'print list of VOC', color=axcolor, hovercolor='0.975')



code1.on_clicked(c_temperature)
code2.on_clicked(c_pressure)
code3.on_clicked(c_relative_humid)
code4.on_clicked(c_condens_sink)
code5.on_clicked(c_shortwave_rad)
code6.on_clicked(c_Ion_Prod_Rate)
code7.on_clicked(c_H2SO4)
code8.on_clicked(c_NH3)
code9.on_clicked(c_DMA)
code10.on_clicked(c_SO2)
code11.on_clicked(c_NO)
code12.on_clicked(c_NO2)
code13.on_clicked(c_CO)
code14.on_clicked(c_H2)
code15.on_clicked(c_O3)
code16.on_clicked(LIST)

#
code.on_clicked(ccode)

plt.show()
