import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import string

# only relevant value you might occasionally have to change in the code
duration_in_hours = 24
#----------------------------------------------------------------------
#----------------------------------------------------------------------

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
    strr = 'MODS(%d) = %d "%s" %fd%d %fd%d %fd0 %0fd0 %fd0 %fd0 %fd0 1d0 0d0'%(
    index,mode, name,
    range.minimum/(10**np.floor(np.log10(range.minimum))),np.floor(np.log10(range.minimum)),
    range.maximum/(10**np.floor(np.log10(range.maximum))),np.floor(np.log10(range.maximum)),
    ssigma.val,speak.val,
    sfreq.val,sphase.val,samp.val)
    print()
    print(strr)

def c_H2SO4(event):
    codeout(event,1,'H2SO4')
def c_NH3(event):
    codeout(event,2,'NH3')
def c_DMA(event):
    codeout(event,3,'DMA')
def c_condens_sink(event):
    codeout(event,4,'condens_sink')
def c_shortwave_rad(event):
    codeout(event,5,'shortwave_rad')
def c_relative_humid(event):
    codeout(event,6,'relative_humid')
def c_pressure(event):
    codeout(event,7,'pressure')
def c_temperature(event):
    codeout(event,8,'temperature')
def c_SO2(event):
    codeout(event,9,'SO2')
def c_NO(event):
    codeout(event,10,'NO')
def c_NO2(event):
    codeout(event,11,'NO2')
def c_CO(event):
    codeout(event,12,'CO')
def c_H2(event):
    codeout(event,13,'H2')
def c_O3(event):
    codeout(event,14,'O3')
def c_Ion_Prod_Rate(event):
    codeout(event,15,'Ion_Prod_Rate')
def c_CH3O(event):
    codeout(event,16,'CH3O')
def c_CH3C(event):
    codeout(event,17,'CH3C')
def c_C2H5(event):
    codeout(event,18,'C2H5')
def c_C5H8(event):
    codeout(event,19,'C5H8')
def c_MVK(event):
    codeout(event,20,'MVK')
def c_MEK(event):
    codeout(event,21,'MEK')
def c_BENZENE(event):
    codeout(event,22,'BENZENE')
def c_ALPHAPINENE(event):
    codeout(event,23,'ALPHAPINENE')
def c_BETAPINENE(event):
    codeout(event,24,'BETAPINENE')
def c_LIMONENE(event):
    codeout(event,25,'LIMONENE')
def c_CARENE(event):
    codeout(event,26,'CARENE')
def c_TOLUENE(event):
    codeout(event,27,'TOLUENE')


compoundaxes =[]
j=27
tbax.annotate('Print input for init file:', wrap=True, xy=(0.8,0.025+0.03*(j+1)),
xycoords=('figure fraction', 'figure fraction'), xytext=(0, 10), textcoords='offset points')
while j>0:
    compoundaxes.append(plt.axes([0.78, 0.025+0.03*j, 0.20, 0.03]))
    j=j-1


i=0  ;code1  = Button(compoundaxes[i], 'H2SO4', color=axcolor, hovercolor='0.975')
i=i+1;code2  = Button(compoundaxes[i], 'NH3', color=axcolor, hovercolor='0.975')
i=i+1;code3  = Button(compoundaxes[i], 'DMA', color=axcolor, hovercolor='0.975')
i=i+1;code4  = Button(compoundaxes[i], 'condens_sink', color=axcolor, hovercolor='0.975')
i=i+1;code5  = Button(compoundaxes[i], 'shortwave_rad', color=axcolor, hovercolor='0.975')
i=i+1;code6  = Button(compoundaxes[i], 'relative_humid', color=axcolor, hovercolor='0.975')
i=i+1;code7  = Button(compoundaxes[i], 'pressure', color=axcolor, hovercolor='0.975')
i=i+1;code8  = Button(compoundaxes[i], 'temperature', color=axcolor, hovercolor='0.975')
i=i+1;code9  = Button(compoundaxes[i], 'SO2', color=axcolor, hovercolor='0.975')
i=i+1;code10 = Button(compoundaxes[i], 'NO', color=axcolor, hovercolor='0.975')
i=i+1;code11 = Button(compoundaxes[i], 'NO2', color=axcolor, hovercolor='0.975')
i=i+1;code12 = Button(compoundaxes[i], 'CO', color=axcolor, hovercolor='0.975')
i=i+1;code13 = Button(compoundaxes[i], 'H2', color=axcolor, hovercolor='0.975')
i=i+1;code14 = Button(compoundaxes[i], 'O3', color=axcolor, hovercolor='0.975')
i=i+1;code15 = Button(compoundaxes[i], 'Ion_Prod_Rate', color=axcolor, hovercolor='0.975')
i=i+1;code16 = Button(compoundaxes[i], 'CH3O', color=axcolor, hovercolor='0.975')
i=i+1;code17 = Button(compoundaxes[i], 'CH3C', color=axcolor, hovercolor='0.975')
i=i+1;code18 = Button(compoundaxes[i], 'C2H5', color=axcolor, hovercolor='0.975')
i=i+1;code19 = Button(compoundaxes[i], 'C5H8', color=axcolor, hovercolor='0.975')
i=i+1;code20 = Button(compoundaxes[i], 'MVK', color=axcolor, hovercolor='0.975')
i=i+1;code21 = Button(compoundaxes[i], 'MEK', color=axcolor, hovercolor='0.975')
i=i+1;code22 = Button(compoundaxes[i], 'BENZENE', color=axcolor, hovercolor='0.975')
i=i+1;code23 = Button(compoundaxes[i], 'ALPHAPINENE', color=axcolor, hovercolor='0.975')
i=i+1;code24 = Button(compoundaxes[i], 'BETAPINENE', color=axcolor, hovercolor='0.975')
i=i+1;code25 = Button(compoundaxes[i], 'LIMONENE', color=axcolor, hovercolor='0.975')
i=i+1;code26 = Button(compoundaxes[i], 'CARENE', color=axcolor, hovercolor='0.975')
i=i+1;code27 = Button(compoundaxes[i], 'TOLUENE', color=axcolor, hovercolor='0.975')


code1.on_clicked(c_H2SO4)
code2.on_clicked(c_NH3)
code3.on_clicked(c_DMA)
code4.on_clicked(c_condens_sink)
code5.on_clicked(c_shortwave_rad)
code6.on_clicked(c_relative_humid)
code7.on_clicked(c_pressure)
code8.on_clicked(c_temperature)
code9.on_clicked(c_SO2)
code10.on_clicked(c_NO)
code11.on_clicked(c_NO2)
code12.on_clicked(c_CO)
code13.on_clicked(c_H2)
code14.on_clicked(c_O3)
code15.on_clicked(c_Ion_Prod_Rate)
code16.on_clicked(c_CH3O)
code17.on_clicked(c_CH3C)
code18.on_clicked(c_C2H5)
code19.on_clicked(c_C5H8)
code20.on_clicked(c_MVK)
code21.on_clicked(c_MEK)
code22.on_clicked(c_BENZENE)
code23.on_clicked(c_ALPHAPINENE)
code24.on_clicked(c_BETAPINENE)
code25.on_clicked(c_LIMONENE)
code26.on_clicked(c_CARENE)
code27.on_clicked(c_TOLUENE)

#
code.on_clicked(ccode)

plt.show()
