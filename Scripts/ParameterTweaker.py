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
fig, ax = plt.subplots(figsize=(8,7), num='ParameterTweaker')
plt.subplots_adjust(left=0.15, bottom=0.4)

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
axsigma = plt.axes([0.15, 0.27, 0.75, 0.02], facecolor=axcolor)
axpeak = plt.axes([0.15, 0.24, 0.75, 0.02], facecolor=axcolor)
axfreq = plt.axes([0.15, 0.21, 0.75, 0.02], facecolor=axcolor)
axphase = plt.axes([0.15, 0.18, 0.75, 0.02], facecolor=axcolor)
axamp = plt.axes([0.15, 0.15, 0.75, 0.02], facecolor=axcolor)

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

resetax = plt.axes([0.8, 0.025, 0.1, 0.03])
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

codeax = plt.axes([0.38, 0.025, 0.24, 0.03])
code = Button(codeax, 'Show code on terminal', color=axcolor, hovercolor='0.975')

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

code.on_clicked(ccode)

plt.show()
