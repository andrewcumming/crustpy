from crust import *
import numpy
import time

from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput, Button
from bokeh.plotting import figure

# Calculate the outburst and plot the lightcurve
def do_outburst():
	# Initialize the crust
	crust = Crust(mass=params['M'],radius=params['R'],ngrid=50,Qimp=params['Qimp'],Tc=params['Tc'])
	print(crust)
	# Set the top temperature and accrete
	crust.set_top_temperature(params['Tb'])
	t, Teff = crust.evolve(time=params['tout']*365.0,mdot=params['mdot'])
	# Now turn off accretion and cool
	t, Teff = crust.evolve(time=10000,mdot=0.0)
	print('Done')
	return t, Teff
	# Final temperature profile to plot
	#rho2, TT2 = crust.temperature_profile()

# default parameters
def set_defaults():
	return dict(Qimp = 3.2, Tb = 4.2e8, Tc = 3.1e7, M=1.62, R=11.2, mdot =0.1, tout=2.5)

# Handlers for the different controls
def Q_slider_handler(attrname, old, new):
	params['Qimp'] = Q_slider.value

def Tb_slider_handler(attrname, old, new):
	params['Tb'] = Tb_slider.value * 1e8

def Tc_slider_handler(attrname, old, new):
	params['Tc'] = Tc_slider.value * 1e7

def R_slider_handler(attrname, old, new):
	params['R'] = R_slider.value

def M_slider_handler(attrname, old, new):
	params['M'] = M_slider.value

def mdot_slider_handler(attrname, old, new):
	params['mdot'] = mdot_slider.value

def tout_slider_handler(attrname, old, new):
	params['tout'] = tout_slider.value
	
def button_handler():
	go_button.label='Running'
	go_button.button_type = 'default'
	# recalculate the cooling curve
	t, Teff = do_outburst()
	source.data = dict(x=t, y=Teff)
	go_button.label='Go'
	go_button.button_type = 'success'

def def_handler():
	# reset to default parameters
	params = set_defaults()
	Q_slider.value = params['Qimp']
	Tc_slider.value = params['Tc']/1e7
	Tb_slider.value = params['Tb']/1e8
	M_slider.value = params['M']
	R_slider.value = params['R']
	tout_slider.value = params['tout']
	mdot_slider.value = params['mdot']

# Setup:
# These are the 1659 temperature measurements
t0 = 52159.5
tobs = numpy.array([52197.8,52563.2,52712.2,52768.9,53560.0,53576.7,54583.8,56113])-t0
Teffobs = numpy.array([121,85,77,73,58,54,56,48.8])

# Initial values of parameters
params = set_defaults()

# Calculate the cooling curve and make a plot
t, Teff = do_outburst()
source = ColumnDataSource(data=dict(x=t, y=Teff))
plot = figure(plot_height=600, plot_width=600,
              tools="crosshair,pan,reset,save,wheel_zoom",
              y_range=[50,150], x_range=[1.0,1e4], x_axis_type="log",
              x_axis_label = 'Time since end of outburst (days)', 
              y_axis_label = 'Teff (eV)')
plot.line('x','y',source=source, line_width=2, line_alpha=0.6)

source2 = ColumnDataSource(data=dict(x=tobs, y=Teffobs))
plot.scatter('x','y',source=source2, size=8, marker='circle', fill_color='red', line_color='black')

# Set up the controls 
go_button = Button(label="Go", button_type="success")
go_button.on_click(button_handler)

def_button = Button(label="Defaults", button_type="default")
def_button.on_click(def_handler)

Q_slider = Slider(title="Qimp", value=params['Qimp'], start=0.0, end=30.0, step=0.1)
Q_slider.on_change('value', Q_slider_handler)

Tb_slider = Slider(title="Top temperature (1e8 K)", value=params['Tb']/1e8, start=0.3, end=10.0, step=0.1)
Tb_slider.on_change('value', Tb_slider_handler)

Tc_slider = Slider(title="Core temperature (1e7 K)", value=params['Tc']/1e7, start=1.0, end=30.0, step=0.1)
Tc_slider.on_change('value', Tc_slider_handler)

M_slider = Slider(title="Mass (solar masses)", value=params['M'], start=1.1, end=2.5, step=0.01)
M_slider.on_change('value', M_slider_handler)

R_slider = Slider(title="Radius (km)", value=params['R'], start=8.0, end=16.0, step=0.1)
R_slider.on_change('value', R_slider_handler)

tout_slider = Slider(title="Outburst duration (yrs)", value=params['tout'], start=0.1, end=30.0, step=0.1)
tout_slider.on_change('value', tout_slider_handler)

mdot_slider = Slider(title="Accretion rate (Eddington)", value=params['mdot'], start=0.0, end=1.0, step=0.01)
mdot_slider.on_change('value', mdot_slider_handler)

# Show everything
inputs = widgetbox(Q_slider,Tb_slider, Tc_slider, M_slider, R_slider, mdot_slider, tout_slider, go_button, def_button)
curdoc().add_root(row(inputs, plot, width=1000))
