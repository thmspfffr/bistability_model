'''

Hierarchical network model of perceptual decision making

@author: Klaus Wimmer

wimmer.klaus@googlemail.com

'''


from brian import *
import numpy
import random as pyrandom
from numpy.random import rand as rand
from numpy.random import randn as randn
from scipy.signal import lfilter
from sys import platform
from integration_circuit import make_integration_circuit
import h5py
import os 
from subprocess import call

if platform == 'darwin':
    root_dir = '~/Dropbox/projects/phd/pupmod/decision_network/proc/'
else:
    root_dir = '/home/tpfeffer/pupmod/decision_network/proc/' 

def get_OU_stim(n, tau):
# UO process in discrete time => AR(1) process
      
    a = numpy.exp(-(1.0 / tau))
    i = lfilter(numpy.ones(1),[1.0, -a], numpy.sqrt(1-a*a)*randn(n))   
         
    return i

# ---------------------
# VERSION 21
# ---------------------
#v = 1
#inh_mod = 1
#all_inp_mod = numpy.array([1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8])
#all_bg_mod  = numpy.array([1, 1.05, 1.1, 1.15, 1.2])
bin_size = 20
# ---------------------
# VERSION 2 
# ---------------------
#v = 2
#inh_mod = 1
#all_inp_mod = numpy.array([1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33])
#all_bg_mod = numpy.array([1.13, 1.14, 1.15, 1.16, 1.17])
bin_size = 20
# ---------------------
# VERSION 3 
# ---------------------
#v = 3
#inh_mod = 1
#all_inp_mod = numpy.array([1.37, 1.38, 1.39, 1.4, 1.41, 1.42, 1.43])
#all_bg_mod = numpy.array([1.18, 1.19, 1.20, 1.21, 1.22])
bin_size = 20
# ---------------------
# VERSION 4 
# ---------------------
#v = 4
#inh_mod = 1
#all_inp_mod = numpy.array([1.37, 1.38, 1.39, 1.4, 1.41, 1.42, 1.43])
#all_bg_mod = numpy.array([1.18, 1.19, 1.20, 1.21, 1.22])
#bin_size = 500
# ---------------------
# VERSION 5 (with inp/bg parameters that looked best from v3)
# ---------------------
v = 5
all_inh_mod = numpy.array([0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05])
all_inp_mod = numpy.array([1.39])
all_bg_mod = numpy.array([1.20])
bin_size = 20
# ---------------------
# VERSION 6 (with inp/bg parameters that looked best from v3)
# ---------------------
v = 6
all_inh_mod = numpy.array([0.995, 0.996, 0.997, 0.998, 0.999, 1, 1.001, 1.002, 1.003, 1.004, 1.005])
all_inp_mod = numpy.array([1.39])
all_bg_mod = numpy.array([1.20])
bin_size = 20
ntrls = 1;
# ---------------------
# VERSION 7: similar to 6, but with voltage and larger bin size
# ---------------------
v = 7
all_inh_mod = numpy.array([0.995, 0.996, 0.997, 0.998, 0.999, 1, 1.001, 1.002, 1.003, 1.004, 1.005])
all_inp_mod = numpy.array([1.39])
all_bg_mod = numpy.array([1.20])
bin_size = 100
ntrls = 20;
# ---------------------

for itr in range(ntrls):
  for i_inh_mod in range(all_inh_mod.size):
    for i_inp_mod in range(all_inp_mod.size):
      for i_bg_mod in range(all_bg_mod.size):

          inp_mod = all_inp_mod[i_inp_mod]
          bg_mod = all_bg_mod[i_bg_mod]
          inh_mod = all_inh_mod[i_inh_mod]

          fn = os.path.expanduser(root_dir + 'pupmod_decision_network_inh%d_inp%d_bg%d_tr%d_v%d_processing.txt') % (i_inh_mod, i_inp_mod, i_bg_mod, itr, v)
          if os.path.isfile(fn)==False:
              call(['touch', fn])
          else:
              continue

          if __name__ == '__main__':

              #  initialize  
              defaultclock.reinit()
              clear(True) 


              #------------------------------------------------------------------------------ 
              # Simulation parameters 
              #------------------------------------------------------------------------------ 

              connect_seed = 1284                          # seed for random number generators (set before establishing network connectivity)
              stim_seed = 123                              # seed for random number generators (set before generating common part of stimulus)
              init_seed = 8190                             # seed for random number generators (set before generating private part of stimulus)

              # Timing 
              # stim_on = 500.0 * ms                         # stimulus onset
              # stim_off = 2500.0 * ms                       # stimulus offset   
              # stim_duration = stim_off - stim_on           # duration of stimulus interval
              runtime = 600000.0 * ms                        # total simulation time


              #------------------------------------------------------------------------------
              # Construct hierarchical network
              #------------------------------------------------------------------------------ 

              # Set the seed of the random number generator
              # numpy.random.seed(connect_seed) 
              # pyrandom.seed(connect_seed)

              # Integration circuit
              Dgroups, Dconnections, Dnetfunctions, Dsubgroups = make_integration_circuit(inh_mod,inp_mod,bg_mod)

              # get populations from the integrations circuit
              decisionE = Dgroups['DE']
              decisionI = Dgroups['DI']
              decisionE1 = Dsubgroups['DE1']
              decisionE2 = Dsubgroups['DE2']
              decisionE3 = Dsubgroups['DE3']

              # Sensory network  
              # Sgroups, Sconnections, Ssubgroups = make_sensory_circuit()

              # # get sensory populations
              # sensoryE = Sgroups['SE']
              # sensoryI = Sgroups['SI']
              # sensoryE1 = Ssubgroups['SE1']
              # sensoryE2 = Ssubgroups['SE2']

              # Feed-forward connections from the sensory to the integration circuit

              # wSD = 0.0036                                 # Synaptic weight of feed-forward connections from the corresponding stimulus-encoding population 
                                                          # from the sensory circuit (E1 -> D1, E2 -> D2); 
                                                          # synaptic weight (0.09 nS) is given in multiples of the leak conductance of the excitatory neurons in the integration circuit  

              # C_SE1_DE1 = Connection(sensoryE1, decisionE1, 'gea', weight=wSD, sparseness=0.2, delay=1.0 * ms)
              # C_SE2_DE2 = Connection(sensoryE2, decisionE2, 'gea', weight=wSD, sparseness=0.2, delay=1.0 * ms)

              # Top-down feedback connections from the integration circuit to the sensory circuit

              # b_FB = 0.0                                   # Feedback strength
              # wDS = 0.004 * b_FB                           # Synaptic weight of feedback connections from the integration circuit to the sensory circuit (D1 -> E1, D2 -> E2);
              #                                              # synaptic weight (0.0668 nS) is given in multiples of the leak conductance of the excitatory neurons in the sensory circuit  

              # C_DE1_SE1 = Connection(decisionE1, sensoryE1, 'xe', weight=wDS, sparseness=0.2, delay=1.0 * ms)
              # C_DE2_SE2 = Connection(decisionE2, sensoryE2, 'xe', weight=wDS, sparseness=0.2, delay=1.0 * ms)


              #------------------------------------------------------------------------------
              # Stimulus 
              #------------------------------------------------------------------------------ 

              # Stimulus parameters    
              # I0 = 0.08 * nA                               # Mean input current for zero-coherence stimulus 
              # c = 0.0                                      # Stimulus coherence (between 0 and 1) 
              # mu_E1 = +0.25                                # Average additional input current to E1 at highest coherence (c = 1)
              # mu_E2 = -0.25                                # Average additional input current to E2 at highest coherence (c = 1)    
              # sigma = 1.0																	 # Amplitude of temporal modulations of the stimulus
              # sigma_stim = 0.212 * sigma                   # S.d. of modulations of stimulus inputs
              # sigma_ind = 0.212 * sigma                    # S.d. of modulations in individual inputs
              # tau_stim = 20.0 * ms                         # Correlation time constant of Ornstein-Uhlenbeck process

              # Generate stimulus
              # set seed of random number generator (in order to generate a specific stimulus each time)
              # numpy.random.seed(stim_seed)    
              # pyrandom.seed(stim_seed)

              # "common part" of the stimulus
              # z1 = get_OU_stim(stim_duration/ms, tau_stim/ms)
              # z1 = numpy.tile(z1,(len(sensoryE1),1))
              # z2 = get_OU_stim(stim_duration/ms, tau_stim/ms)
              # z2 = numpy.tile(z2,(len(sensoryE2),1))

              # set seed of random number generator (in order to generate a specific stimulus each time)
              # numpy.random.seed(init_seed)    
              # pyrandom.seed(init_seed)

              # "private part" - part of the stimulus for each neuron, different in each trial               
              # zk1 = get_OU_stim(int(stim_duration/ms * len(sensoryE1)), tau_stim/ms)  
              # zk1 = numpy.asarray(zk1).reshape(len(sensoryE1), stim_duration/ms)
              # zk2 = get_OU_stim(int(stim_duration/ms * len(sensoryE2)), tau_stim/ms)  
              # zk2 = numpy.asarray(zk2).reshape(len(sensoryE2), stim_duration/ms)

              # stimulus (time series with dt = 1ms)
              # most general case: different input to each neuron in each time step 
              # i1 = I0 * (1 + c * mu_E1 + sigma_stim * z1 + sigma_ind * zk1)     
              # i2 = I0 * (1 + c * mu_E2 + sigma_stim * z2 + sigma_ind * zk2)     
              # ii = numpy.zeros((len(sensoryI),stim_duration/ms))     # no current input to inh population

              # Stimulus-related external current input     
              # myclock=Clock(dt=1*ms)
              # @network_operation(myclock)
              # def update_input():  
              #     if myclock.t >= stim_on and myclock.t < stim_off:           
              #         sensoryE1.I = i1[:,int( (myclock.t - stim_on) / (1 * ms))] * amp
              #         sensoryE2.I = i2[:,int( (myclock.t - stim_on) / (1 * ms))] * amp
              #         sensoryI.I = ii[:,int( (myclock.t - stim_on) / (1 * ms))] * amp
              #     else:
              #         sensoryE1.I = 0 * nA
              #         sensoryE2.I = 0 * nA
              #         sensoryI.I = 0 * nA

              #------------------------------------------------------------------------------
              # Initial conditions and Monitors
              #------------------------------------------------------------------------------

              # --- set seed of random number generator to a different value in each run
              # np_seed = int(1587.47)
              # numpy.random.seed(np_seed)
              # py_seed = int(4736.28) 
              # pyrandom.seed(py_seed)

              # ---- set initial conditions (random)
              decisionE.gen = decisionE.gen * (1 + 0.2 * rand(decisionE.__len__()))
              decisionI.gen = decisionI.gen * (1 + 0.2 * rand(decisionI.__len__()))
              decisionE.V = decisionE.V + rand(decisionE.__len__()) * 2 * mV
              decisionI.V = decisionI.V + rand(decisionI.__len__()) * 2 * mV

              # ---- set initial conditions (random)
              # sensoryE.V = -50.0 * mV - 2 * mV + rand(sensoryE.__len__()) * 2 * mV
              # sensoryI.V = -50.0 * mV - 2 * mV + rand(sensoryI.__len__()) * 2 * mV
              # sensoryE.gea = 0.05 * (1 + rand(sensoryE.__len__()) * 0.2)
              # sensoryI.gea = 0.05 * (1 + rand(sensoryI.__len__()) * 0.2)

              # record spikes of excitatory neurons
              # S_DE1 = SpikeMonitor(decisionE1, record=True)
              # S_DE2 = SpikeMonitor(decisionE2, record=True)
              # S_SE1 = SpikeMonitor(sensoryE1, record=True)
              # S_SE2 = SpikeMonitor(sensoryE2, record=True)

              # record instantaneous populations activity
              R_DE1 = PopulationRateMonitor(decisionE1, bin=bin_size*ms)
              R_DE2 = PopulationRateMonitor(decisionE2, bin=bin_size*ms)
              # R_SE1 = PopulationRateMonitor(sensoryE1, bin=bin_size*ms)
              # R_SE2 = PopulationRateMonitor(sensoryE2, bin=bin_size*ms)

              V_DE1 = StateMonitor(decisionE1,'V',record=True,timestep=bin_size)
              V_DE2 = StateMonitor(decisionE2,'V',record=True,timestep=bin_size)

              #------------------------------------------------------------------------------
              # Run the simulation
              #------------------------------------------------------------------------------

              # construct network       
              net = Network(Dgroups.values(), Dconnections.values(), Dnetfunctions,
                          R_DE1, R_DE2, V_DE1, V_DE2)
              net.prepare()
              net.run(runtime)        

              rate_D1 = []
              rate_D2 = []

              for ineuron in range(0,len(R_DE2.rate)):
                  rate_D1 = numpy.append(rate_D1,R_DE1.rate[ineuron], axis=None)
                  rate_D2 = numpy.append(rate_D2,R_DE2.rate[ineuron], axis=None)

              #print rate_D1.mean()
              #print rate_D2.mean()

              # print("Saving output...")
              hf = h5py.File(os.path.expanduser(root_dir + 'pupmod_decision_network_inh%d_inp%d_bg%d_tr%d_v%d.h5' % (i_inh_mod,i_inp_mod,i_bg_mod,itr,v) ), 'w')
              hf.create_dataset('rate_D1', data=rate_D1)
              hf.create_dataset('rate_D2', data=rate_D2)
              hf.close()

              volt_D1 = V_DE1.values
              volt_D2 = V_DE2.values

              hf = h5py.File(os.path.expanduser(root_dir + 'pupmod_decision_network_inh%d_inp%d_bg%d_tr%d_voltage_v%d.h5' % (i_inh_mod,i_inp_mod,i_bg_mod,itr,v) ), 'w')
              hf.create_dataset('volt_D1', data=volt_D1)
              hf.create_dataset('volt_D2', data=volt_D2)
              hf.close()
