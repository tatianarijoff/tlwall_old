import scipy.constants as const
import numpy as np
from scipy.special import i0, i1 ,k0,  k1, iv
import cmath
import sys
if sys.version_info.major == 2:
    import ConfigParser as configparser
else:
    import configparser

accuracy_factor = 0.3
np.seterr(over='ignore')
np.seterr(invalid='ignore')

def welcome_message():
    print ('')
    print ('*************************************************************')
    print ('*     PyTLWall                                              *')
    print ('*                                                           *')
    print ('*  C. Zannini and T. Rijoff                                 *')
    print ('*************************************************************')
    print (' usage: ')
    print ('       exec_pytlwall    --   for interactive mode            ')
    print ('       exec_pytlwall cfg_file  -- to  read information from config file  ')
    print (' ')

def safe_float(x):
    try:
        y = float(x)
        return y
    except ValueError:
        return 0.
    
def safe_int(x):
    try:
        y = int(x)
        return y
    except ValueError:
        return 0
        
def createdir(savename):
    import os
    dir_arr = savename.split('/')
    mydir = './'
    for i in range(len(dir_arr)-1):
        mydir += dir_arr[i] + '/'
    os.makedirs(mydir)
                        
class Freq(object):
    def __init__(self, fmin = 0, fmax = 0, dfexp = 0 , filename = None):
        self.default_fmin = 0 
        self.default_fmax = 8 
        self.default_dfexp = 2
        self.freq_file = None
        
        self.fmin = safe_int(fmin) 
        self.fmax = safe_int(fmax) 
        self.dfexp = safe_int(dfexp) 
        if filename != None:
            self.freq_file = filename
            self.freq = self.read_freq_from_file(self.freq_file)
        elif self.fmin != 0 and self.fmax != 0 and self.dfexp != 0:
            self.freq = self.calc_freq_array(self.fmin, self.fmax, self.dfexp)
        else:
            self.freq = np.array([])
    
    def read_freq_from_file( self, freq_file):
        f = np.array([])
        for line in open(freq_file,'r'):
            newf = map (safe_float, line.split())
            f = np.append(f, newf)
        f = (np.unique(f))
        return f

    def calc_freq_array( self, fmin, fmax, dfexp):
        f = np.array([])
        for p in np.arange(1,fmax-fmin +1):
            v1 =(1+(10**(1-dfexp)))* 10.**(fmin-1+p)
            v2 = 10.**(fmin+p)
            v3 = 10.**(fmin-1+p-(dfexp-1))
            f = np.append(f, np.arange( v1 , v2 +v3, v3 ))

        return f
        

    #
    #  Input
    #
    def ask_freq_info(self):
        choice = raw_input('Which is the minimum frequency exponent? (default %d)  ' %self.default_fmin)
        if choice == '':
            self.fmin = self.default_fmin
        else:
            try:
                self.fmin = int(choice)
            except ValueError:
                print ('%s is not allowed %d is used' %(choice, self.default_fmin))
                self.fmin = self.default_fmin
        choice = raw_input('Which is the maximum frequency exponent? (default %d)  ' %self.default_fmax)
        if choice == '':
            self.fmax = self.default_fmax
        else:
            try:
                self.fmax = int(choice)
            except ValueError:
                print ('%s is not allowed %d is used' %(choice, self.default_fmax))
                self.fmax = self.default_fmax
        choice = raw_input('Which is the frequency exponent step? (default %d)  ' % self.default_dfexp)
        if choice == '':
            self.dfexp = self.default_dfexp
        else:
            try:
                self.dfexp = int(choice)
            except ValueError:
                print ('%s is not allowed %d is used' %(choice, self.default_dfexp))
                self.dfexp = self.default_dfexp
        
        self.freq = self.calc_freq_array( self.fmin, self.fmax, self.dfexp)
    
    
    def read_cfg(self, cfg_file):
        config = configparser.ConfigParser()
        try:
            config.readfp(open(cfg_file))
            self.read_freq_cfg(config)
        except IOError:
            print ('')
            print ('The file %s does not exist, try again!' %cfg_file)


    def read_freq_cfg(self, config):
        try:
            choice = config.get('frequency_info', 'freq_file')
            self.freq_file = choice
        except configparser.NoSectionError:
            print ('')
            print ("The given configurator does not contain the section frequency_info, try with another configurator file")
            return
        except configparser.MissingSectionHeaderError:
            print ('')
            print ("The given configurator does not contain the section frequency_info, try with another configurator file")
            return
        except configparser.NoOptionError:
            try:
                choice = config.get('frequency_info', 'fmin')
                self.fmin = safe_int(choice)
            except configparser.NoOptionError:
                print ('')
                print ("The given configurator does not contain the needed frequency_info, try with another configurator file")
                return
        try:
            choice = config.get('frequency_info', 'fmax')
            self.fmax = safe_int(choice)
        except configparser.NoOptionError:
            print ("Max frequency needed default value is used %d" %self.default_fmax)
            self.fmax = self.default_fmax 
        
        try:
            choice = config.get('frequency_info', 'dfexp')
            self.dfexp = safe_int(choice)
        except configparser.NoOptionError:
            print ("Frequency step needed default value is used %d" %self.default_dfexp)
            self.dfexp = self.default_dfexp 
        self.freq = self.calc_freq_array( self.fmin, self.fmax, self.dfexp)
    #
    #  Output 
    #
    def generate_freq_cfg(self, filename, flag_write_append):
        fd = open(filename, flag_write_append)
        fd.write('[frequency_info] \n')
        fd.write('; you can put the frequency exponent information (minimum frequency "fmin", maximum frequency "fmax", number of points exponent per decade "dfexp") \n')
        fd.write('; or a file with the wanted frequency (freq_file) \n')
        if self.freq_file != None:
            fd.write('freq_file = %s  \n' %self.freq_file)
        else:
            fd.write('fmin = %d  \n' %self.fmin)
            fd.write('fmax = %d  \n' %self.fmax)
            fd.write('dfexp = %d  \n' %self.dfexp)
            
class Acc_El():
    def __init__(self, pipe_len = 0. , pipe_radius = 0., betax = 0., betay = 0., form_factor = 'ROUND', component_name = None ):
        self.default_component_name = 'el'
        self.default_pipe_len = 1.
        self.default_pipe_radius = 0.01
        self.default_form_factor='ROUND'
        self.default_betax=1.
        self.default_betay=1.
        self.default_nbr_layers = 0
        self.default_thick = 0.1
        self.default_type = 'CW'
        self.default_boundary_type = 'V'  # considered if nbr_layers >= 1 otherwise default_type is used
        self.default_sigmaDC=1.e6        
        self.default_epsr = 1.
        self.default_muinf = 0.
        self.default_k = float('inf')
        self.default_tau = 0.
        
        self.pipe_len = pipe_len
        self.pipe_radius = pipe_radius
        self.component_name = component_name
    
        self.default_cfg='res_wall_component.cfg'
        #
        #   define default value for optic betas
        #
        self.betax = betax
        self.betay = betay
        self.form_factor = form_factor

        #
        #  define the array for the form factor
        #
        self._dipx_form_factor={}
        self._dipy_form_factor={}
        self._quadx_form_factor={}
        self._quady_form_factor={}
        
        self._dipx_form_factor['ROUND'] = 1.
        self._dipy_form_factor['ROUND'] = 1.
        self._quadx_form_factor['ROUND'] = 0.
        self._quady_form_factor['ROUND'] = 0.

        self._dipx_form_factor['RECTANGULAR'] = np.pi * np.pi / 24.
        self._dipy_form_factor['RECTANGULAR'] = np.pi * np.pi / 12.
        self._quadx_form_factor['RECTANGULAR'] = - np.pi * np.pi / 24.
        self._quady_form_factor['RECTANGULAR'] = np.pi * np.pi / 24.

        diam_fact = 1.39
        #diam_fact = 1.
        self._dipx_form_factor['DIAMOND'] = diam_fact * self._dipx_form_factor['RECTANGULAR']
        self._dipy_form_factor['DIAMOND'] = diam_fact * self._dipy_form_factor['RECTANGULAR']
        self._quadx_form_factor['DIAMOND'] = diam_fact * self._quadx_form_factor['RECTANGULAR']
        self._quady_form_factor['DIAMOND'] = diam_fact * self._quady_form_factor['RECTANGULAR']

        self._dipx_form_factor['OBLONG'] = self._dipx_form_factor['RECTANGULAR']
        self._dipy_form_factor['OBLONG'] = self._dipy_form_factor['RECTANGULAR']
        self._quadx_form_factor['OBLONG'] = self._quadx_form_factor['RECTANGULAR']
        self._quady_form_factor['OBLONG'] = self._quady_form_factor['RECTANGULAR']

    @property
    def form_factor(self):
        return self._form_factor
        
    @form_factor.setter
    def form_factor(self, newform):
        
        if newform.upper() ==  'RECTANGULAR' or newform.upper() ==  'DIAMOND' or newform.upper() ==  'OBLONG':
            self._form_factor = newform.upper()
        else:
            self._form_factor = 'ROUND'

    def calc_sigmaAC(self, f, sigmaDC, tau):
        """ sigmaAC = sigmaDC / (1 + j 2 pi tau f ) """
        return  sigmaDC / ( 1. + 2.j *const.pi  * tau * f)
        
    def calc_eps(self, f, sigma_ac, epsr):
        """ eps = const.epsilon_0 * epsr - j (sigma_ac / (2 pi f))   """
        return (const.epsilon_0 *  epsr) - 1.j * ( sigma_ac / (2 * const.pi * f))

    def calc_sigmaPM(self, f, eps, sigma_ac):
        """ sigmaPM = sqrt((2 pi  f  eps)^2 + sigma_ac^2 ) """
        return np.sqrt((2 * const.pi * f * eps)**2 + sigma_ac**2 ) 
        
    def calc_delta(self, f, mu, eps, sigma_ac):
        """ delta = sqrt(2 / (2  pi f  mu sigma_ac + j  mu eps ( 2 pi f)^2 )) """
        return  np.sqrt(2 / (2 * const.pi * f * mu * sigma_ac + 1.j * mu * eps *( 2 * const.pi * f)**2 )) 

    def calc_deltaM(self, f, mu, eps, sigma_ac):
        """ deltaM = sqrt(2 / (2 pi f mu sigma_ac - j mu eps ( 2 pi f)^2)) """
        return  np.sqrt(2 / (2 * const.pi * f * mu * sigma_ac -1.j * mu * eps *( 2 * const.pi * f)**2)) 

    def calc_mur(self, f, k, muinf):
        """ mur = 1 +  (muinf / ( 1 + j ( f / k)))  
                   if k = inf  """
        return  1 +  (muinf / ( 1 + 1.j * (f / k )))

    #
    #  Input 
    #
    def ask_acc_el_info(self):
        choice = raw_input('Which is the component name (default "%s")?  ' %self.default_component_name)
        self.save_component_name(choice)
        choice = raw_input('Which is the %s length (default %f [m])?  ' %(self.component_name, self.default_pipe_len))
        self.save_pipe_length(choice)
        choice = raw_input('Which is the %s radius (default %f [m])?  ' %(self.component_name, self.default_pipe_radius))
        self.save_pipe_radius(choice)
        choice = raw_input('Which is the %s shape (default "%s")?  ' %(self.component_name, self.default_form_factor))
        self.save_form_factor(choice)
        choice = raw_input('Which is the %s horiz betatron function (default %f [m])?  '%(self.component_name, self.default_betax))
        self.save_betax(choice)
        choice = raw_input('Which is the %s horiz betatron function (default %f [m])? '%(self.component_name, self.default_betay))
        self.save_betay(choice)
        self.ask_layers()
        choice = raw_input('Choose a name to save the data (default %s , digit NO if you don\'t want a file) '%self.default_cfg )
        if choice.upper() != 'NO':
            if choice != '':
                filename = choice
            else:
                filename = self.default_cfg
            flag_write_append='w'
            self.generate_acc_component_cfg(filename, flag_write_append)

        print ("=========================================================================")
        print ("accelerator component definition: done ")
        print ("=========================================================================")


    def save_component_name(self, choice):
        if choice == '':
            self.component_name = self.default_component_name        
        else:
            self.component_name = choice

    def save_pipe_length(self, choice):
        if choice == '':
            self.pipe_len = self.default_pipe_len        
        else:
            try:
                self.pipe_len = float(choice)
            except ValueError:
                print ('%s is not a good length default is used' %choice)
                self.pipe_len = self.default_pipe_len

    def save_pipe_radius(self, choice):
        if choice == '':
            self.pipe_radius = self.default_pipe_radius        
        else:
            try:
                self.pipe_radius = float(choice)
            except ValueError:
                print ('%s is not a good radius default is used' %choice)
                self.pipe_radius = self.default_pipe_radius
                            
    def save_form_factor(self, choice):
        if choice.upper() == 'RECTANGULAR' or choice.upper() == 'DIAMOND' or choice.upper() == 'OBLONG':
            self.form_factor = choice.upper()        
        else:
            self.form_factor = self.default_form_factor
        
    def save_betax(self, choice):
        if choice == '':
            self.betax = self.default_betax        
        else:
            try:
                self.betax = float(choice)
            except ValueError:
                print ('%s is not a good horiz betatron function default is used' %choice)
                self.betax = self.default_betax

    def save_betay(self, choice):
        if choice == '':
            self.betay = self.default_betay        
        else:
            try:
                self.betay = float(choice)
            except ValueError:
                print ('%s is not a good horiz betatron function default is used' %choice)
                self.betay = self.default_betay

    def ask_layers(self):    
        choice = raw_input('Which is the number of layers for %s (default %d)?  '%(self.component_name, self.default_nbr_layers))
        try:
            self.nbr_layers = int(choice)
        except ValueError:
                print ('%s is not a good number of layers default is used' %choice)    
                self.nbr_layers = self.default_nbr_layers        
        self.layers = []
        for i in xrange(self.nbr_layers + 1):
            self.layers.append({})
            if i < self.nbr_layers :
                string_ask = 'layer %d'%i+1
                default_type = self.default_type
                
                # thickness is needed only for normal layers (not for the boundary)
                choice = raw_input('Which is the thickness of the %s (default %f [m])?  '%(string_ask, self.default_thick))
                self.save_thikness(i, choice)
                
            else:
                string_ask = 'boundary'
                if i == 0 :  # boundary is the only layer default is cw 
                    default_type = self.default_type
                else:  
                    default_type = self.default_boundary_type
                
            choice = raw_input('Which is the material of the %s(default %s, possible values "CW" (conductive wall)  "V" (Vacuum) "PEC" (Perfect electric conductor)  )?  '%(string_ask, default_type))
            self.save_type(i, choice, default_type)
            if self.layers[i]['type'].upper() != 'V' and self.layers[i]['type'].upper() != 'PEC':
                choice = raw_input('Which is the DC conductivity  of the %s (default %e [S/m])?  '%(string_ask, self.default_sigmaDC))
                self.save_sigmaDC(i , choice)
                self.save_epsr(i ,  choice)
                choice = raw_input('Which is the real relative permittivity %s (default %e)?  '%(string_ask, self.default_epsr))
                self.save_muinf(i ,  choice)
                choice = raw_input('Which is the real relative permeability %s (default %e)?  '%(string_ask, self.default_muinf))
                self.save_k(i ,  choice)
                choice = raw_input('Which is the relaxation frequency for permeability %s (default %e [Hz])?  '%(string_ask, self.default_k))
                self.save_tau(i ,  choice)    
                choice = raw_input('Which is the relaxation time for permittivity %s (default %e [s])?  '%(string_ask, self.default_tau))

    def save_thikness(self, i , choice):
            try:
                self.layers[i]['thick'] = float(choice)
            except ValueError:
                print ('%s is not a good thickness default is used' %choice)
                self.layers[i]['thick'] = self.default_thick
                        
    def save_type(self, i, choice, default_type):

        if choice.upper() == 'CW' or choice.upper() == 'V' or choice.upper() == 'PEC' :
            self.layers[i]['type'] = choice.upper()
        else:
            print ('%s is not a good material default is used' %choice)
            self.layers[i]['type'] = default_type

    def save_sigmaDC(self, i , choice):
        if choice.lower() == 'inf' or choice.lower() == 'infinity':
            print ("Infinity is not allowed for sigma DC default value is used")
            self.layers[i]['sigmaDC'] = self.default_sigmaDC
        else:
            try:
                self.layers[i]['sigmaDC'] = float(choice)
            except ValueError:
                print ('%s is not a good sigma DC value default is used' %choice)
                self.layers[i]['sigmaDC'] = self.default_sigmaDC

    def save_epsr(self, i , choice):
        if choice.lower() == 'inf' or choice.lower() == 'infinity':
            self.layers[i]['epsr'] = float("inf")
        else:
            try:
                self.layers[i]['epsr'] = float(choice)
            except ValueError:
                print ('%s is not a good epsr value default is used' %choice)
                self.layers[i]['epsr'] = self.default_epsr
        
    def save_muinf(self, i , choice):
        if choice.lower() == 'inf' or choice.lower() == 'infinity':
            self.layers[i]['muinf'] = float("inf")
        else:
            try:
                self.layers[i]['muinf'] = float(choice)
            except ValueError:
                print ('%s is not a good muinf value default is used' %choice)
                self.layers[i]['muinf'] = self.default_muinf

    def save_k(self, i , choice):
        if choice.lower() == 'inf' or choice.lower() == 'infinity':
            self.layers[i]['k'] = float("inf")
        else:
            try:
                self.layers[i]['k'] = float(choice)
            except ValueError:
                print ('%s is not a good k value default is used' %choice)
                self.layers[i]['k'] = self.default_k
        if self.layers[i]['k'] == 0:
            print ("0 is not allowed for the relaxation frequency default is used")
            self.layers[i]['k'] = self.default_k

    def save_tau(self, i , choice):
        try: 
            self.layers[i]['tau'] = float(choice)
        except ValueError:
            try:
                if choice.lower() == 'inf' or choice.lower() == 'infinity':
                    self.layers[i]['tau'] = float("inf")
            except ValueError:
                print ('%s is not a good tau value default is used' %choice)
                self.layers[i]['tau'] = self.default_tau                                                
    
    def read_cfg(self, cfg_file):

        config = configparser.ConfigParser()
        try:
            config.readfp(open(cfg_file))
            self.read_acc_component_cfg(config)
        except IOError:
            print ('')
            print ('The file %s does not exist, try again!' %cfg_file)


    def read_acc_component_cfg(self, config):
        try:
            choice = config.get('base_info', 'component_name')
            self.save_component_name(choice)
        except configparser.NoSectionError:
            print ('')
            print ("The given configurator does not contain the section base_info, try with another configurator file")
            return
        except configparser.MissingSectionHeaderError:
            print ('')
            print ("The given configurator does not contain the section base_info, try with another configurator file")
            return
        except configparser.NoOptionError:
            self.save_component_name(self.default_component_name)
        
        try:
            choice = config.get('base_info', 'pipe_len')
            self.save_pipe_length(choice)
        except configparser.NoOptionError:
            self.save_pipe_length(self.default_pipe_len)
        try:
            choice = config.get('base_info', 'pipe_radius')
            self.save_pipe_radius(choice)
        except configparser.NoOptionError:
            self.save_pipe_radius(self.default_pipe_radius)
        try:
            choice = config.get('base_info', 'betax')
            self.save_betax(choice)
        except configparser.NoOptionError:
            self.save_betax(self.default_betax)        
        try:
            choice = config.get('base_info', 'betay')
            self.save_betay(choice)
        except configparser.NoOptionError:
            self.save_betay(self.default_betay)                
        try:
            choice = config.get('base_info', 'form_factor')
            self.save_form_factor(choice)
        except configparser.NoOptionError:
            self.save_form_factor(self.default_form_factor)

        i = 0
        self.layers = []
        while config.has_section('layer%d'%i):
            self.layers.append({})
            try:
                choice = config.get('layer%d'%i, 'thick')
                self.save_thikness(i, choice)
            except configparser.NoOptionError:
                self.save_thikness(i, self.default_thick)
            try:
                choice = config.get('layer%d'%i, 'type')
                self.save_type(i, choice, self.default_type)
            except configparser.NoOptionError:
                self.save_type(i, self.default_type)
            if self.layers[i]['type'].upper() != 'V' and self.layers[i]['type'].upper() != 'PEC':
                try:
                    choice = config.get('layer%d'%i, 'sigmaDC')
                    self.save_sigmaDC(i, choice)
                except configparser.NoOptionError:
                    self.save_sigmaDC(i, self.default_sigmaDC)
                try:
                    choice = config.get('layer%d'%i, 'epsr')
                    self.save_epsr(i, choice)
                except configparser.NoOptionError:
                    self.save_epsr(i, self.default_epsr)                
                try:
                    choice = config.get('layer%d'%i, 'muinf')
                    self.save_muinf(i, choice)
                except configparser.NoOptionError:
                    self.save_muinf(i, self.default_muinf)            
                try:
                    choice = config.get('layer%d'%i, 'k')
                    self.save_k(i, choice)
                except configparser.NoOptionError:
                    self.save_k(i, self.default_k)        
                try:
                    choice = config.get('layer%d'%i, 'tau')
                    self.save_tau(i, choice)
                except configparser.NoOptionError:
                    self.save_tau(i, self.default_tau)                                                
                i = i+1
            
        #
        # Add boundary
        #    
        self.layers.append({})
        if i == 0:
            default_type = self.default_type
        else:
            default_type = self.default_boundary_type
        try:
            choice = config.get('boundary', 'type')
            self.save_type(i, choice, default_type)
        except configparser.NoOptionError:
            self.save_type(i, default_type)
        if self.layers[i]['type'].upper() != 'V' and self.layers[i]['type'].upper() != 'PEC':
            try:
                choice = config.get('boundary', 'sigmaDC')
                self.save_sigmaDC(i, choice)
            except configparser.NoOptionError:
                self.save_sigmaDC(i, self.default_sigmaDC)
            try:
                choice = config.get('boundary', 'epsr')
                self.save_epsr(i, choice)
            except configparser.NoOptionError:
                self.save_epsr(i, self.default_epsr)                
            try:
                choice = config.get('boundary', 'muinf')
                self.save_muinf(i, choice)
            except configparser.NoOptionError:
                self.save_muinf(i, self.default_muinf)            
            try:
                choice = config.get('boundary', 'k')
                self.save_k(i, choice)
            except configparser.NoOptionError:
                self.save_k(i, self.default_k)        
            try:
                choice = config.get('boundary', 'tau')
                self.save_tau(i, choice)
            except configparser.NoOptionError:
                self.save_tau(i, self.default_tau)
                
        self.nbr_layers=len(self.layers)-1
        
        
                                    
    #
    #  Output 
    #
    def generate_acc_component_cfg(self, filename, flag_write_append):
        fd = open(filename, flag_write_append)
        fd.write('[base_info] \n')
        fd.write('; Component name (a label to identify the element, not necessary) \n')
        fd.write('component_name = %s  \n' %self.component_name)
        fd.write('; Pipe radius in meter \n')
        fd.write('pipe_radius = %e  \n' %self.pipe_radius)
        fd.write('; Length in meter          \n')
        fd.write('pipe_len = %e  \n' %self.pipe_len)
        fd.write('; average twiss horizontal beta in the component (if not indicated default is %e) \n' %self.default_betax)
        fd.write('betax = %e \n' %self.betax)
        fd.write('; average twiss vertical beta in the component (if not indicated default is %e) \n' %self.default_betay)
        fd.write('betay = %e \n' %self.betay)
        fd.write('; shape of the vacuum chamber, allowed values: ROUND(default), RECTANGULAR, DIAMOND, OBLONG \n')
        fd.write('form_factor = %s \n' %self.form_factor)
        fd.write('[layers_info] \n')
        fd.write('; number of layers before the boundary \n')
        fd.write('nbr_layers = %d \n' %self.nbr_layers)
        for i in xrange(self.nbr_layers):
            fd.write('[layer%d] \n' %i)
            fd.write('; layer type \n')
            fd.write('; if type is "V"(vacuum), "PEC"(Perfect electric conductor) or \n')
            fd.write('; "PMC" (Perfect magnetic conductor) only the thikness must be specified \n')
            fd.write('type : %s             \n' %self.layers[i]['type'] )
            fd.write('; layer thickness \n')
            fd.write('thick : %e           \n' %self.layers[i]['thick'] )
            fd.write('; mu infinity \n')
            fd.write('muinf : %e               \n' %self.layers[i]['muinf'] )
            fd.write('; Relaxation frequency for permeability in Hz ZERO NOT ALLOWED \n')
            fd.write('k : %f  \n' %self.layers[i]['k'] )
            fd.write('; Real relative permittivity   \n')
            fd.write('epsr : %e        \n' %self.layers[i]['epsr'])
            fd.write('; DC conductivity in S/m INFINITY NOT ALLOWED      \n')
            fd.write('sigmaDC :%e       \n'%self.layers[i]['sigmaDC'])
            fd.write('; Relaxation time for permittivity in second  \n')
            fd.write('tau : %e              \n' %self.layers[i]['tau'])
        i = self.nbr_layers
        fd.write('[boundary] \n')
        fd.write('type: %s  \n' %self.layers[i]['type'])
        if self.layers[i]['type'] != 'V' and self.layers[i]['type'] != 'PEC':    
            fd.write('; mu infinity \n')
            fd.write('muinf: %e  \n' %self.layers[i]['muinf'])
            fd.write('; Relaxation frequency for permeability in Hz ZERO NOT ALLOWED \n')
            fd.write('k: %e      \n' %self.layers[i]['k'])
            fd.write('; DC conductivity in S/m INFINITY NOT ALLOWED      \n')
            fd.write('sigmaDC :%e  \n' %self.layers[i]['sigmaDC'])
            fd.write('; Real relative permittivity   \n')
            fd.write('epsr : %e \n' %self.layers[i]['epsr'])
            fd.write('; Relaxation time for permittivity in second  \n')
            fd.write('tau : %e   \n' %self.layers[i]['tau'])

        fd.close()

class Beam(object):
    def __init__(self, Ekin = 0., betarel = 1., gammarel = float('inf'), test_beam_shift = 0., charge = 0., mass = 0., particle_name = None):
        #~ self.default_betarel=1. - 1.e-16
        self.default_betarel = 1.
        self.default_gammarel=1.
        self.default_Ekin = 0.
        self.default_test_beam_shift = 0.001
        self.default_particle_name = 'proton'
        
        self.particle_name = None
        self._m = const.m_p
        if particle_name != None:
            self.save_particle_name(particle_name)
            
        self.Ekin = Ekin 
        if betarel != 1:
            self.betarel = betarel
        if gammarel != float('inf'):
            self.gammarel = gammarel
            
        self.test_beam_shift = test_beam_shift
    #
    #   Kinematics characterization (gammarel, betarel)
    #
    
    @property
    def betarel(self):
        return self._betarel        
    
    @betarel.setter
    def betarel(self, newbeta):
        try:
            self._betarel= float(newbeta)
        except ValueError:
            print ('%s is not a good value for relativistic beta default is used %.2f' %(newbeta, self.default_betarel))
            self._betarel = self.default_betarel
        
        if self._betarel <= 0:
            print ('%s is not a good value for relativistic beta default is used %.2f ' %(self._betarel, self.default_betarel))
            self._betarel = self.default_betarel
            
        try:
            self._gammarel =  np.sqrt(1./(1.- (self._betarel*self._betarel)))
        except ZeroDivisionError:
            self._betarel = self.default_betarel    
            self._gammarel =  float('inf')
        self._Ekin = self._m * self._betarel**2 / 2
        
        
    @property
    def gammarel(self):
        return self._gammarel
        
    @gammarel.setter
    def gammarel(self, newgamma):
        try:
            self._gammarel = float(newgamma)
        except ValueError:
            print ('%s is not a good value for relativistic gamma default is used ' %newgamma)
            self._gammarel = self.default_gammarel

        self._Ekin = self.gammarel * self._m
        self._betarel = np.sqrt(1 - (1/(self._gammarel*self._gammarel)))
        
    @property
    def Ekin(self):
        return self._Ekin
        
    @Ekin.setter
    def Ekin(self, newEkin):
        try:
            self._Ekin = float(newEkin)
        except ValueError:
            print ('%s is not a good value for kinetic energy default is used ' %newEkin)
            self._Ekin = self.default_Ekin

        self._gammarel = (self._Ekin/  self._m) + 1.
        self._betarel = np.sqrt(1 - (1/(self._gammarel*self._gammarel)))

    #
    # Input
    #
    def ask_beam_info(self):
        choice = raw_input('Which are the beam particle (default %s)' %self.default_particle_name)
        self.save_particle_name(choice)
        
        while choice.lower() != 'ekin' and choice.lower() != 'beta' and choice.lower() != 'gamma':
            choice = raw_input('What do you want to insert?  kinetic energy (Ekin), beta or gamma? ' )
        if choice.lower() == 'ekin':
            choice = raw_input('Which is the beam kinetic energy? (default %e [MeV])  ' %self.default_kinetic_energy)
            self.Ekin = choice
        elif choice.lower() == 'beta':
            choice = raw_input('Which is the beam relativistic beta? (default %e) ' %self.default_betarel )
            self.betarel = choice
        else:
            choice = raw_input('Which is the beam relativistic gamma? (default %e) ' %self.default_gammarel )
            self.gammarel= choice
        choice = raw_input('Which is the source test separation (default %e [m])?  ' %self.default_test_beam_shift)
        try:
            self.test_beam_shift = float(choice)
        except ValueError:
            if choice != '':
                print ('%s is not a good value for test beam separation default is used' %choice)
            self.test_beam_shift = self.default_test_beam_shift
        
    def save_particle_name(self, choice):
        self.particle_name = choice
        if choice == 'proton' or choice == 'antiproton' or choice == 'pbar' or choice == 'p' :            
            self._m = const.physical_constants['proton mass energy equivalent in MeV'][0]
        elif choice == 'electron':        
            self._m = const.physical_constants['electron mass energy equivalent in MeV'][0]
        elif choice=='':
            self.particle_name = self.default_particle_name
            self._m = const.physical_constants['proton mass energy equivalent in MeV'][0]
        else:
            print ('%s not managed yet, supposed %s' %(choice, self.default_particle_name))
            self._m = const.physical_constants['proton mass energy equivalent in MeV'][0]
        
        
            
    def read_cfg(self, cfg_file):

        config = configparser.ConfigParser()
        try:
            config.readfp(open(cfg_file))
            self.read_beam_cfg(config)
        except IOError:
            print ('')
            print ('The file %s does not exist, try again!' %cfg_file)


    def read_beam_cfg(self, config):
        
        try:
            choice = config.get('beam_info', 'particle_name')
            self.save_particle_name(choice)
        except configparser.NoSectionError:
            print ('')
            print ("The given configurator does not contain the section beam_info")
            self.save_particle_name(self.default_particle_name)
            self.betarel = self.default_betarel        
        except configparser.NoOptionError:
            self.save_particle_name(self.default_particle_name)

        try:
            choice = config.get('beam_info', 'betarel')
            self.betarel = choice
            
        except configparser.NoOptionError:
            try:
                choice = config.get('beam_info', 'gammarel')
                self.gammarel = choice
            except configparser.NoOptionError:
                try:
                    choice = config.get('beam_info', 'Ekin')
                    self.Ekin = choice    
                except configparser.NoOptionError:
                    print ('')
                    print ("In the given file there are not informations about the beam dynamic, default values are used")
                    self.betarel = self.default_betarel
        try:
            choice = config.get('beam_info', 'test_beam_shift')
        except configparser.NoOptionError:
            print ('')
            print ("Source test separation not inserted, default is used")
            self.test_beam_shift = self.default_test_beam_shift
        try:
            self.test_beam_shift = float(choice)
        except ValueError:
            print ('%s is not a good value for test beam separation default is used' %choice)
            self.test_beam_shift = self.default_test_beam_shift        
        
    #
    #  output
    #
    def generate_beam_cfg(self, filename, flag_write_append):
        fd = open(filename, flag_write_append)

        fd.write('[beam_info] \n')
        fd.write('particle_name : %s \n' %self.particle_name)
        fd.write('; You can use relativistic beta (betarel), relativistic gamma (gammarel) or kinetic energy (Ekin [MeV])  \n')
        fd.write('; Relativistic beta  \n')
        fd.write('betarel : %e \n' %self.betarel)
        fd.write(';gammarel : %e \n' %self.gammarel)
        fd.write(';Ekin : %e \n' %self.Ekin)        
        fd.write('; test source separation in meters \n')
        fd.write('test_beam_shift = %e \n' %self.test_beam_shift)
        fd.close()
        
            
class pyTlWall(object):
    def __init__(self, component, beam, freq, cfg_file = None):
        self._ZLong = None
        self._ZDip = None
        self._ZTrans = None
        self._ZLongDSC = None
        self._ZLongISC = None
        self._ZDipDSC  = None
        self._ZTransDSC  = None
        self._ZDipISC  = None
        self._ZTransISC  = None        
        
        self.file_ZLong = None
        self.file_ZDip = None
        self.file_ZTrans = None
        self.file_ZLongDSC = None
        self.file_ZLongISC = None
        self.file_ZDipDSC  = None
        self.file_ZTransDSC  = None
        self.file_ZDipISC  = None
        self.file_ZTransISC  = None
        self.img_ZLong = None
        self.img_ZDip = None
        self.img_ZTrans = None
        self.img_ZLongDSC = None
        self.img_ZLongISC = None
        self.img_ZDipDSC  = None
        self.img_ZTransDSC  = None
        self.img_ZDipISC  = None
        self.img_ZTransISC  = None        
        #
        #  define internal function used to calculate kprop and KZ
        #
        self.calc_kprop = self.calc_kprop_method2
        #~ self.calc_kprop = self.calc_kprop_method1_withbeta
        self.calc_KZ = self.calc_KZ_method2
        
        
        self.component = component
        self.beam = beam
        self.f = freq.freq
        self.freq = freq

        
        if cfg_file != None:
            self.read_cfg(cfg_file)

        self.corr_impedance_factor = self.calc_corr_impedance_factor()

        if self.file_ZLong != None:
            self.gen_and_save_ZLong(self.file_ZLong)
        if self.file_ZDip != None:
            self.gen_and_save_ZDip(self.file_ZDip)
        if self.file_ZTrans != None:
            self.gen_and_save_ZTrans(self.file_ZTrans)
        if self.file_ZLongDSC != None:
            self.gen_and_save_ZLongDSC(self.file_ZLongDSC)
        if self.file_ZLongISC != None:
            self.gen_and_save_ZLongISC(self.file_ZLongISC)
        if self.file_ZDipDSC != None:
            self.gen_and_save_ZDipDSC(self.file_ZDipDSC)
        if self.file_ZTransDSC != None:
            self.gen_and_save_ZTransDSC(self.file_ZTransDSC)
        if self.file_ZDipISC != None:
            self.gen_and_save_ZDipISC(self.file_ZDipISC)
        if self.file_ZTransISC != None:
            self.gen_and_save_ZTransISC(self.file_ZTransISC)
        if self.img_ZLong != None:
            self.plot_ZLong(self.img_ZLong)
        if self.img_ZDip != None:
            self.plot_ZDip(self.img_ZDip)
        if self.img_ZDipX != None:
            self.plot_ZDipX(self.img_ZDipX)
        if self.img_ZDipY!= None:
            self.plot_ZDipY(self.img_ZDipY)
        if self.img_ZQuadX != None:
            self.plot_ZQuadX(self.img_ZQuadX)
        if self.img_ZQuadY!= None:
            self.plot_ZQuadY(self.img_ZQuadY)

        if self.img_ZLongDSC != None:
            self.plot_ZLongDSC(self.img_ZLongDSC)
        if self.img_ZLongISC != None:
            self.plot_ZLongISC(self.img_ZLongISC)
        if self.img_ZDipDSC != None:
            self.plot_ZDipDSC(self.img_ZDipDSC)
        if self.img_ZTransDSC != None:
            self.plot_ZTransDSC(self.img_ZTransDSC)
        if self.img_ZDipISC != None:
            self.plot_ZDipISC(self.img_ZDipISC)
        if self.img_ZTransISC != None:
            self.plot_ZTransISC(self.img_ZTransISC)
                        
    @property
    def ZLong(self):
        """ Zlong = pipe_len  KZeff / (2 pi pipe_radius) """
        KZeff = self.calc_KZeff()
        ZLong = self.component.pipe_len * KZeff /(2. * const.pi * self.component.pipe_radius)
        ZLong_corr = ZLong / self.corr_impedance_factor
        return ZLong_corr
    
    def exist(self):
        return True
    @property
    def ZDip(self):
        """  Z_dip =  2* Zlongin * bypass /( pipe_radius**2 * beta ) """
        KZeffin  = self.calc_KZeffin()
        Zlongin = self.component.pipe_len * KZeffin /(2 * const.pi * self.component.pipe_radius)
        
        beta = 2 * const.pi * self.f * np.sqrt(const.epsilon_0 * const.mu_0) / self.beam.betarel
        if ( self.component.layers[-1]['type'].upper() == 'V' or self.component.layers[-1]['type'].upper() == 'PEC') :
            Ind = const.mu_0 / 2
        else:
            Ind = const.mu_0
            
        
        Zind = 1.j * self.f * (Ind ) * self.component.pipe_len 
        bypass = Zind / (Zlongin + Zind) 
        
        ZDip =  2* Zlongin * bypass /( self.component.pipe_radius**2 * beta )
        ZDip_corr = ZDip / self.corr_impedance_factor
        return ZDip_corr
    
    #
    #  Space charge functions
    #
    @property
    def ZLongDSC(self):
        if self.beam.gammarel == float('inf'):
            print ('Longitudinal direct space charge not calculated ')
            return   np.ones(len(self.f) )* float('nan') + 1.j * np.ones(len(self.f) )* float('nan')
        kbess = 2 * const.pi * self.f / (self.beam.betarel * const.c)
        argbess0 = kbess * self.beam.test_beam_shift/ self.beam.gammarel

        BessBeamL = (i0(argbess0))**2    
        BessBeamLDSC = k0(argbess0) / i0(argbess0)

        # longitudinal direct space charge
        Z_long_DSC = - 1.j * const.pi * self.f * self.component.pipe_len * BessBeamL  * BessBeamLDSC / \
                            (const.epsilon_0 * const.pi * (self.beam.gammarel * self.beam.betarel * const.c)**2 )
        return Z_long_DSC

    @property
    def ZLongISC(self):
        if self.beam.gammarel == float('inf'):
            print ('Longitudinal indirect space charge not calculated ')
            return   np.ones(len(self.f) )* float('nan') + 1.j * np.ones(len(self.f) )* float('nan')
            
        kbess = 2 * const.pi * self.f / (self.beam.betarel * const.c)
        argbess0 = kbess * self.beam.test_beam_shift/ self.beam.gammarel
        argbess1 = kbess * self.component.pipe_radius / self.beam.gammarel

        BessBeamL = (i0(argbess0))**2
        BessBeamLISC = - k0(argbess1) / i0(argbess1)
        # longitudinal indirect space charge
        Z_long_ISC = - 1.j * const.pi * self.f * self.component.pipe_len * BessBeamL  * BessBeamLISC / \
                            (const.epsilon_0 * const.pi * (self.beam.gammarel * self.beam.betarel * const.c)**2 )
        
        return Z_long_ISC  
        
    @property
    def ZDipISC(self):
        if self.beam.gammarel == float('inf'):
            print ('Dipolar indirect space charge not calculated ')
            return   np.ones(len(self.f) )* float('nan') + 1.j * np.ones(len(self.f) )* float('nan')
        kbess = 2 * const.pi * self.f / (self.beam.betarel * const.c)

        argbess0 = kbess * self.beam.test_beam_shift/ self.beam.gammarel
        argbess1 = kbess * self.component.pipe_radius / self.beam.gammarel

        BessBeamTISC = - k1(argbess1) / i1(argbess1)

        # transverse indirect space charge
        if ( self.beam.test_beam_shift == 0. ) :    
            BessA = kbess * kbess / (2 * self.beam.gammarel * self.beam.gammarel)  \
                    *  k1( argbess1) / i1(argbess1 )
            Z_trans_ISC = 1.j * const.physical_constants['characteristic impedance of vacuum'][0] * self.component.pipe_len * \
                          BessA / (2 * const.pi * self.beam.gammarel * self.beam.gammarel * self.beam.betarel )
        else:
            BessBeamT = (i1(argbess0) / self.beam.test_beam_shift)**2
            Z_trans_ISC = - 1.j * const.physical_constants['characteristic impedance of vacuum'][0] * self.component.pipe_len * \
                           BessBeamT  * BessBeamTISC / (const.pi * self.beam.gammarel**2 * self.beam.betarel)
                     
        return Z_trans_ISC
    
    @property
    def ZDipDSC(self):
        if self.beam.gammarel == float('inf'):
            print ('Dipolar direct space charge not calculated ')
            return   np.ones(len(self.f) )* float('nan') + 1.j * np.ones(len(self.f) )* float('nan')
            
        kbess = 2 * const.pi * self.f / (self.beam.betarel * const.c)
        argbess0 = kbess * self.beam.test_beam_shift/ self.beam.gammarel
        argbess1 = kbess * self.component.pipe_radius / self.beam.gammarel

        BessBeamT = (i1(argbess0) / self.beam.test_beam_shift)**2
        
        BessBeamTDSC = k1(argbess0) / i1(argbess0)

        Z_trans_DSC = - 1.j * const.physical_constants['characteristic impedance of vacuum'][0] * self.component.pipe_len * BessBeamT  * \
                      BessBeamTDSC / (const.pi * self.beam.gammarel**2 * self.beam.betarel)

        return Z_trans_DSC
    
    #
    #  Impedance functions
    #
    def calc_corr_impedance(self, Z):
        """ Z_corr = Z/ reduct_factor  """
        index = iv(1, 2 * const.pi * self.f * self.component.pipe_radius /(self.beam.betarel * const.c * self.beam.gammarel)) > accuracy_factor * iv(0, 2 * const.pi * self.f * self.component.pipe_radius /(self.beam.betarel * const.c * self.beam.gammarel))
        if len(self.f[index]) > 0.:
            print ("Warning: results not accured for f > %.2e MHz" %(min(self.f[index])*1.e-6))
        reduct_factor =  iv(0, 2 * const.pi * self.f * self.component.pipe_radius /(self.beam.betarel * const.c * self.beam.gammarel))
        Z_corr = Z/ reduct_factor

        return Z_corr
    
    def calc_corr_impedance_factor(self):
        """ Z_corr = Z/ reduct_factor  """
        index = iv(1, 2 * const.pi * self.f * self.component.pipe_radius /(self.beam.betarel * const.c * self.beam.gammarel)) > accuracy_factor * iv(0, 2 * const.pi * self.f * self.component.pipe_radius /(self.beam.betarel * const.c * self.beam.gammarel))
        if len(self.f[index]) > 0.:
            print ("Warning: results not accured for f > %.2e MHz" %(min(self.f[index])*1.e-6))
        reduct_factor =  iv(0, 2 * const.pi * self.f * self.component.pipe_radius /(self.beam.betarel * const.c * self.beam.gammarel))
    
        return reduct_factor

        
        
        
    #
    #   KZeff and KZeffin functions
    #
    def calc_KZeff(self):
        """ KZeff is calculated recursively with the formula
        KZeff = KZ ( (KZeff  + j KZ  tan(kprop t)) / (KZ  + 1.j  KZeff  tan(kprop t) )
        for boundary 
            if PEC   KZeff = 0
            elif V   KZeff = const.physical_constants['characteristic impedance of vacuum'][0] (1 - (1 /BesselI(0,kprop * pipe_radius )))
            else     KZeff = sqrt(mu/eps)
            """
        if self.component.layers[-1]['type'].upper() == 'PEC':
            KZeff = 0
        elif  self.component.layers[-1]['type'].upper() == 'V':
            #kprop = self.calc_kprop(self.f, 0., 0. ,1., float('inf'), 1.)
            kprop =  2 * const.pi * self.f / (  const.c) 
            Scil = (1 /iv(0,kprop * self.component.pipe_radius ))
            KZeff = const.physical_constants['characteristic impedance of vacuum'][0] * (1 -  Scil) 
        elif  self.component.layers[-1]['type'].upper() == 'PMC':
            print ('Error not implemented method')
        else:
            sigmaAC =  self.component.calc_sigmaAC(self.f, self.component.layers[-1]['sigmaDC'], self.component.layers[-1]['tau'])
            eps = self.component.calc_eps(self.f, sigmaAC, self.component.layers[-1]['epsr'])
            mu  = const.mu_0 * self.component.calc_mur(self.f, self.component.layers[-1]['k'], self.component.layers[-1]['muinf'])
            KZ = np.sqrt( mu / eps )
            #~ kprop = self.calc_kprop(self.f, self.component.layers[-1]['sigmaDC'], self.component.layers[-1]['tau'], self.component.layers[-1]['epsr'], 
                                   #~ self.component.layers[-1]['k'], self.component.layers[-1]['muinf'])

            #~ Scil = (1 /iv(0,np.abs(kprop) * self.component.pipe_radius ))
            #~ KZeff = KZ * (1 - Scil)
            KZeff = KZ 

        for i in range(len(self.component.layers)-2 , -1 ,-1 ):
            if self.component.layers[i]['type'].upper() == 'PEC':
                KZ = 0
            elif self.component.layers[i]['type'].upper() == 'V':
                KZ = const.physical_constants['characteristic impedance of vacuum'][0]
                #kprop = self.calc_kprop(self.f, 0., 0. ,1., float('inf'), 1.)
                kprop = 2 * const.pi * self.f / const.c 
            elif  self.component.layers[-1]['type'].upper() == 'PMC':
                print ('Error not implemented method')
            else: 
                KZ = self.calc_KZ( self.f, self.component.layers[i]['sigmaDC'], self.component.layers[i]['tau'], self.component.layers[i]['epsr'], 
                                   self.component.layers[i]['k'], self.component.layers[i]['muinf'])
                kprop = self.calc_kprop(self.f, self.component.layers[i]['sigmaDC'], self.component.layers[i]['tau'], self.component.layers[i]['epsr'], 
                                   self.component.layers[i]['k'], self.component.layers[i]['muinf'])
                
            try:
                tan_t_kprop = np.array(list(map(lambda  currkprop: cmath.tan( currkprop  * self.component.layers[i]['thick']), kprop)))
            except  TypeError:
                tan_t_kprop = cmath.tan( kprop  * self.component.layers[i]['thick'])
                
            KZeff = KZ * ( (KZeff  + 1.j * KZ * tan_t_kprop ) / (KZ  + 1.j * KZeff  * tan_t_kprop) )
            

        return KZeff
        
    def calc_KZeffin(self):
        """ KZeffin is calculated recursively with the formula
        KZeffin = KZ ( (KZeff  + j KZ  tan(kprop t)) / (KZ  + 1.j  KZeffin  tan(kprop t) )
        for boundary 
            if PEC   KZeff = 0
            elif V   KZeff = const.physical_constants['characteristic impedance of vacuum'][0] 
            else     KZeff = sqrt(mu/eps)
            """
        
        # boundary layer
        if self.component.layers[-1]['type'].upper() == 'PEC':
            KZeffin = 0
        elif  self.component.layers[-1]['type'].upper() == 'V':
            KZeffin = const.physical_constants['characteristic impedance of vacuum'][0] 
        elif  self.component.layers[-1]['type'].upper() == 'PMC':
            print ('Error not implemented method')
        else:
            sigmaAC =  self.component.calc_sigmaAC(self.f, self.component.layers[-1]['sigmaDC'], self.component.layers[-1]['tau'])
            eps = self.component.calc_eps(self.f, sigmaAC, self.component.layers[-1]['epsr'])
            mu  = const.mu_0 * self.component.calc_mur(self.f, self.component.layers[-1]['k'], self.component.layers[-1]['muinf'])
            KZ = np.sqrt( mu / eps )
            KZeffin = KZ
        
        for i in range(len(self.component.layers)-2 , -1 ,-1 ):
            if self.component.layers[i]['type'].upper() == 'PEC':
                KZ = 0
            elif self.component.layers[i]['type'].upper() == 'V':
                KZ = const.physical_constants['characteristic impedance of vacuum'][0]
                #kprop = self.calc_kprop(self.f, 0., 0. ,1., float('inf'), 1.)
                kprop = 2 * const.pi * self.f / const.c 
            elif  self.component.layers[-1]['type'].upper() == 'PMC':
                print ('Error not implemented method')
            else: 
                KZ = self.calc_KZ( self.f, self.component.layers[i]['sigmaDC'], self.component.layers[i]['tau'], self.component.layers[i]['epsr'], 
                                   self.component.layers[i]['k'], self.component.layers[i]['muinf'])
                kprop = self.calc_kprop(self.f, self.component.layers[i]['sigmaDC'], self.component.layers[i]['tau'], self.component.layers[i]['epsr'], 
                                   self.component.layers[i]['k'], self.component.layers[i]['muinf'])

            try:
                tan_t_kprop = np.array(list(map(lambda  currkprop: cmath.tan( currkprop  * self.component.layers[i]['thick']), kprop)))
            except  TypeError:
                tan_t_kprop = cmath.tan( kprop  * self.component.layers[i]['thick'])
                
            KZeffin = KZ * ( (KZeffin  + 1.j * KZ * tan_t_kprop ) / (KZ  + 1.j * KZeffin  * tan_t_kprop) )
            
        return     KZeffin 

    #
    #  Single layer function
    #
    def calc_KZ_method1(self, f, sigmaDC, tau, epsr, k, muinf):
        """ KZ = sqrt(mu / eps) """
        sigma_ac = self.component.calc_sigmaAC(f, sigmaDC, tau)
        eps = self.component.calc_eps( f, sigma_ac, epsr)
        mur =  self.component.calc_mur(f, k, muinf)
        mu = const.mu_0 * mur
        KZ = np.sqrt( mu / eps  )    
        
        return KZ

    def calc_KZ_method2(self,  f, sigmaDC, tau, epsr, k, muinf ):
        """   KZ = (1 + j) / (sigmapm + deltam) """
        eps = const.epsilon_0 * epsr
        sigma_ac = self.component.calc_sigmaAC(f, sigmaDC, tau)
        mur =  self.component.calc_mur(f, k, muinf)
        mu = const.mu_0 * mur
        
        sigmapm = self.component.calc_sigmaPM( f, eps, sigma_ac)
        deltam = self.component.calc_deltaM(f, mu, eps, sigma_ac)
        
        KZ = (1. + 1.j) / (sigmapm * deltam)
        return KZ
        
    def calc_kprop_method1_withbeta(self, f, sigmaDC, tau , epsr, k, muinf):
        """ kprop = 2 pi f sqrt( mu eps - (betarel * c)^-2) """
        sigma_ac = self.component.calc_sigmaAC(f, sigmaDC, tau)
        eps = self.component.calc_eps(f,  sigma_ac, epsr)
        mur =  self.component.calc_mur(f, k, muinf)
        mu = const.mu_0 * mur
        kprop = 2  * const.pi * f * np.sqrt((mu * eps)  - (self.beam.betarel  * const.c )**(-2) ) 
        return kprop

    def calc_kprop_method1_nobeta(self, f,sigmaDC, tau , epsr, k, muinf):
        """ kprop = 2 pi f sqrt( const.mu_0 mur eps ) """
        sigma_ac = self.component.calc_sigmaAC(f, sigmaDC, tau)
        eps = self.component.calc_eps( f, sigma_ac, epsr)
        mur =  self.component.calc_mur(f, k, muinf)
        mu = const.mu_0 * mur
        kprop = 2 * const.pi * f * np.sqrt(mu * eps  ) 
        return kprop

    def calc_kprop_method2(self, f, sigmaDC, tau, epsr, k, muinf):
        """  kprop = (1 - 1.j) / delta """
        sigma_ac = self.component.calc_sigmaAC(f, sigmaDC, tau)
        eps = const.epsilon_0 * epsr
        mur =  self.component.calc_mur(f, k, muinf)
        mu = const.mu_0 * mur
        delta = self.component.calc_delta(f, mu, eps, sigma_ac)
        kprop =   (1 - 1.j ) / delta 
        return kprop
        
        
    #
    #  Input functions
    #
    def ask_impedance_calc(self):
        choice = ''
        while choice.upper() != 'X':
            print ('What do you want to calculate? ')
            print ('-- ZLong -- Longitudinal impedance')
            print ('-- ZDip -- Transverse impedance considering round chamber and no twiss beta')
            print ('-- ZTrans -- Dipolar and Quadrupolar impedance in X and Y direction taking care of the %s shape and betatron functions betax= %.3f [m] betay= %.3f [m]'%(self.component.form_factor, self.component.betax, self.component.betay))
            print ('-- ZLongISC -- Longitudinal indirect space charge impedance')
            print ('-- ZLongDSC -- Longitudinal direct space charge impedance')
            print ('-- ZDipISC -- Transverse indirect space charge impedance considering round chamber and no twiss beta')
            print ('-- ZTransISC -- Dipolar and Quadrupolar indirect space charge impedance in X and Y direction taking care of the %s shape and betatron functions betax= %.3f [m] betay= %.3f [m]'%(self.component.form_factor, self.component.betax, self.component.betay))
            print ('-- ZDipDSC -- Transverse direct space charge impedance considering round chamber and no twiss beta')
            print ('-- ZTransDSC -- Dipolar and Quadrupolar direct space charge impedance in X and Y direction taking care of the %s shape and betatron functions betax= %.3f [m] betay= %.3f [m]'%(self.component.form_factor, self.component.betax, self.component.betay))
            choice = raw_input('digit x to close this menu  ')
            
            if choice.lower() == 'zlong':
                self.sav_ZLong = self.ZLong
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZLong = choice
                    self.save_ZLong(choice, self.sav_ZLong, 'ZLong')
                
            if choice.lower() == 'zdip':
                self.sav_ZDip = self.ZDip
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZDip = choice
                    self.save_ZDip(choice, self.sav_ZDip, 'ZDip')

            if choice.lower() == 'ztrans':
                self.sav_ZDipX = self.ZDip * self.component._dipx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZDipY = self.ZDip * self.component._dipy_form_factor[self.component.form_factor] * self.component.betay
                self.sav_ZQuadX = self.ZDip * self.component._quadx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZQuadY = self.ZDip * self.component._quady_form_factor[self.component.form_factor] * self.component.betay

                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZTrans = choice
                    self.save_ZTrans(choice, self.sav_ZDipX,self.sav_ZDipY , self.sav_ZQuadX,self.sav_ZQuadY,  'ZTrans')
                    
            if choice.lower() == 'zlongisc':
                self.sav_ZLongISC = self.ZLongISC
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZLongISC = choice
                    self.save_ZLong(choice, self.sav_ZLongISC, 'ZLongISC')
                
            if choice.lower() == 'zdipisc':
                self.sav_ZDipISC = self.ZDipISC
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZDipISC = choice
                    self.save_ZDip(choice, self.sav_ZDip, 'ZDipISC')

            if choice.lower() == 'ztransisc':
                self.sav_ZDipXISC = self.ZDipISC * self.component._dipx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZDipYISC = self.ZDipISC * self.component._dipy_form_factor[self.component.form_factor] * self.component.betay
                self.sav_ZQuadXISC = self.ZDipISC * self.component._quadx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZQuadYISC = self.ZDipISC * self.component._quady_form_factor[self.component.form_factor] * self.component.betay

                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZTransISC = choice
                    self.save_ZTrans(choice, self.sav_ZDipXISC,self.sav_ZDipYISC , self.sav_ZQuadXISC,self.sav_ZQuadYISC,  'ZTransISC')
                    
            if choice.lower() == 'zlongdsc':
                self.sav_ZLongDSC = self.ZLongDSC
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZLongDSC = choice
                    self.save_ZLong(choice, self.sav_ZLongDSC, 'ZLongDSC')
                
            if choice.lower() == 'zdipdsc':
                self.sav_ZDipDSC = self.ZDipDSC
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZDipDSC = choice
                    self.save_ZDip(choice, self.sav_ZDip, 'ZDipDSC')

            if choice.lower() == 'ztransdsc':
                self.sav_ZDipXDSC = self.ZDipDSC * self.component._dipx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZDipYDSC = self.ZDipDSC * self.component._dipy_form_factor[self.component.form_factor] * self.component.betay
                self.sav_ZQuadXDSC = self.ZDipDSC * self.component._quadx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZQuadYDSC = self.ZDipDSC * self.component._quady_form_factor[self.component.form_factor] * self.component.betay

                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZTransDSC = choice
                    self.save_ZTrans(choice, self.sav_ZDipXDSC,self.sav_ZDipYDSC , self.sav_ZQuadXDSC,self.sav_ZQuadYDSC,  'ZTransDSC')
                    
    def ask_impedance_plot(self):
        choice = ''
        while choice.upper() != 'X':
            print ('What do you want to plot? ')
            print ('-- ZLong ')
            print ('-- ZDip ')
            print ('-- ZTrans ')
            print ('-- ZLongISC ')
            print ('-- ZLongDSC ')
            print ('-- ZDipISC ')
            print ('-- ZTransISC  ')
            print ('-- ZDipDSC ')
            print ('-- ZTransDSC ')
            choice = raw_input('digit x to close this menu  ')
            
            if choice.lower() == 'zlong':
                self.sav_ZLong = self.ZLong
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZLong = choice
                    self.save_ZLong(choice, self.sav_ZLong, 'ZLong')
                
            if choice.lower() == 'zdip':
                self.sav_ZDip = self.ZDip
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZDip = choice
                    self.save_ZDip(choice, self.sav_ZDip, 'ZDip')

            if choice.lower() == 'ztrans':
                self.sav_ZDipX = self.ZDip * self.component._dipx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZDipY = self.ZDip * self.component._dipy_form_factor[self.component.form_factor] * self.component.betay
                self.sav_ZQuadX = self.ZDip * self.component._quadx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZQuadY = self.ZDip * self.component._quady_form_factor[self.component.form_factor] * self.component.betay

                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZTrans = choice
                    self.save_ZTrans(choice, self.sav_ZDipX,self.sav_ZDipY , self.sav_ZQuadX,self.sav_ZQuadY,  'ZTrans')
                    
            if choice.lower() == 'zlongisc':
                self.sav_ZLongISC = self.ZLongISC
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZLongISC = choice
                    self.save_ZLong(choice, self.sav_ZLongISC, 'ZLongISC')
                
            if choice.lower() == 'zdipisc':
                self.sav_ZDipISC = self.ZDipISC
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZDipISC = choice
                    self.save_ZDip(choice, self.sav_ZDip, 'ZDipISC')

            if choice.lower() == 'ztransisc':
                self.sav_ZDipXISC = self.ZDipISC * self.component._dipx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZDipYISC = self.ZDipISC * self.component._dipy_form_factor[self.component.form_factor] * self.component.betay
                self.sav_ZQuadXISC = self.ZDipISC * self.component._quadx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZQuadYISC = self.ZDipISC * self.component._quady_form_factor[self.component.form_factor] * self.component.betay

                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZTransISC = choice
                    self.save_ZTrans(choice, self.sav_ZDipXISC,self.sav_ZDipYISC , self.sav_ZQuadXISC,self.sav_ZQuadYISC,  'ZTransISC')
                    
            if choice.lower() == 'zlongdsc':
                self.sav_ZLongDSC = self.ZLongDSC
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZLongDSC = choice
                    self.save_ZLong(choice, self.sav_ZLongDSC, 'ZLongDSC')
                
            if choice.lower() == 'zdipdsc':
                self.sav_ZDipDSC = self.ZDipDSC
                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZDipDSC = choice
                    self.save_ZDip(choice, self.sav_ZDip, 'ZDipDSC')

            if choice.lower() == 'ztransdsc':
                self.sav_ZDipXDSC = self.ZDipDSC * self.component._dipx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZDipYDSC = self.ZDipDSC * self.component._dipy_form_factor[self.component.form_factor] * self.component.betay
                self.sav_ZQuadXDSC = self.ZDipDSC * self.component._quadx_form_factor[self.component.form_factor] * self.component.betax
                self.sav_ZQuadYDSC = self.ZDipDSC * self.component._quady_form_factor[self.component.form_factor] * self.component.betay

                choice = raw_input('Digit choice to save the data (ENTER to not save) ')
                if choice != '':
                    self.file_ZTransDSC = choice
                    self.save_ZTrans(choice, self.sav_ZDipXDSC,self.sav_ZDipYDSC , self.sav_ZQuadXDSC,self.sav_ZQuadYDSC,  'ZTransDSC')
    def read_cfg(self, cfg_file):

        config = configparser.ConfigParser()
        config.readfp(open(cfg_file))
        
        # Component informations
        self.component.read_acc_component_cfg(config)
        
        # Beam informations
        self.beam.read_beam_cfg(config)

        # frequency informations
        self.freq.read_freq_cfg(config)
        self.f = self.freq.freq
        
        #
        #  Outputs Information
        #
        try:
            self.file_ZLong = config.get('output', 'ZLong')
        except configparser.NoOptionError:
            self.file_ZLong = None
        try:
            self.file_ZDip = config.get('output', 'ZDip')
        except configparser.NoOptionError:
            self.file_ZDip = None
        try:
            self.file_ZTrans = config.get('output', 'ZTrans')
        except configparser.NoOptionError:
            self.file_ZTrans = None
        try:
            self.file_ZLongDSC = config.get('output', 'ZLongDSC')
        except configparser.NoOptionError:
            self.file_ZLongDSC = None
        try:
            self.file_ZLongISC = config.get('output', 'ZLongISC')
        except configparser.NoOptionError:
            self.file_ZLongISC = None
        try:
            self.file_ZDipDSC = config.get('output', 'ZDipDSC')
        except configparser.NoOptionError:
            self.file_ZDipDSC = None
        try:
            self.file_ZTransDSC = config.get('output', 'ZTransDSC')
        except configparser.NoOptionError:
            self.file_ZTransDSC = None
        try:
            self.file_ZDipISC = config.get('output', 'ZDipISC')
        except configparser.NoOptionError:
            self.file_ZDipISC = None
        try:
            self.file_ZTransISC = config.get('output', 'ZTransISC')
        except configparser.NoOptionError:
            self.file_ZTransISC = None

        #
        #  Image Information
        #
        if config.has_section('graphical_output'):
            
            try:
                self.img_ZLong = config.get('graphical_output', 'ZLong')
            except configparser.NoOptionError:
                self.img_ZLong = None
            if self.img_ZLong != None: 
                try:
                    self.img_ZLong_title = config.get('graphical_output', 'ZLong_title')
                except configparser.NoOptionError:
                    self.img_ZLong_title = None            
                try:
                    self.img_ZLong_freq_scale = config.get('graphical_output', 'ZLong_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZLong_freq_scale = None            
                try:
                    self.img_ZLong_scale = config.get('graphical_output', 'ZLong_scale')
                except configparser.NoOptionError:
                    self.img_ZLong_scale = None            
            try:
                self.img_ZDip = config.get('graphical_output', 'ZDip')
            except configparser.NoOptionError:
                self.img_ZDip = None
            if self.img_ZDip != None: 
                try:
                    self.img_ZDip_title = config.get('graphical_output', 'ZDip_title')
                except configparser.NoOptionError:
                    self.img_ZDip_title = None            
                try:
                    self.img_ZDip_freq_scale = config.get('graphical_output', 'ZDip_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZDip_freq_scale = None            
                try:
                    self.img_ZDip_scale = config.get('graphical_output', 'ZDip_scale')
                except configparser.NoOptionError:
                    self.img_ZDip_scale = None            
            try:
                self.img_ZDipX = config.get('graphical_output', 'ZDipX')
            except configparser.NoOptionError:
                self.img_ZDipX = None
            if self.img_ZDipX != None: 
                try:
                    self.img_ZDipX_title = config.get('graphical_output', 'ZDipX_title')
                except configparser.NoOptionError:
                    self.img_ZDipX_title = None            
                try:
                    self.img_ZDipX_freq_scale = config.get('graphical_output', 'ZDipX_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZDipX_freq_scale = None            
                try:
                    self.img_ZDipX_scale = config.get('graphical_output', 'ZDipX_scale')
                except configparser.NoOptionError:
                    self.img_ZDipX_scale = None        
            try:
                self.img_ZDipY = config.get('graphical_output', 'ZDipY')
            except configparser.NoOptionError:
                self.img_ZDipY = None
            if self.img_ZDipY != None: 
                try:
                    self.img_ZDipY_title = config.get('graphical_output', 'ZDipY_title')
                except configparser.NoOptionError:
                    self.img_ZDipY_title = None            
                try:
                    self.img_ZDipY_freq_scale = config.get('graphical_output', 'ZDipY_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZDipY_freq_scale = None            
                try:
                    self.img_ZDipY_scale = config.get('graphical_output', 'ZDipY_scale')
                except configparser.NoOptionError:
                    self.img_ZDipY_scale = None        
            try:
                self.img_ZQuadX = config.get('graphical_output', 'ZQuadX')
            except configparser.NoOptionError:
                self.img_ZQuadX = None
            if self.img_ZQuadX != None: 
                try:
                    self.img_ZQuadX_title = config.get('graphical_output', 'ZQuadX_title')
                except configparser.NoOptionError:
                    self.img_ZQuadX_title = None            
                try:
                    self.img_ZQuadX_freq_scale = config.get('graphical_output', 'ZQuadX_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZQuadX_freq_scale = None            
                try:
                    self.img_ZQuadX_scale = config.get('graphical_output', 'ZQuadX_scale')
                except configparser.NoOptionError:
                    self.img_ZQuadX_scale = None        
            try:
                self.img_ZQuadY = config.get('graphical_output', 'ZQuadY')
            except configparser.NoOptionError:
                self.img_ZQuadY = None
            if self.img_ZQuadY != None: 
                try:
                    self.img_ZQuadY_title = config.get('graphical_output', 'ZQuadY_title')
                except configparser.NoOptionError:
                    self.img_ZQuadY_title = None            
                try:
                    self.img_ZQuadY_freq_scale = config.get('graphical_output', 'ZQuadY_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZQuadY_freq_scale = None            
                try:
                    self.img_ZQuadY_scale = config.get('graphical_output', 'ZQuadY_scale')
                except configparser.NoOptionError:
                    self.img_ZQuadY_scale = None        

            try:
                self.img_ZLongISC = config.get('graphical_output', 'ZLongISC')
            except configparser.NoOptionError:
                self.img_ZLongISC = None
            if self.img_ZLongISC != None: 
                try:
                    self.img_ZLongISC_title = config.get('graphical_output', 'ZLongISC_title')
                except configparser.NoOptionError:
                    self.img_ZLongISC_title = None            
                try:
                    self.img_ZLongISC_freq_scale = config.get('graphical_output', 'ZLongISC_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZLongISC_freq_scale = None            
                try:
                    self.img_ZLongISC_scale = config.get('graphical_output', 'ZLongISC_scale')
                except configparser.NoOptionError:
                    self.img_ZLongISC_scale = None            
            try:
                self.img_ZLongDSC = config.get('graphical_output', 'ZLongDSC')
            except configparser.NoOptionError:
                self.img_ZLongDSC = None
            if self.img_ZLongDSC != None: 
                try:
                    self.img_ZLongDSC_title = config.get('graphical_output', 'ZLongDSC_title')
                except configparser.NoOptionError:
                    self.img_ZLongDSC_title = None            
                try:
                    self.img_ZLongDSC_freq_scale = config.get('graphical_output', 'ZLongDSC_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZLongDSC_freq_scale = None            
                try:
                    self.img_ZLongDSC_scale = config.get('graphical_output', 'ZLongDSC_scale')
                except configparser.NoOptionError:
                    self.img_ZLongDSC_scale = None                            
            try:
                self.img_ZDipISC = config.get('graphical_output', 'ZDipISC')
            except configparser.NoOptionError:
                self.img_ZDipISC = None
            if self.img_ZDipISC != None: 
                try:
                    self.img_ZDipISC_title = config.get('graphical_output', 'ZDipISC_title')
                except configparser.NoOptionError:
                    self.img_ZDipISC_title = None            
                try:
                    self.img_ZDipISC_freq_scale = config.get('graphical_output', 'ZDipISC_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZDipISC_freq_scale = None            
                try:
                    self.img_ZDipISC_scale = config.get('graphical_output', 'ZDipISC_scale')
                except configparser.NoOptionError:
                    self.img_ZDipISC_scale = None            
            try:
                self.img_ZDipDSC = config.get('graphical_output', 'ZDipDSC')
            except configparser.NoOptionError:
                self.img_ZDipDSC = None
            if self.img_ZDipDSC != None: 
                try:
                    self.img_ZDipDSC_title = config.get('graphical_output', 'ZDipDSC_title')
                except configparser.NoOptionError:
                    self.img_ZDipDSC_title = None            
                try:
                    self.img_ZDipDSC_freq_scale = config.get('graphical_output', 'ZDipDSC_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZDipDSC_freq_scale = None            
                try:
                    self.img_ZDipDSC_scale = config.get('graphical_output', 'ZDipDSC_scale')
                except configparser.NoOptionError:
                    self.img_ZDipDSC_scale = None
            try:
                self.img_ZTransISC = config.get('graphical_output', 'ZTransISC')
            except configparser.NoOptionError:
                self.img_ZTransISC = None
            if self.img_ZTransISC != None: 
                try:
                    self.img_ZTransISC_title = config.get('graphical_output', 'ZTransISC_title')
                except configparser.NoOptionError:
                    self.img_ZTransISC_title = None            
                try:
                    self.img_ZTransISC_freq_scale = config.get('graphical_output', 'ZTransISC_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZTransISC_freq_scale = None            
                try:
                    self.img_ZTransISC_scale = config.get('graphical_output', 'ZTransISC_scale')
                except configparser.NoOptionError:
                    self.img_ZTransISC_scale = None            
            try:
                self.img_ZTransDSC = config.get('graphical_output', 'ZTransDSC')
            except configparser.NoOptionError:
                self.img_ZTransDSC = None
            if self.img_ZTransDSC != None: 
                try:
                    self.img_ZTransDSC_title = config.get('graphical_output', 'ZTransDSC_title')
                except configparser.NoOptionError:
                    self.img_ZTransDSC_title = None            
                try:
                    self.img_ZTransDSC_freq_scale = config.get('graphical_output', 'ZTransDSC_freq_scale')
                except configparser.NoOptionError:
                    self.img_ZTransDSC_freq_scale = None            
                try:
                    self.img_ZTransDSC_scale = config.get('graphical_output', 'ZTransDSC_scale')
                except configparser.NoOptionError:
                    self.img_ZTransDSC_scale = None
    #
    #  Output functions
    #
    def generate_cfg(self, filename):
        self.component.generate_acc_component_cfg(filename, 'w')
        self.freq.generate_freq_cfg(filename, 'a')
        self.beam.generate_beam_cfg(filename, 'a')
        self.generate_output_cfg(filename, 'a')
    
    def generate_output_cfg(self, filename, flag_write_append):
        fd = open(filename, flag_write_append)

        fd.write('[output] \n')
        if self.file_ZLong != None:
            fd.write('ZLong = %s \n' %self.file_ZLong)
        if self.file_ZDip != None:
            fd.write('ZDip = %s \n' %self.file_ZDip)
        if self.file_ZTrans != None:
            fd.write('ZTrans = %s \n' %self.file_ZTrans)
        if self.file_ZLongISC != None:
            fd.write('ZLongISC = %s \n' %self.file_ZLongISC)
        if self.file_ZDipISC != None:
            fd.write('ZDipISC = %s \n' %self.file_ZDipISC)
        if self.file_ZTransISC != None:
            fd.write('ZTransISC = %s \n' %self.file_ZTransISC)
        if self.file_ZLongDSC != None: 
            fd.write('ZLongDSC = %s \n' %self.file_ZLongDSC)
        if self.file_ZDipDSC != None:
            fd.write('ZDipDSC = %s \n' %self.file_ZDipDSC)
        if self.file_ZTransDSC != None:
            fd.write('ZTransDSC = %s \n' %self.file_ZTransDSC)
        
        
    def gen_and_save_ZLong(self,filename):
        ZLong = self.ZLong
        self.save_ZLong(filename, ZLong, 'ZLong')
        

    def gen_and_save_ZDip(self,filename):
        ZDip = self.ZDip
        self.save_ZDip(filename, ZDip, 'ZDip')

    def gen_and_save_ZTrans(self,filename):
        self.ZDipX = self.ZDip * self.component._dipx_form_factor[self.component.form_factor] * self.component.betax
        self.ZDipY = self.ZDip * self.component._dipy_form_factor[self.component.form_factor] * self.component.betay
        self.ZQuadX = self.ZDip * self.component._quadx_form_factor[self.component.form_factor] * self.component.betax
        self.ZQuadY = self.ZDip * self.component._quady_form_factor[self.component.form_factor] * self.component.betay
        self.save_ZTrans(filename, self.ZDipX, self.ZDipY,self.ZQuadX, self.ZQuadY , 'ZTrans')

    def gen_and_save_ZLongISC(self,filename):
        ZLongISC = self.ZLongISC
        self.save_ZLong(filename, ZLongISC, 'ZLongISC')

    def gen_and_save_ZDipISC(self,filename):
        ZDipISC = self.ZDipISC
        self.save_ZDip(filename, ZDipISC, 'ZDipISC')

    def gen_and_save_ZTransISC(self,filename):
        ZDipXISC = self.ZDipISC * self.component._dipx_form_factor[self.component.form_factor] * self.component.betax
        ZDipYISC = self.ZDipISC * self.component._dipy_form_factor[self.component.form_factor] * self.component.betay
        ZQuadXISC = self.ZDipISC * self.component._quadx_form_factor[self.component.form_factor] * self.component.betax
        ZQuadYISC = self.ZDipISC * self.component._quady_form_factor[self.component.form_factor] * self.component.betay
        self.save_ZTrans(filename, ZDipXISC, ZDipYISC,ZQuadXISC, ZQuadYISC , 'ZTransISC')
        
    def gen_and_save_ZLongDSC(self,filename):
        ZLongDSC = self.ZLongDSC
        self.save_ZLong(filename, ZLongDSC, 'ZLongDSC')

    def gen_and_save_ZDipDSC(self,filename):
        ZDipDSC = self.ZDipDSC
        self.save_ZDip(filename, ZDipDSC, 'ZDipDSC')

    def gen_and_save_ZTransDSC(self,filename):
        ZDipXDSC = self.ZDipDSC * self.component._dipx_form_factor[self.component.form_factor] * self.component.betax
        ZDipYDSC = self.ZDipDSC * self.component._dipy_form_factor[self.component.form_factor] * self.component.betay
        ZQuadXDSC = self.ZDipDSC * self.component._quadx_form_factor[self.component.form_factor] * self.component.betax
        ZQuadYDSC = self.ZDipDSC * self.component._quady_form_factor[self.component.form_factor] * self.component.betay
        self.save_ZTrans(filename, ZDipXDSC, ZDipYDSC,ZQuadXDSC, ZQuadYDSC , 'ZTransDSC')
        
    def save_ZLong(self,filename, ZLong, label):
        print ('Saving %s data in %s' %(label ,filename))
        try: 
            fd = open(filename, 'w')
        except IOError:
            createdir(filename)
            fd = open(filename, 'w')
        fd.write('f           ZLong.real            ZLong.imaginary \n')
        fd.write('[Hz]           [Ohm]            [Ohm] \n')
        for i in range(len(self.f)):
            fd.write('%e           %e           %e \n' %(self.f[i], ZLong[i].real, ZLong[i].imag))
        fd.close()

    def save_ZDip(self,filename, ZDip, label):
        print ('Saving %s data in %s' %(label, filename))
        try: 
            fd = open(filename, 'w')
        except IOError:
            createdir(filename)
            fd = open(filename, 'w')
        fd.write('f            ZDip.real            ZDip.imaginary \n')
        fd.write('[Hz]           [Ohm/m]            [Ohm/m] \n')
        for i in range(len(self.f)):
            fd.write('%e      %e      %e \n' %(self.f[i], ZDip[i].real, ZDip[i].imag))
        fd.close()

    def save_ZTrans(self,filename, ZDipX, ZDipY, ZQuadX, ZQuadY, label):
        print ('Saving %s data in %s' %(label, filename))
        try: 
            fd = open(filename, 'w')
        except IOError:
            createdir(filename)
            fd = open(filename, 'w')
            
        fd.write('f            ZDipX.real        ZDipX.imaginary      ZDipY.real        ZDipY.imaginary       ZQuadX.real        ZQuadX.imaginary      ZQuadY.real        ZQuadY.imaginary\n')
        fd.write('[Hz]           [Ohm/m]            [Ohm/m]           [Ohm/m]            [Ohm/m]           [Ohm/m]            [Ohm/m]           [Ohm/m]            [Ohm/m] \n')
        for i in range(len(self.f)):
            fd.write('%e           %e           %e           %e           %e           %e           %e           %e           %e\n' %(self.f[i], ZDipX[i].real, ZDipX[i].imag, ZDipY[i].real, ZDipY[i].imag, ZQuadX[i].real, ZQuadX[i].imag, ZQuadY[i].real, ZQuadY[i].imag))
        fd.close()

    #
    #  Input functions
    #
    def read_pytlwall_output(self, filename):
        print ('Loading %s' %(filename))
        fd = open(filename, 'r')
        arr = fd.readline().split()
        n_col = len(arr)
        unit_meas = fd.readline().split()
        data = np.fromfile(fd, sep='\n')
        fd.close()
        n_row = len(data)/ n_col
        data = data.reshape(n_row, n_col)

        data_dict = {}
        for i  in range(n_col):
            print ('Reading %s in %s'%(arr[i], unit_meas[i]))
            data_dict[arr[i]] =[ data[:,i], unit_meas[i]]
    
        return data_dict

    #
    # plot functions
    #

    def plot_ZLong(self, filename):
        self.plot_Z_vs_f_simple(self.f, self.ZLong,  'L', self.img_ZLong_title, filename,  xscale= self.img_ZLong_freq_scale, yscale=self.img_ZLong_scale)
        
    def plot_ZDip(self, filename):
        self.plot_Z_vs_f_simple(self.f, self.ZDip,  'D', self.img_ZDip_title, filename,  xscale= self.img_ZDip_freq_scale, yscale=self.img_ZDip_scale)

    def plot_ZDipX(self, filename):
        self.plot_Z_vs_f_simple(self.f, self.ZDipX,  'D', self.img_ZDipX_title, filename,  xscale= self.img_ZDipX_freq_scale, yscale=self.img_ZDipX_scale)

    def plot_ZDipY(self, filename):
        self.plot_Z_vs_f_simple(self.f, self.ZDipY,  'D', self.img_ZDipY_title, filename,  xscale= self.img_ZDipY_freq_scale, yscale=self.img_ZDipY_scale)

    def plot_ZQuadX(self, filename):
        self.plot_Z_vs_f_simple(self.f, self.ZQuadX,  'D', self.img_ZQuadX_title, filename,  xscale= self.img_ZQuadX_freq_scale, yscale=self.img_ZQuadX_scale)

    def plot_ZQuadY(self, filename):
        self.plot_Z_vs_f_simple(self.f, self.ZQuadY,  'D', self.img_ZQuadY_title, filename,  xscale= self.img_ZQuadY_freq_scale, yscale=self.img_ZQuadY_scale)


    def plot_ZLongISC(self, filename):
        self.plot_Z_vs_f_simple_single(self.f, self.ZLongISC.imag,  'D', self.img_ZLongISC_title, filename,  xscale= self.img_ZLongISC_freq_scale, yscale=self.img_ZLongISC_scale)

    def plot_ZLongDSC(self, filename):
        self.plot_Z_vs_f_simple_single(self.f, self.ZLongDSC.imag,  'D', self.img_ZLongDSC_title, filename,  xscale= self.img_ZLongDSC_freq_scale, yscale=self.img_ZLongDSC_scale)


    def plot_ZDipISC(self, filename):
        self.plot_Z_vs_f_simple_single(self.f, self.ZDipISC.imag,  'D', self.img_ZDipISC_title, filename,  xscale= self.img_ZDipISC_freq_scale, yscale=self.img_ZDipISC_scale)

    def plot_ZDipDSC(self, filename):
        self.plot_Z_vs_f_simple_single(self.f, self.ZDipDSC.imag,  'D', self.img_ZDipDSC_title, filename,  xscale= self.img_ZDipDSC_freq_scale, yscale=self.img_ZDipDSC_scale)
        
    def plot_Z_vs_f_simple(self, f, Z,  imped_type, title, savename,  xscale= 'lin', yscale='lin',):
        import matplotlib.pyplot as pl
        
        mask_real = np.logical_not(np.isnan(Z.real))
        fRe = f[mask_real]
        ZRe = Z[mask_real].real
        mask_imag = np.logical_not(np.isnan(Z.imag))
        fIm = f[mask_imag]
        ZIm = Z[mask_imag].imag
        
        if (len(fRe) == 0. and len(fIm) == 0.):
            return
        print (savename)
        if imped_type == 'L':
            Z_unit = '$\Omega$'
        else :
            Z_unit = '$\Omega$/m'
            
            
        fig = pl.figure()
        ax = pl.subplot(1,1,1)
        if title != None:
            pl.title(title, fontsize=24)
        ax.plot(fRe, ZRe, linewidth= 5, label='real '  )
        ax.plot(fIm, ZIm, linewidth= 5, label='imaginary '  )
        ax.set_ylabel('Z ['+ Z_unit + ']', fontsize=20)
        ax.set_xlabel('f [Hz]', fontsize=20)
        pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if xscale == 'log':
            pl.xscale('log')
        if xscale == 'symlog':
            pl.xscale('symlog')    
        if yscale == 'symlog':
            pl.yscale('symlog')
        if yscale == 'log':
            pl.yscale('log')        
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(18)         
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(18)     

        ax.grid(True)
        
        ax.legend(loc = 'best', fontsize=20)
        pl.tight_layout()
        if savename !=None:
            try: 
                pl.savefig(savename)
            except IOError:
                createdir(savename)

                pl.savefig(savename)
        return fig

    def plot_Z_vs_f_simple_single(self, f, Z,  imped_type, title, savename,  xscale= 'lin', yscale='lin'):
        import matplotlib.pyplot as pl
        
        mask = np.logical_not(np.isnan(Z))
        f = f[mask]
        Z = Z[mask]

        if (len(f) == 0. ):
            return
            
        if imped_type == 'L':
            Z_unit = '$\Omega$'
        else :
            Z_unit = '$\Omega$/m'
            
            
        fig = pl.figure()
        ax = pl.subplot(1,1,1)
        if title != None:
            pl.title(title, fontsize=24)
        ax.plot(f, Z, linewidth= 5 )
        ax.set_ylabel('Z ['+ Z_unit + ']', fontsize=20)
        ax.set_xlabel('f [Hz]', fontsize=20)
        pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if xscale == 'log':
            pl.xscale('log')
        if xscale == 'symlog':
            pl.xscale('symlog')    
        if yscale == 'symlog':
            pl.yscale('symlog')
        if yscale == 'log':
            pl.yscale('log')        
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(18)         
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(18)     

        ax.grid(True)
        pl.tight_layout()
        if savename !=None:
            try: 
                pl.savefig(savename)
            except IOError:
                createdir(savename)

                pl.savefig(savename)
        return fig

    def plot_Z_vs_f_simple_single_compare(self, list_f, list_Z, list_label,  imped_type, title, savename,  xscale= 'lin', yscale='lin'):
        import matplotlib.pyplot as pl
        

        if imped_type == 'L':
            Z_unit = '$\Omega$'
        else :
            Z_unit = '$\Omega$/m'
            
            
        fig = pl.figure()
        ax = pl.subplot(1,1,1)
        if title != None:
            pl.title(title, fontsize=24)
            
        for i in range(len(list_label)):
            ax.plot(list_f[i], list_Z[i], linewidth= 6-2* i , label = list_label[i])
        ax.set_ylabel('Z ['+ Z_unit + ']', fontsize=20)
        ax.set_xlabel('f [Hz]', fontsize=20)
        pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if xscale == 'log':
            pl.xscale('log')
        if xscale == 'symlog':
            pl.xscale('symlog')    
        if yscale == 'symlog':
            pl.yscale('symlog')
        if yscale == 'log':
            pl.yscale('log')        
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(18)         
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(18)     

        ax.grid(True)
        ax.legend(loc = 'best', fontsize=20)
        pl.tight_layout()
        if savename !=None:
            try: 
                pl.savefig(savename)
            except IOError:
                createdir(savename)

                pl.savefig(savename)
        pl.show()
        return fig
