import numpy as np
import sys

import pytlwall as tw

tw.welcome_message()


curr_el = tw.Acc_El( )
beam = tw.Beam()
freq = tw.Freq()

choice = ''
if len(sys.argv) == 2:
    filename = sys.argv[1]
    w = tw.pyTlWall(curr_el, beam, freq, cfg_file = filename)
        
else:
    while choice.upper() != 'X':
        #~ if curr_el.component_name != None:
            #~ print 'Accelerator component  %s:  ' %curr_el.component_name
            #~ print 'len:%f [m] radius:%f [m] shape: %s beta_x [m]: %f beta_y: %f [m]' %( curr_el.pipe_len, curr_el.pipe_radius, curr_el.form_factor, curr_el.betax, curr_el.betay)
            #~ print 'nbr layers :%d  boundary: %s' %(curr_el.nbr_layers,  curr_el.layers[curr_el.nbr_layers-1]['type'])
            #~ print curr_el.layers
        #~ if beam.particle_name != None:
            #~ print 'Beam of %s ' %beam.particle_name
            #~ print 'Ekin = %e betarel = %e gammarel = %e test beam separation = %e' %( beam.Ekin, beam.betarel, beam.gammarel, beam.test_beam_shift )  
        if curr_el.component_name != None and beam.particle_name != None and len(freq.freq) > 0 :
            try:
                w.exist()
            except NameError:
                w = tw.pyTlWall(curr_el, beam, freq)
        print('=====================================================================')
        print('What do you want to do?')
        print('1 Define accelerator component ')
        print('2 Define beam characteristics ')
        print('3 Define frequency range ')
        print('4 Read all the info from a config')
        try:
            w.exist()
            print('5 Calculate impedance')
        except NameError:
            pass
        try:
            w.exist()
            print('6 Plot impedance')
        except NameError:
            pass
        print('X Exit')
        choice = raw_input('your choice: ')
        if choice == '1':
            choice = raw_input('Press ENTER if you want insert interactively all the informations or digit the name of the configurator file to use     ')
            if choice == '':
                curr_el.ask_acc_el_info()
            else:
                curr_el.read_cfg(choice)
                        
        elif choice == '2':
            beam.ask_beam_info()
            
        elif choice == '3':
            choice = raw_input('Press ENTER if you want insert interactively the frequency range or digit the name of the frequency file to use     ')
            if choice == '':
                freq.ask_freq_info()
            else:
                freq.read_freq_from_file(choice)
        elif choice == '4':
            choice = raw_input('Digit the filename     ')
            w = tw.pyTlWall(curr_el, beam, freq, choice)
        elif choice == '5':
            w.ask_impedance_calc()    
        elif choice == '6':
            w.ask_impedance_plot()    
        elif choice.upper() == 'X':
            try:
                w.exist()
            except NameError:
                print('Goodbye')
                exit(1)
            choice = raw_input('If you want to generate a configurator file digit here the name     ')
            if choice != '':
                try:
                    w.generate_cfg(choice)
                except NameError:
                    print("You cannot save data if you haven't defined element component, beam characteristic and frequency ")
            print('Goodbye')
            exit(1)
        else:
            print('I cannot understand %s, please try again' %choice)
        #w = tw.pyTlWall(curr_el, beam, freq)



