
'''
@authors: Tatiana Rijoff,
          Carlo Zannini
@date:    01/03/2013
@copyright CERN
'''
import configparser
import tlwall


class CfgIo(object):
    def __init__(self, cfg_file=None):
        if cfg_file is not None:
            self.config = self.read_cfg(cfg_file)
        return

    def read_chamber(self, cfg_file=None):
        if cfg_file is not None:
            config = self.read_cfg(cfg_file)
        else:
            config = self.config
        pipe_len_m = config.getfloat('base_info', 'pipe_len_m')
        # pipe radius is equal to vertical dimension so in the cfg file or
        # it is defined the radius or the vertical dimension
        if config.has_option('base_info', 'pipe_radius_m'):
            pipe_rad_m = config.getfloat('base_info', 'pipe_radius_m')
            pipe_hor_m = pipe_rad_m
            pipe_ver_m = pipe_rad_m
        else:
            pipe_ver_m = config.getfloat('base_info', 'pipe_ver_m')
            pipe_rad_m = pipe_ver_m
        if config.has_option('base_info', 'pipe_hor_m'):
            pipe_hor_m = config.getfloat('base_info', 'pipe_hor_m')

        chamber_shape = config.get('base_info', 'chamber_shape')
        if config.has_option('base_info', 'component_name'):
            component_name = config.get('base_info', 'component_name')
        else:
            component_name = 'chamber'
        betax = config.getfloat('base_info', 'betax')
        betay = config.getfloat('base_info', 'betay')
        nbr_layers = config.getint('layers_info', 'nbr_layers')
        layers = []

        for i in range(nbr_layers):
            layer_type = config.get('layer' + str(i), 'type')
            thick_m = config.getfloat('layer' + str(i), 'thick_m')
            if layer_type == 'CW':
                muinf_Hz = config.getfloat('layer' + str(i), 'muinf_Hz')
                epsr = config.getfloat('layer' + str(i), 'epsr')
                sigmaDC = config.getfloat('layer' + str(i), 'sigmaDC')
                k_Hz = config.getfloat('layer' + str(i), 'k_Hz')
                tau = config.getfloat('layer' + str(i), 'tau')
                RQ = config.getfloat('layer' + str(i), 'RQ')
                layers.append(tlwall.Layer(layer_type=layer_type,
                                           thick_m=thick_m,
                                           muinf_Hz=muinf_Hz,
                                           epsr=epsr,
                                           sigmaDC=sigmaDC,
                                           k_Hz=k_Hz,
                                           tau=tau,
                                           RQ=RQ,
                                           boundary=False))
            else:
                layers.append(tlwall.Layer(layer_type=layer_type,
                                           thick_m=thick_m,
                                           boundary=False))

        layer_type = config.get('boundary', 'type')
        if layer_type == 'CW':
            muinf_Hz = config.getfloat('boundary', 'muinf_Hz')
            epsr = config.getfloat('boundary', 'epsr')
            sigmaDC = config.getfloat('boundary', 'sigmaDC')
            k_Hz = config.getfloat('boundary', 'k_Hz')
            tau = config.getfloat('boundary', 'tau')
            RQ = config.getfloat('boundary', 'RQ')
        layers.append(tlwall.Layer(layer_type=layer_type,
                                   thick_m=thick_m,
                                   muinf_Hz=muinf_Hz,
                                   epsr=epsr,
                                   sigmaDC=sigmaDC,
                                   k_Hz=k_Hz,
                                   tau=tau,
                                   RQ=RQ,
                                   boundary=True))
        chamber = tlwall.Chamber(pipe_len_m=pipe_len_m,
                                 pipe_rad_m=pipe_rad_m,
                                 pipe_hor_m=pipe_hor_m,
                                 pipe_ver_m=pipe_ver_m,
                                 chamber_shape=chamber_shape,
                                 betax=betax,
                                 betay=betay,
                                 layers=layers,
                                 component_name=component_name)
        return chamber

    def read_freq(self, cfg_file=None):
        if cfg_file is not None:
            config = self.read_cfg(cfg_file)
        else:
            config = self.config
        fmin = config.getfloat('frequency_info', 'fmin')
        fmax = config.getfloat('frequency_info', 'fmax')
        fstep = config.getfloat('frequency_info', 'fstep')
        freq = tlwall.Frequencies(fmin=fmin, fmax=fmax, fstep=fstep)
        return freq

    def read_beam(self, cfg_file=None):
        if cfg_file is not None:
            config = self.read_cfg(cfg_file)
        else:
            config = self.config
        test_beam_shift = config.getfloat('beam_info', 'test_beam_shift')
        if config.has_option('beam_info', 'betarel'):
            betarel = config.getfloat('beam_info', 'betarel')
            if config.has_option('beam_info', 'mass_MeV_c2'):
                mass_MeV_c2 = config.getfloat('beam_info', 'mass_MeV_c2')
                beam = tlwall.Beam(betarel=betarel,
                                   test_beam_shift=test_beam_shift,
                                   mass_MeV_c2=mass_MeV_c2)
            else:
                beam = tlwall.Beam(betarel=betarel,
                                   test_beam_shift=test_beam_shift)
            return beam
        if config.has_option('beam_info', 'gammarel'):
            gammarel = config.getfloat('beam_info', 'gammarel')
            if config.has_option('beam_info', 'mass_MeV_c2'):
                mass_MeV_c2 = config.getfloat('beam_info', 'mass_MeV_c2')
                beam = tlwall.Beam(gammarel=gammarel,
                                   test_beam_shift=test_beam_shift,
                                   mass_MeV_c2=mass_MeV_c2)
            else:
                beam = tlwall.Beam(gammarel=gammarel,
                                   test_beam_shift=test_beam_shift)
            return beam
        if config.has_option('beam_info', 'Ekin_MeV'):
            Ekin_MeV = config.getfloat('beam_info', 'Ekin_MeV')
            if config.has_option('beam_info', 'mass_MeV_c2'):
                mass_MeV_c2 = config.getfloat('beam_info', 'mass_MeV_c2')
                beam = tlwall.Beam(Ekin_MeV=Ekin_MeV,
                                   test_beam_shift=test_beam_shift,
                                   mass_MeV_c2=mass_MeV_c2)
            else:
                beam = tlwall.Beam(Ekin_MeV=Ekin_MeV,
                                   test_beam_shift=test_beam_shift)
            return beam
        if config.has_option('beam_info', 'p_MeV_c'):
            p_MeV_c = config.getfloat('beam_info', 'p_MeV_c')
            if config.has_option('beam_info', 'mass_MeV_c2'):
                mass_MeV_c2 = config.getfloat('beam_info', 'mass_MeV_c2')
                beam = tlwall.Beam(p_MeV_c=p_MeV_c,
                                   test_beam_shift=test_beam_shift,
                                   mass_MeV_c2=mass_MeV_c2)
            else:
                beam = tlwall.Beam(p_MeV_c=p_MeV_c,
                                   test_beam_shift=test_beam_shift)
            return beam

    def read_tlwall(self, cfg_file=None):
        if cfg_file is not None:
            config = self.read_cfg(cfg_file)
        else:
            config = self.config
        chamber = self.read_chamber()
        beam = self.read_beam()
        freq = self.read_freq()
        mywall = tlwall.TlWall(chamber, beam, freq)
        return mywall

    def read_cfg(self, cfg_file):
        config = configparser.ConfigParser()
        try:
            config.read(cfg_file)
        except IOError:
            print('The file %s does not exist, try again!' % cfg_file)

        return config
