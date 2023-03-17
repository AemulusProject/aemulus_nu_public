import sys, os
import ctypes
import numpy as np
from classy import Class

def get_growth_expansion_info(params, outbase, sysname):

    sys.path.append('../exec/{}/monofonic/build/_deps/zwindstroom-src/'.format(sysname))
    import zwindstroom.units as units
    import zwindstroom.cosmology as cosmology
    import zwindstroom.backend as backend
    import zwindstroom.fluid as fluid
    import zwindstroom.class_link as class_link
    

    M_nu = [params['nu_mass_ev_ov_3']]
    deg_nu = [3.0] # degeneracies
    N_nu = len(M_nu)

    # Initialise a unit system (default uses Mpc lengths and km/s velocities)
    unit_system, physical_consts = units.init_units()

    # We want to integrate the cosmological tables starting at this scale factor
    a_start = 1e-3
    h = params['h']

    # Set up a cosmological model
    params_zw = {"h": h,
                 "Omega_b": params['OmegaB'],
                 "Omega_c": params['OmegaCDM'],
                 "N_ur": 0.00641000,
                 "N_nu": N_nu,
                 "M_nu": M_nu,
                 "deg_nu": deg_nu,
                 "T_nu_0": 1.95,
                 "T_CMB_0": 2.7255,
                 "w0": params['w0'],
                 "wa": params['wa']}

    model = cosmology.MODEL()
    model.set(params_zw)
    # Let zwindstroom know how neutrinos are treated in the cosmological simulation
    # The options are:
    # + 0 (fully relativistic)
    # + 1 (relativistic Hubble rate, but constant particle masses)
    # + 2 (fully non-relativistic)
    model.set_sim_type(1)

    print("Integrating cosmological tables.")

    # Integrate tables of cosmological background quantities
    model.compute(unit_system, physical_consts, a_start)

    # Tolerance and step parameter for the fluid integration
    tol = 1e-10
    hstart = 1e-10

    # We want to integrate the fluid equations between these two points
    # For rescaling, this would be the starting and final times of the simulation
    zini = params['zini']
    a_start_fl = 1.0 / (zini+1)
    a_final_fl = 1.0

    # Prepare the fluid integrator
    fluid.prepare_fluid_integrator(model, unit_system, physical_consts, model.tables, tol, hstart)

    print("Running CLASS.")

    # Maximum wavenumber for the growth factors (default units = 1 / Mpc)
    k_max = 10.0

    # Run CLASS on this model
    cosmo = class_link.run_class(model, unit_system, a_start_fl, k_max)

    # Extract growth factors, growth rates, and wavenumbers from CLASS
    k = class_link.get_wavenumbers(model, cosmo, unit_system)
    g = class_link.get_growth_rates(model, cosmo, a_start_fl)
    D = class_link.get_growth_factors(model, cosmo, a_start_fl)
    nk = len(k)

    print("Integrating fluid equations.")

    # Compute Newtonian growth factors with the fluid approximation
    D_cdm = np.zeros(nk)
    D_b = np.zeros(nk)
    D_nu = np.zeros((nk, N_nu))

    f_cdm = np.zeros(nk)
    f_b = np.zeros(nk)
    f_nu = np.zeros((nk, N_nu))

    for i in range(nk):
        delta_n = [D["d_ncdm[%d]" % j][i] for j in range(N_nu)]
        gn = [g["d_ncdm[%d]" % j][i] for j in range(N_nu)]
        Dn = [0.0 for j in range(N_nu)]

        growth_factors = fluid.GROWTH_FACTORS()
        growth_factors.k = k[i]
        growth_factors.delta_c = D["d_cdm"][i]
        growth_factors.delta_b = D["d_b"][i]
        growth_factors.delta_n = (ctypes.c_double * N_nu)(*delta_n)
        growth_factors.gc = g["d_cdm"][i]
        growth_factors.gb = g["d_b"][i]
        growth_factors.gn = (ctypes.c_double * N_nu)(*gn)
        growth_factors.Dc = 0.
        growth_factors.Db = 0.
        growth_factors.Dn = (ctypes.c_double * N_nu)(*Dn)

        fluid.integrate_fluid_equations(model, unit_system, physical_consts, model.tables, growth_factors, a_start_fl, a_final_fl)

        D_cdm[i] = growth_factors.Dc
        D_b[i] = growth_factors.Db
        D_nu[i,:] = np.array([growth_factors.Dn[i] for i in range(N_nu)])
    
        f_cdm[i] = growth_factors.gc
        f_b[i] = growth_factors.gb
        f_nu[i,:] = np.array([growth_factors.gn[i] for i in range(N_nu)])    

        # Clean up the integrator
    fluid.free_fluid_integrator()


    print("Done with integrations.")

    params_class = {'h': h,
                    'Omega_b': params['OmegaB'],
                    'Omega_cdm': params['OmegaCDM'],
                    'N_ur': 0.00641,
                    'N_ncdm': 1,
                    'output': 'mPk',
                    'z_pk': '0.0,99',
                    'P_k_max_h/Mpc': 20.,
                    'm_ncdm': params['nu_mass_ev_ov_3'],
                    'deg_ncdm': 3,
                    'T_cmb': 2.7255,
                    'A_s': params['As'],
                    'n_s': params['SpectralIndex']}


    print('Running class again')

    pkclass = Class()
    pkclass.set(params_class)
    pkclass.compute()

    P0 = np.array( [pkclass.pk_lin(ki*h, 0.0) * h**3 for ki in k] )
    ocb = params_class['Omega_b'] + params_class['Omega_cdm']
    fr_b = params_class['Omega_b'] / ocb
    fr_c = params_class['Omega_cdm'] / ocb

    omega_nu = params['OmegaNu']
    fnu = omega_nu / (params_class['Omega_cdm'] + params_class['Omega_b'] + omega_nu)
    Dm = fnu*D_nu[:,0]+(1.-fnu)*(fr_b * D_b + fr_c * D_cdm)

    Pb = D_b*D_b*P0;
    Pc = D_cdm*D_cdm*P0;
    Pn = D_nu[:,0]*D_nu[:,0]*P0;
    Pm = Dm*Dm*P0;

    Tb = np.sqrt(Pb/Pm);
    Tc = np.sqrt(Pc/Pm);
    Tn = np.sqrt(Pn/Pm);
    Tm = np.sqrt(Pm/Pm);
    
    T  = np.zeros((len(k), 7))
    T[:,0] = k
    T[:,1] = Tc 
    T[:,2] = Tb
    T[:,5] = Tn
    T[:,6] = Tm
    
    P_ini = np.zeros((len(k),2))
    P_ini[:,0] = k
    P_ini[:,1] = Pm

    P_final = np.zeros((len(k),2))
    P_final[:,0] = k
    P_final[:,1] = P0
    
    f_cb = (fr_b * f_b + fr_c * f_cdm)
    
    f_b_out = np.zeros((len(k),2))
    f_cb_out = np.zeros((len(k),2))
    f_n_out = np.zeros((len(k),2))

    f_b_out[:,0] = k
    f_cb_out[:,0] = k
    f_n_out[:,0] = k
    
    f_b_out[:,1] = f_b
    f_cb_out[:,1] = f_cb
    f_n_out[:,1] = f_nu[:,0]
    

    print('Creating H(z) table')
    
    a_tab_start = a_start_fl
    a_tab_final = a_final_fl
    na = 1000
    avec = np.exp(np.linspace(np.log(a_tab_start), np.log(a_tab_final), na))
    Hvec = np.zeros_like(avec)
    zvec = 1/avec - 1

    # Extract Hubble rate as a function of time
    for i in range(na):
        Hvec[i] = model.get_H_of_a(avec[i])

    cosmology.free_cosmology_tables(model.tables)
    # Write the output to a text file
    z_max = np.max(zvec)
    z_min = 0.0
    bins = 1000
    z_new = np.logspace(np.log10(1.0+z_min),np.log10(1.0+z_max),bins)
    Hz_new = np.interp(z_new-1,zvec[::-1],Hvec[::-1])/model.h #/1000.0
    wvec = np.ones(bins) * params['w0']
    H_table = np.array([z_new[::-1], wvec[::-1], Hz_new[::-1]]).T

    try:
        os.mkdir(outbase)
    except:
        pass

    fname = "{}/zwind_hubble.txt".format(outbase)
    
    np.savetxt(fname, H_table)
    np.savetxt('{}/P_mm_rescaled_z{}_zwind.txt'.format(outbase, zini), P_ini)
    np.savetxt('{}/P_mm_rescaled_z{}_zwind.txt'.format(outbase, 0.0), P_final)    
    np.savetxt('{}/transfer_rescaled_z{}_zwind.txt'.format(outbase, zini), T)
    np.savetxt('{}/f_b_z{}_zwind.txt'.format(outbase, zini), f_b_out)
    np.savetxt('{}/f_cb_z{}_zwind.txt'.format(outbase, zini), f_cb_out)
    np.savetxt('{}/f_n_z{}_zwind.txt'.format(outbase, zini), f_n_out)

