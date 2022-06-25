# importing used modules
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
# functions searching for location of maximal values in two-dimensional array


def max_location_phi(array_2d):                         # phi is fuel–air equivalence ratio
    result = np.where(array_2d == np.amax(array_2d))
    return result[0][0]


def max_location_x(array_2d):                           # x is methane's share in fuel (in loop named also 'fuels')
    result = np.where(array_2d == np.amax(array_2d))
    return result[1][0]


# used mechanism in Cantera and defining of arrays (analysed: methane shares, fuel-air ratios, max. temp.,
# max. pres. and max. pres. rate in each case
gas = ct.Solution('gri30.yaml')
points_fuels = 21
points_phis = 21
fuels = np.linspace(0., 1., points_fuels)
phis = np.linspace(0., 10., points_phis)
max_temperatures = np.zeros((points_fuels, points_phis))
max_pressures = np.zeros((points_fuels, points_phis))
max_pres_rates = np.zeros((points_fuels, points_phis))

# list of analysed cases
switch_case = [[890., ct.one_atm], [950., ct.one_atm], [1000., ct.one_atm],  [1050., ct.one_atm],  [1100., ct.one_atm],
               [950., ct.one_atm], [950., 1.5*ct.one_atm], [950., 2.*ct.one_atm],  [950., 3*ct.one_atm]]
for k in range(len(switch_case)):
    temperature_0 = switch_case[k][0]
    pressure_0 = switch_case[k][1]
    for i in range(points_phis):
        for j in range(points_fuels):
            # defining shares of methane and ethane in fuel and needed amount of oxygen for stoichiometric mixture
            methane_share = fuels[j]                    # methane's share in fuel
            ethane_share = 1.-methane_share  # ethane's share in fuel
            needed_O2 = (methane_share+2*ethane_share)+(2*methane_share+3*ethane_share)/2

            # mixture properties
            phi = phis[i]                               # fuel–air equivalence ratio
            if phi != 0:
                gas.X = {'CH4': methane_share, 'C2H6': ethane_share, 'O2': needed_O2/phi, 'N2': needed_O2*3.76/phi}
            else:
                gas.X = {'CH4': 0, 'C2H6': 0, 'O2': needed_O2, 'N2': needed_O2 * 3.76}
            gas.TP = temperature_0, pressure_0

            # simulation of autoignition process - preparation
            max_temp = gas.T
            max_pres = gas.P
            max_d_pres = 0.0
            previous_pres = ct.one_atm
            d_time = 0.001
            r = ct.IdealGasReactor(gas)
            sim = ct.ReactorNet([r])
            time = 0.0

            # simulation of autoignition process - calculations and saving data
            # (max. temp., max. pres. and max. pres. rate in each case)
            for t in range(1000):
                time += d_time
                sim.advance(time)
                if gas.T > max_temp:
                    max_temp = gas.T
                if gas.P > max_pres:
                    max_pres = gas.P
                if (gas.P-previous_pres)/d_time > max_d_pres:
                    max_d_pres = (gas.P-previous_pres)/d_time
                previous_pres = gas.P

            max_temperatures[i][j] = max_temp
            max_pressures[i][j] = max_pres
            max_pres_rates[i][j] = max_d_pres

    # plots of max. temp., max. pres. and max. pres. rate
    plot_dif = 11
    cm = 1/2.54                                         # centimeters in inches

    fig1, ax = plt.subplots(constrained_layout=True, figsize=(30*cm, 20*cm))
    FIG1 = ax.contourf(fuels, phis, max_temperatures, plot_dif-1, cmap=plt.cm.bone, origin='lower')
    CS1 = ax.contour(FIG1, levels=FIG1.levels[::1], colors='w', origin='lower')
    ax.set_title('Maximal temperature during explosion (T=' + str(int(temperature_0)) + 'K, p='
                 + str(pressure_0/ct.one_atm) + 'atm)')
    ax.set_xlabel('Methane\'s share in fuel [1]')
    ax.set_ylabel('Fuel–air equivalence ratio \u03A6 [1]')
    cbar = fig1.colorbar(FIG1)
    cbar.ax.set_ylabel('Temperature [K]')
    ax.clabel(CS1, fmt='%.0f', levels=FIG1.levels[::2], colors='w', fontsize=10)
    cbar.add_lines(CS1)
    name_temp1 = 'Figures/plot_temp_T' + str(int(temperature_0)) + '_p' + str(int(pressure_0)) + '.png'
    plt.savefig(name_temp1, dpi=1000)

    fig2, ax = plt.subplots(constrained_layout=True, figsize=(30*cm, 20*cm))
    FIG2 = ax.contourf(fuels, phis, max_pressures/1e5, plot_dif-1, cmap=plt.cm.bone, origin='lower')
    CS2 = ax.contour(FIG2, levels=FIG2.levels[::1], colors='w', origin='lower')
    ax.set_title('Maximal pressure during explosion (T=' + str(int(temperature_0)) + 'K, p='
                 + str(pressure_0/ct.one_atm) + 'atm)')
    ax.set_xlabel('Methane\'s share in fuel [1]')
    ax.set_ylabel('Fuel–air equivalence ratio \u03A6 [1]')
    cbar = fig2.colorbar(FIG2)
    cbar.ax.set_ylabel('Pressure [bar]')
    ax.clabel(CS2, fmt='%.1f', levels=FIG2.levels[::2], colors='w', fontsize=10)
    cbar.add_lines(CS2)
    name_temp2 = 'Figures/plot_pres_T' + str(int(temperature_0)) + '_p' + str(int(pressure_0)) + '.png'
    plt.savefig(name_temp2, dpi=1000)

    fig3, ax = plt.subplots(constrained_layout=True, figsize=(30*cm, 20*cm))
    FIG3 = ax.contourf(fuels, phis, max_pres_rates/1e5, plot_dif-1, cmap=plt.cm.bone, origin='lower')
    CS3 = ax.contour(FIG3, levels=FIG3.levels[::1], colors='w', origin='lower')
    ax.set_title('Maximal pressure rate during explosion (T=' + str(int(temperature_0)) + 'K, p='
                 + str(pressure_0/ct.one_atm) + 'atm)')
    ax.set_xlabel('Methane\'s share in fuel [1]')
    ax.set_ylabel('Fuel–air equivalence ratio \u03A6 [1]')
    cbar = fig3.colorbar(FIG3)
    cbar.ax.set_ylabel('Pressure rate [bar/s]')
    ax.clabel(CS3, fmt='%.1f', levels=FIG3.levels[::2], colors='w', fontsize=10)
    cbar.add_lines(CS3)
    name_temp3 = 'Figures/plot_d_pres_T' + str(int(temperature_0)) + '_p' + str(int(pressure_0)) + '.png'
    plt.savefig(name_temp3, dpi=1000)

    # additional data - specific phis and fuel compositions with the highest max. temp., max. pres. and max. pres. rate
    # in each case
    print('Starting parameters:\tT_0=' + str(int(temperature_0)) + 'K, p_0=' + str(pressure_0/ct.one_atm) + 'atm')
    print('Maximal parameters:\t\tT_max=' + str(int(np.amax(max_temperatures))) + 'K (for (\u03A6,x_CH4)=('
          + format(phis[max_location_phi(max_temperatures)], ".1f") + ','
          + format(fuels[max_location_x(max_temperatures)], ".2f")
          + ')), p=' + str(int(np.amax(max_pressures)/1e5)) + 'bar (for (\u03A6,x_CH4)=('
          + format(phis[max_location_phi(max_pressures)], ".1f") + ','
          + format(fuels[max_location_x(max_pressures)], ".2f")
          + ')), dp=' + str(int(np.amax(max_pres_rates)/1e5)) + 'bar/s (for (\u03A6,x_CH4)=('
          + format(phis[max_location_phi(max_pres_rates)], ".1f") + ','
          + format(fuels[max_location_x(max_pres_rates)], ".2f")
          + '))\n')
