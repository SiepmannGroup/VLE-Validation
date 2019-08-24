import argparse, csv
import numpy as np

# gas constant in J/(K mol)
R = 8.3144598

# interpolate the y value given x using
# the straight line between tuples p1(x1, y1) and p2(x2, y2)
def linear_interpolate(x, p1, p2):
    x1, y1 = p1
    x2, y2 = p2 
    return (y2 - y1) / (x2 - x1) * (x - x1) + y1

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze outliers from vapor--liquid coexistence data')
    parser.add_argument('input', type=str, help="path to input file")
    parser.add_argument('weight', type=float, help='molecular weight in g/mol')
    parser.add_argument('critical_temperature', type=float, help='critical temperature of the molecule in K')
    parser.add_argument('-o', '--output', type=str, default='outliers.csv', help='output file')
    return parser.parse_args()

class VLEData:

    # static method to get expected compressibility factor
    get_z_exd = lambda t: 0.2732 * np.log(1.070 - t ** 2.921) + 1.007

    # t_r bounds and ranges to dectect outliers 
    tr_bounds = [0.70, 0.85]
    z_ranges = {'low': [-0.02, 0.05], 'mid': [-0.06, 0.05], 'high': [-0.13, 0.12]}
    p_ranges = {'cc_low': {'low': [-0.06, 0.03], 'mid': [-0.06, 0.05], 'high': [-0.07, 0.07]},
                'cc_mid': {'low': [-0.02, 0.05], 'mid': [-0.03, 0.04], 'high': [-0.03, 0.04]},
                'cc_high': {'low': [-0.06, 0.01], 'mid': [-0.06, 0.04], 'high': [-0.05, 0.05]}
               }

    def __init__(self, mol_weight, critical_temp, data, columns):
        self.mol_weight = mol_weight
        self.critical_temp = critical_temp
        # read and convert temperature to K, vapor density 
        # to kg/m^3, pressure to KPa, and reduced temperature
        self.temps = data[:, columns['temp']]
        self.rho_v = data[:, columns['rho']] * 1000
        self.p_sim = data[:, columns['pres']]
        self.t_r = self.temps / self.critical_temp
        # calculate simulated and expected compressibility factor
        self.z_avg = self.p_sim * self.mol_weight / (self.rho_v * R * self.temps)
        if 'z' in columns:
            self.z_sim = data[:, columns['z']]
        else:
            print("Warning: compressibility factor was not found in input file. \
                   Compressibility will be calculated using average pressure and density"
                    )
        self.z_exd = VLEData.get_z_exd(self.t_r)
        # calculated log expected pressure
        self.log_p_cc = self.__get_log_p_cc()
        self.ln_zdiff = None
        self.ln_pdiff = {}

    def __get_log_p_cc(self):
        inv_temp = 1 / self.temps
        log_press = np.log(self.p_sim)
        log_pcc_low = linear_interpolate(
                            inv_temp[:-2], 
                            (inv_temp[1:-1], log_press[1:-1]),
                            (inv_temp[2:], log_press[2:]))
        log_pcc_mid = linear_interpolate(
                            inv_temp[1:-1], 
                            (inv_temp[:-2], log_press[:-2]),
                            (inv_temp[2:], log_press[2:]))
        log_pcc_high = linear_interpolate(
                            inv_temp[2:], 
                            (inv_temp[1:-1], log_press[1:-1]),
                            (inv_temp[:-2], log_press[:-2]))
        log_p_cc = np.tile(log_press, [3, 1])
        log_p_cc[0, :-2] = log_pcc_low
        log_p_cc[1, 1:-1] = log_pcc_mid
        log_p_cc[2, 2:] = log_pcc_high
        return {'cc_low': log_p_cc[0], 'cc_mid': log_p_cc[1], 'cc_high': log_p_cc[2]}

    def __split_indices(self):
        ind_low = np.where(self.t_r <= VLEData.tr_bounds[0])
        ind_mid = np.where(np.logical_and(self.t_r > VLEData.tr_bounds[0],
                                          self.t_r <= VLEData.tr_bounds[1]))
        ind_high = np.where(self.t_r > VLEData.tr_bounds[1])
        return {'low': ind_low, 'mid': ind_mid, 'high': ind_high}

    def find_z_outliers(self):
        outliers = {}
        data_col = [''] * len(self.temps)
        indices = self.__split_indices()
        prefixsum = 0
        self.ln_zdiff = np.log(self.z_sim / self.z_exd)
        for key in ['low', 'mid', 'high']:
            sliced = self.ln_zdiff[indices[key]]
            zrange = VLEData.z_ranges[key]
            outliers[key] = prefixsum + np.where(np.logical_or(sliced < zrange[0], 
                                                            sliced > zrange[1]))[0]
            data_col = ['X' if i in outliers[key] else x for i, x in enumerate(data_col)]
            prefixsum += len(sliced)
        return outliers, data_col

    def find_p_outliers(self, cc_type):
        outliers = {}
        indices = self.__split_indices()
        prefixsum = 0
        data_col = [''] * len(self.temps)

        self.ln_pdiff[cc_type] = np.log(self.p_sim) - self.log_p_cc[cc_type]
        for key in ['low', 'mid', 'high']:
            sliced = self.ln_pdiff[cc_type][indices[key]]
            prange = VLEData.p_ranges[cc_type][key]
            outliers[key] = prefixsum + np.where(np.logical_or(sliced < prange[0], 
                                                            sliced > prange[1]))[0]
            data_col = ['X' if i in outliers[key] else x for i, x in enumerate(data_col)]
            prefixsum += len(sliced)
        return outliers, data_col

    def write_data(self, writer):
        _, z_outliers = self.find_z_outliers()
        _, p_outliers_low = self.find_p_outliers('cc_low')
        _, p_outliers_mid = self.find_p_outliers('cc_mid')
        _, p_outliers_high = self.find_p_outliers('cc_high')

        columns = ['Temperature (K)', 'Tr', 'Vapor density (kg/m^3)',
                   'Zsim', 'Zexd', 'ln(Zsim/Zexd)',
                   'Psim',
                   'Pcc_LOW', 'ln(Psim/Pcc_LOW)',
                   'Pcc_MID', 'ln(Psim/Pcc_MID)',
                   'Pcc_HIGH', 'ln(Psim/Pcc_HIGH)',
                   'Z outliers', 'P_CC_LOW outliers',
                   'P_CC_MID outliers', 'P_CC_HIGH outliers',
                  ]
        fields = [self.temps, self.t_r, self.rho_v,
                  self.z_sim, self.z_exd, self.ln_zdiff,
                  self.p_sim, 
                  np.exp(self.log_p_cc['cc_low']), self.ln_pdiff['cc_low'],
                  np.exp(self.log_p_cc['cc_mid']), self.ln_pdiff['cc_mid'],
                  np.exp(self.log_p_cc['cc_high']), self.ln_pdiff['cc_high'],
                  z_outliers, p_outliers_low, p_outliers_mid, p_outliers_high
                 ]
        writer.writerow(columns)
        indices = {k: v[0].tolist() for k, v in self.__split_indices().items()}
        ranges_seen = {k: len(indices[k]) == 0 for k in ['low', 'mid', 'high']}
        for i in range(len(self.temps)):
            row = [x[i] for x in fields]
            key = 'low' if i in indices['low'] else ('mid' if i in indices['mid'] else 'high')
            if not ranges_seen[key]:
                ranges_seen[key] = True
                upper_bound_row = ['T_r range: %s' % key,'Upper bound'] + [''] * (len(columns) - 2)
                lower_bound_row = ['', 'Lower bound'] + [''] * (len(columns) - 2)
                upper_bound_row[columns.index('ln(Zsim/Zexd)')] = VLEData.z_ranges[key][1]
                lower_bound_row[columns.index('ln(Zsim/Zexd)')] = VLEData.z_ranges[key][0]
                for cc_type in ['cc_LOW', 'cc_MID', 'cc_HIGH']:
                    upper_bound_row[columns.index('ln(Psim/P%s)' % cc_type)] \
                                = VLEData.p_ranges[cc_type.lower()][key][1]
                    lower_bound_row[columns.index('ln(Psim/P%s)' % cc_type)] \
                                = VLEData.p_ranges[cc_type.lower()][key][0]
                writer.writerow(lower_bound_row)
                writer.writerow(upper_bound_row)
                
            if i < 2:
                row[columns.index('Pcc_HIGH')] = '-'
                row[columns.index('ln(Psim/Pcc_HIGH)')] = '-'
                row[columns.index('P_CC_HIGH outliers')] = '-'
            if i < 1 or i >= len(self.temps) - 1:
                row[columns.index('Pcc_MID')] = '-'
                row[columns.index('ln(Psim/Pcc_MID)')] = '-'
                row[columns.index('P_CC_MID outliers')] = '-'
            if i >= len(self.temps) - 2:
                row[columns.index('Pcc_LOW')] = '-'
                row[columns.index('ln(Psim/Pcc_LOW)')] = '-'
                row[columns.index('P_CC_LOW outliers')] = '-'
            writer.writerow(row)

    def report(self, f): 
        csvwriter = csv.writer(f, delimiter=',', quotechar='"', 
                               lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC)
        csvwriter.writerow(['Molecular weight (g/mol)', self.mol_weight])
        csvwriter.writerow(['Gas Constant [J/(K mol)]', R])
        csvwriter.writerow('')
        self.write_data(csvwriter)
        csvwriter.writerow('')

if __name__ == '__main__':
    args = parse_args()
    # columns: T, vapor density, error, liquid density, error, pressure, error
    data = np.loadtxt(args.input, skiprows=1)
    columns = {}
    with open(args.input, 'r') as f:
        header = f.readline().strip().lower().split()
        if 't' in header:
            columns['temp'] = header.index('t')
        elif 'temperature' in header:
            columns['temp'] = header.index('temperature')
        else:
            raise ValueError("Temperature not found in input file!")

        if 'rhog' in header:
            columns['rho'] = header.index('rhog')
        elif 'density' in header:
            columns['rho'] = header.index('density')
        else:
            raise ValueError("Vapor density not found in input file!")
        
        if 'p' in header:
            columns['pres'] = header.index('p')
        elif 'pressure' in header:
            columns['pres'] = header.index('pressure')
        else:
            raise ValueError("Pressure not found in input file!")
        if 'z' in header:
            columns['z'] = header.index('z')
        elif 'compressibility' in header:
            columns['z'] = header.index('compressibility')


    vle = VLEData(args.weight, args.critical_temperature, data, columns)
    with open(args.output, 'w+') as f:
        vle.report(f)