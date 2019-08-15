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
    parser = argparse.ArgumentParser(description='Analyze outlers from vapor--liquid coexistence data')
    parser.add_argument('input', type=str, help="path to input file")
    parser.add_argument('weight', type=float, help='molecular weight in g/mol')
    parser.add_argument('critical_temperature', type=float, help='critical temperature of the molecule in K')
    parser.add_argument('-o', '--output', type=str, default='outliers.csv', help='output file')
    parser.add_argument('--skip', type=int, default=4, help='rows to skip in input file')
    return parser.parse_args()

class VLEData:

    # static method to get expected compressibility factor
    get_z_exd = lambda t: 0.2732 * np.log(1.070 - t ** 2.921) + 1.007

    # t_r bounds and ranges to dectect outliers 
    tr_bounds = [0.70, 0.85]
    z_ranges = {'low': [-0.02, 0.05], 'mid': [-0.06, 0.05], 'high': [-0.13, 0.12]}
    p_ranges = {'low': [-0.06, 0.03], 'mid': [-0.06, 0.05], 'high': [-0.07, 0.07]}

    def __init__(self, mol_weight, critical_temp, data):
        self.mol_weight = mol_weight
        self.critical_temp = critical_temp
        # read and convert temperature to K, vapor density 
        # to kg/m^3, pressure to KPa, and reduced temperature
        self.temps = data[:, 0]
        self.rho_v = data[:, 1] * 1000
        self.p_sim = data[:, 5]
        self.t_r = self.temps / self.critical_temp
        # calculate simulated and expected compressibility factor
        self.z_sim = self.p_sim * self.mol_weight / (self.rho_v * R * self.temps)
        self.z_exd = VLEData.get_z_exd(self.t_r)
        # calculated log expected pressure
        self.log_p_cc = self.__get_log_p_cc()
        self.ln_zdiff = None
        self.ln_pdiff = None

    def __get_log_p_cc(self):
        inv_temp = 1 / self.temps
        log_press = np.log(self.p_sim)
        log_pcc_low = linear_interpolate(
                            inv_temp[:-2], 
                            (inv_temp[1:-1], log_press[1:-1]),
                            (inv_temp[2:], log_press[2:]))
        log_pcc_high = linear_interpolate(
                            inv_temp[-2:], 
                            (inv_temp[-3:-1], log_press[-3:-1]),
                            (inv_temp[-4:-2], log_press[-4:-2]))
        return np.concatenate([log_pcc_low, log_pcc_high])

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
            data_col = [key if i in outliers[key] else x for i, x in enumerate(data_col)]
            prefixsum += len(sliced)
        return outliers, data_col

    def find_p_outliers(self):
        outliers = {}
        indices = self.__split_indices()
        prefixsum = 0
        data_col = [''] * len(self.temps)
        self.ln_pdiff = np.log(self.p_sim) - self.log_p_cc
        for key in ['low', 'mid', 'high']:
            sliced = self.ln_zdiff[indices[key]]
            prange = VLEData.p_ranges[key]
            outliers[key] = prefixsum + np.where(np.logical_or(sliced < prange[0], 
                                                            sliced > prange[1]))[0]
            data_col = [key if i in outliers[key] else x for i, x in enumerate(data_col)]
            prefixsum += len(sliced)
        return outliers, data_col

    def write_data(self, writer):
        _, z_outliers = self.find_z_outliers()
        _, p_outliers = self.find_p_outliers()
        columns = ['Temperature (K)', 'Tr', 'Vapor density (kg/m^3)',
                   'Zsim', 'Zexd', 'ln(Zsim/Zexd)',
                   'Psim', 'Pcc', 'ln(Psim/Pcc)',
                   'Z outliers', 'P outliers'
                  ]
        fields = [self.temps, self.t_r, self.rho_v,
                  self.z_sim, self.z_exd, self.ln_zdiff,
                  self.p_sim, np.exp(self.log_p_cc), self.ln_pdiff,
                  z_outliers, p_outliers
                 ]
        writer.writerow(columns)
        for i in range(len(self.temps)):
            writer.writerow([x[i] for x in fields])

    def write_analysis(self, field, writer):
        if field is self.ln_zdiff:
            field_range = VLEData.z_ranges
            name = 'ln(Zsim/Zexd)'
        elif field is self.ln_pdiff:
            field_range = VLEData.p_ranges
            name = 'ln(Psim/Pcc)'
        else:
            raise ValueError('Invalid field to analyze!')
        writer.writerow([name,'Tr <= 0.7', '','','0.7 < Tr <= 0.85', '','','Tr > 0.85','',''])
        mean_row = []
        std_row = []
        range_row = []
        indices = self.__split_indices()
        for key in ['low', 'mid', 'high']:
            mean_row.extend(['', 'Mean', np.mean(field[indices[key]])])
            std_row.extend(['', 'Standard Deviation', np.std(field[indices[key]])])
            range_row.extend(['Lower/upper bounds'] + field_range[key])
        writer.writerow(mean_row)
        writer.writerow(std_row)
        writer.writerow(range_row)

    def report(self, f): 
        csvwriter = csv.writer(f, delimiter=',', quotechar='"', 
                               lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC)
        csvwriter.writerow(['Molecular weight (g/mol)', self.mol_weight])
        csvwriter.writerow(['Gas Constant [J/(K mol)]', R])
        csvwriter.writerow('')
        self.write_data(csvwriter)
        csvwriter.writerow('')
        self.write_analysis(self.ln_zdiff, csvwriter)
        self.write_analysis(self.ln_pdiff, csvwriter)

if __name__ == '__main__':
    args = parse_args()
    # columns: T, vapor density, error, liquid density, error, pressure, error
    data = np.loadtxt(args.input, skiprows=args.skip)
    vle = VLEData(args.weight, args.critical_temperature, data)
    with open(args.output, 'w+') as f:
        vle.report(f)
