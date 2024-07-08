




of = open('Parameters/combined.txt', 'w')



out_lines_mag = []
out_lines_ph = []


with open('Parameters/20x20x20_ph_disp.txt', 'r') as pf:
    pf.readline()
    for line in pf:

        line = line.strip()
        line_split = line.split(' ')

        data = line_split[:6]
        out_lines_ph.append(data)


with open('Outputs/mag_disp_20x20x20.txt', 'r') as mf:
    mf.readline()
    for line in mf:

        line = line.strip()
        line_split = line.split(',')
        mag_energy = line_split[3].replace('(', '')
        out_lines_mag.append(mag_energy)


assert len(out_lines_ph) == len(out_lines_mag)


for i in range(len(out_lines_ph)):
    out_line = out_lines_ph[i] + [out_lines_mag[i]]
    out_line = ' '.join(out_line) + '\n'
    of.write(out_line)