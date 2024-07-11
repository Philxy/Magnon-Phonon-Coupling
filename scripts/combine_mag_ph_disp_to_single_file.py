




of = open('Parameters/combined.txt', 'w')



out_lines_ev = []
out_lines_k = []


with open('Parameters/20x20x20_ph_disp.txt', 'r') as pf:
    pf.readline()
    for line in pf:

        line = line.strip()
        line_split = line.split(' ')

        data = line_split[:3]
        out_lines_k.append(data)


with open('Outputs/20x20x20/ev.txt', 'r') as mf:
    for line in mf:

        line = line.strip()
        line_split = line.split(',')
        ev = line_split[4:]
        out_lines_ev.append(ev)

print(len(out_lines_ev), len(out_lines_k))
assert len(out_lines_ev) == len(out_lines_k)



for i in range(len(out_lines_ev)):
    out_line = out_lines_k[i] + out_lines_ev[i]
    out_line = ' '.join(out_line) + '\n'
    of.write(out_line)