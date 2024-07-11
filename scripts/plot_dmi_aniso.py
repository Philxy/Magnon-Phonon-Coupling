import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use('science')

eigenenergies_all = []
eigenenergies_dmi = []

with open("Outputs/dmi_and_aniso/GHz_ev.txt", "r") as file:
    for line in file:
        line = line.strip()
        line = line.split(",")
        energies = [float(e) for e in line[4:]]
        eigenenergies_all.append(energies)


with open("Outputs/dmi_and_aniso/GHz_ev_dmi_only.txt", "r") as file:
    for line in file:
        line = line.strip()
        line = line.split(",")
        energies = [float(e) for e in line[4:]]
        eigenenergies_dmi.append(energies)




assert len(eigenenergies_all) == len(eigenenergies_dmi)

x = np.linspace(0, 1, len(eigenenergies_all))


plt.figure(figsize=(8/2.52, 10/2.52))


ax_main = plt.gca()


ax_main.plot([], [], color='tab:red', label='DM only')
ax_main.plot([], [], color='tab:blue' , label='DM and anisotropy', linestyle='--')


for mode in range(4):
    # Plot main data
    ax_main.plot(x, [np.abs(ev[mode]) for ev in eigenenergies_all], color='tab:blue', linestyle='--')
    ax_main.plot(x, [np.abs(ev[mode]) for ev in eigenenergies_dmi], color='tab:red')


# Create inset axes
ax_zoom = plt.axes([0.3, 0.35, 0.4, 0.4])  # [left, bottom, width, height]

for mode in range(4):

    # Plot zoomed-in data

    ax_zoom.plot(x, [np.abs(ev[mode]) for ev in eigenenergies_dmi], color='tab:red', lw=1.5, alpha=0.85)
    ax_zoom.plot(x, [np.abs(ev[mode]) for ev in eigenenergies_all], color='tab:blue', linestyle='--' , lw=1.5, alpha=0.85)



# Set limits for the zoomed-in plot
ax_zoom.set_xlim(0.04, 0.07)
ax_zoom.set_ylim(1.51, 3.45)

# Remove tick labels on x-axis for the zoomed-in plot
#ax_zoom.set_xticklabels([])


#ax_main.indicate_inset_zoom(ax_zoom, edgecolor="black")

ax_main.set_xticks([0, 1])
ax_main.set_xticklabels([r'$\Gamma$', r'$H$'])

ax_main.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
ax_main.set_ylabel('dispersion relation (meV)')
plt.tight_layout()
plt.savefig("scripts/Figures/dmi_and_aniso.png", dpi=500)
plt.show() 
plt.clf()



for mode in range(4):

    plt.scatter(x, [np.abs(eigenenergies_all[i][mode]-eigenenergies_dmi[i][mode]) for i in range(len(eigenenergies_all))], s=1, c="tab:blue")
    #plt.plot(x, [np.abs(eigenenergies_all[i][mode]-eigenenergies_dmi[i][mode]) for i in range(len(eigenenergies_all))])


plt.ylabel('difference in eigenenergies (meV)')
plt.xticks([0, 1], labels=[r'$\Gamma$', r'$H$'])
plt.tight_layout()
plt.savefig("scripts/Figures/dmi_and_aniso_diff.png", dpi=500)
plt.show()






