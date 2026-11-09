#
#   plot_guess                       
#                               
#   Script to plot an initial flowfield guess created using the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_guess.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the guess data
    filename = 'out_guess_' + sys.argv[-1] + '.bin'
    imagename = 'out_guess_' + sys.argv[-1] + '.jpg'
    inname = 'input_' + sys.argv[-1] + '.txt'

    # Read the case and calculate Mach number using its gas properties
    av = read_settings(inname)
    g = read_case(filename)
    g = calc_secondary(av,g)

    # Open figure window for the four primary variables and Mach number
    fig,ax = plt.subplots(3,2,sharex=True,sharey=True,figsize=[14.4,10.8]);

    # Set subplot aspect ratios as equal and remove axes labels
    ax = ax.flatten()
    for a in ax:
        a.set_aspect('equal',adjustable='box'); a.axis('off')

    # Plot the primary flow variables to show the guess
    fieldnames = ['ro','roe','rovx','rovy','mach']
    colnames = ['ro','roe','rovx','rovy','Mach number']
    for n,name in enumerate(fieldnames):
 
        # Plot filled contour levels
        hc = ax[n].pcolormesh(g['x'],g['y'],g[name],shading='gouraud')

        # Add colorbar with variable name
        colorbar(hc,colnames[n])

        # Draw the walls of the block
        plot_wall(ax[n],g)

    # Hide the unused sixth panel
    ax[-1].set_visible(False)
    fig.tight_layout()

    # Show all the plots
    plt.savefig(imagename)
    plt.show()

    
main()


