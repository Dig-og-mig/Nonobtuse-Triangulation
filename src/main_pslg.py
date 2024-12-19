from pslg import PSLG
from pslg_drawer import PSLG_Drawer
import plotting
import mpmath as mp
import json


def main():

    mp.mp.dps = 200

    plotting.plotting_setup(xlim=(0, 20), ylim=(0, 20), delay_amnt=0.01)
    # To run existing instances, you can uncomment the following lines
    # f = open("pslgs/filename.json")
    # js = json.load(f)
    # f.close()
    # p = pslg.PSLG(**js)
    # PSLG_Drawer(p=p)

    # To run the drawer and be able to create your own instances, use the following code.
    # When starting a triangulation, your current instance will be saved as a json file in the pslgs folder.
    PSLG_Drawer()


main()
