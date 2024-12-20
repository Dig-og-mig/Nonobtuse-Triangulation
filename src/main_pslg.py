from models.pslg import PSLG
from misc.pslg_drawer import PSLG_Drawer
import misc.plotting
import misc.import_pslg
import mpmath as mp
import json


def main():

    # If running into areas with no solutions, you may need to increase the digit precision
    mp.mp.dps = 200

    # The setup for the scale of the axis and the amount of delay between each line/disk drawn
    # The existing instances are built on a wide variety of scales, so you may need to adjust the xlim and ylim
    misc.plotting.plotting_setup(xlim=(0, 20), ylim=(0, 20), delay_amnt=0.01)

    # To run cgshop_instances, you can uncomment the following lines
    # More instances can be found in the cgshop_instances folder
    # Some of these instances are quite complex, and therefore take a long time to run.
    # p = misc.import_pslg.import_pslg(
    #     "../cgshop_instances/simple-polygon-exterior_100_686dd044.instance.json")
    # PSLG_Drawer(p=p)

    # To run existing instances, you can uncomment the following lines
    # f = open("../pslgs/3.json")
    # js = json.load(f)
    # f.close()
    # p = PSLG(**js)
    # PSLG_Drawer(p=p)

    # To run the drawer and be able to create your own instances, use the following code.
    # When starting a triangulation, your current instance will be saved as a json file in the pslgs folder.
    PSLG_Drawer()


main()
