from models.polygon import Polygon
from misc.plotting import plotting_setup
from misc.polygon_drawer import Polygon_Drawer
import mpmath as mp
import json


def main():

    # If running into areas with no solutions, you may need to increase the digit precision
    mp.mp.dps = 100

    # The setup for the scale of the axis and the amount of delay between each line/disk drawn
    # The existing instances are built on a wide variety of scales, so you may need to adjust the xlim and ylim
    plotting_setup(xlim=(0, 20), ylim=(0, 20), delay_amnt=0.01)

    # To run existing instances, you can uncomment the following lines
    # f = open("../polygons/pumpkin.json")
    # js = json.load(f)
    # f.close()
    # p = Polygon(**js)
    # Polygon_Drawer(polygon=p)

    # To run the drawer and be able to create your own instances, use the following code.
    # When starting a triangulation, your current instance will be saved as a json file in the polygons folder.
    Polygon_Drawer()


main()
