#!/usr/bin/env python3
"""Class for parsing BAGEL (.bar) output."""

from .MsvUtil import Table

from math import log, pi, sqrt, atan2

class BagelData:
    def __init__(self, fp = None):
        if(fp is None):
            # Default = Diane's G217B Y/M/C direct hybes
            fp = open("/home/mvoorhie/data/Diane/Arrays_4_24_2009/BAGEL2/G217B4.29.09.default.42.bar")

        line = ""
        while(not line.startswith("Unique ID")):
            line = next(fp)

        header = line.rstrip().split("\t")
        rows = [i.rstrip().split("\t") for i in fp]

        self.data = Table(header = header,
                          rows = [i for i in rows if(len(i) == len(header))])
        self.uid2row = dict((i["Unique ID"], i) for i in self.data)

    def __getitem__(self, uid):
        return self.uid2row[uid]

    def __len__(self):
        return len(self.uid2row)

    def polar(self, uid):
        """Return morphological phase in polar coordinates:
        theta: angle, in radians, from Y = 0 through M = (2/3)pi,
            C = (4/3)pi, and Y = 2pi.
        r: radius, equal to log2 of the geometric mean of relative ratios
        """
        j = self.uid2row[uid]
        log2 = log(2)
        sqrt3 = sqrt(3)
        Y = log(float(j["Y"]))/log2
        M = log(float(j["M"]))/log2
        C = log(float(j["C"]))/log2
        x = sqrt3*(M-C)/2.0
        y = Y - (M+C)/2.0
        theta = (pi/2.0 - atan2(y,x)) % (2*pi)
        r = sqrt(x**2+y**2)

        return (theta, r)

    def hexcode(self, uid):
        (theta, r) = self.polar(uid)
        
        s = min(1.0, r/2.0)

        # HSL->RGB conversion based on
        #   http://en.wikipedia.org/wiki/HSV_color_space#From_HSL

        hprime = theta*3.0/pi
        lightness = .5
        # assume lightness = .5
        chroma = s
        X = chroma*(1.0 - abs((hprime % 2) - 1))

        if(hprime < 1):
            rgb1 = (chroma, X, 0)
        elif(hprime < 2):
            rgb1 = (X, chroma, 0)
        elif(hprime < 3):
            rgb1 = (0, chroma, X)
        elif(hprime < 4):
            rgb1 = (0, X, chroma)
        elif(hprime < 5):
            rgb1 = (X, 0, chroma)
        else:
            rgb1 = (chroma, 0, X)

        m = lightness - .5*chroma
        rgb2 = tuple(c+m for c in rgb1)

        hexcode = "#%02x%02x%02x" % (
            tuple(min(255,int(255*c)) for c in rgb2))

        return hexcode

if(__name__ == "__main__"):
    print("Hello, world")
