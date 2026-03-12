from matplotlib.pyplot import figure, plot, subplots
from decorator import decorator
@decorator
def pca_plot(f, self, *args, **kw):
    ax = args[0]
    fig = None
    if(ax is None):
        (fig, ax) = subplots(1)
    f(self,ax,*args[1:],**kw)
    return fig