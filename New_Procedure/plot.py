
def fitsubplotter(x, y, colplot, rowplot, xlabel="x", ylabel="y", formatfile="png", dpi=1200, linestyle="solid", marker="2"):
    """
    function that produce a subplot of a figure
    with resepct to the number of rows and columns
    passed (rowplot and colplot). dpi is the 
    resulution of the image (default=1200dpi).
    x and y are the data to be plotted.
    they mut be lists and y length must be 
    two times x length.
    """
    import matplotlib.pyplot as plt

    plt.figure()

    plt.subplots_adjust(
        left   = 0.1,
        bottom = 0.1,
        right  = 0.95,
        top    = 0.95,
        wspace = 0.3,
        hspace = 0.3
        )

    for i in range(0, (rowplot+colplot)):
        
        plt.subplot(rowplot, colplot,i+1)
        plt.plot(x[i], y[i], label = 'data', color = 'black', marker = marker, linestyle = 'None')
        plt.plot(x[i], y[i+1], label = 'fit line', color = 'grey', linestyle = linestyle, linewidth = 1)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(gridbool)

    plt.savefig(filename, format=formatfile, dpi=resolution)
    plt.close()


