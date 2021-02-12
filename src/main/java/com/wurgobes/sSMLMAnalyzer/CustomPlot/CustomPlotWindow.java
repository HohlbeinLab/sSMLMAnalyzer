package com.wurgobes.sSMLMAnalyzer.CustomPlot;


import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.awt.datatransfer.*;
import java.util.*;
import ij.*;
import ij.gui.*;
import ij.process.*;
import ij.util.*;
import ij.measure.*;
import ij.io.SaveDialog;

/** This class implements the Analyze/Plot Profile command.
 * @author Michael Schmid
 * @author Wayne Rasband
 */
public class CustomPlotWindow extends PlotWindow {

    /** @deprecated */
    public static final int CIRCLE = Plot.CIRCLE;
    /** @deprecated */
    public static final int X = Plot.X;
    /** @deprecated */
    public static final int BOX = Plot.BOX;
    /** @deprecated */
    public static final int TRIANGLE = Plot.TRIANGLE;
    /** @deprecated */
    public static final int CROSS = Plot.CROSS;
    /** @deprecated */
    public static final int LINE = Plot.LINE;
    /** Write first X column when listing or saving. */
    public static boolean saveXValues = true;
    /** Automatically close window after saving values. To set, use Edit/Options/Plots. */
    public static boolean autoClose;
    /** Display the XY coordinates in a separate window. To set, use Edit/Options/Plots. */
    public static boolean listValues;
    /** Interpolate line profiles. To set, use Edit/Options/Plots or setOption("InterpolateLines",boolean). */
    public static boolean interpolate = true;
    // default values for new installations; values will be then saved in prefs
    private static final int WIDTH = 600;
    private static final int HEIGHT = 340;
    private static int defaultFontSize = 14;
    /** The width of the plot (without frame) in pixels. */
    public static int plotWidth = WIDTH;
    /** The height of the plot in pixels. */
    public static int plotHeight = HEIGHT;
    /** The plot text size, can be overridden by Plot.setFont, Plot.setFontSize, Plot.setXLabelFont etc. */
    public static int fontSize = defaultFontSize;
    /** Have axes with no grid lines. If both noGridLines and noTicks are true,
     *	only min&max value of the axes are given */
    public static boolean noGridLines;
    /** Have axes with no ticks. If both noGridLines and noTicks are true,
     *	only min&max value of the axes are given */
    public static boolean noTicks;


    private static final String PREFS_WIDTH = "pp.width";
    private static final String PREFS_HEIGHT = "pp.height";
    private static final String PREFS_FONT_SIZE = "pp.fontsize";
    private static final String OPTIONS = "pp.options";
    private static final int SAVE_X_VALUES = 1;
    private static final int AUTO_CLOSE = 2;
    private static final int LIST_VALUES = 4;
    private static final int INTERPOLATE = 8;
    private static final int NO_GRID_LINES = 16;
    private static final int NO_TICKS = 32;
    private static String moreButtonLabel = "More "+'\u00bb';
    private static String dataButtonLabel = "Data "+'\u00bb';

    boolean wasActivated;			// true after window has been activated once, needed by PlotCanvas

    private Button list, data, more, live;
    private PopupMenu dataPopupMenu, morePopupMenu;
    private static final int NUM_MENU_ITEMS = 22; //how many menu items we have in total
    private MenuItem[] menuItems = new MenuItem[NUM_MENU_ITEMS];
    private Label statusLabel;
    private String userStatusText;
    private static String defaultDirectory = null;
    private static int options;
    private int defaultDigits = -1;
    private int markSize = 5;
    private static CustomPlot staticPlot;
    private CustomPlot plotCustom;

    private PlotMaker plotMaker;
    private ImagePlus srcImp;		// the source image for live plotting
    private Thread bgThread;		// thread for plotting (in the background)
    private boolean doUpdate;		// tells the background thread to update

    private Roi[] rangeArrowRois;	// the overlays (arrows etc) for changing the range. Note: #10-15 must correspond to CustomPlotDialog.dialogType!
    private boolean rangeArrowsVisible;
    private int activeRangeArrow = -1;
    private static Color inactiveRangeArrowColor = Color.GRAY;
    private static Color inactiveRangeRectColor = new Color(0x20404040, true); //transparent gray
    private static Color activeRangeArrowColor = Color.RED;
    private static Color activeRangeRectColor = new Color(0x18ff0000, true); //transparent red

    // static initializer
    static {
        options = Prefs.getInt(OPTIONS, SAVE_X_VALUES);
        autoClose = (options&AUTO_CLOSE)!=0;
        plotWidth = Prefs.getInt(PREFS_WIDTH, WIDTH);
        plotHeight = Prefs.getInt(PREFS_HEIGHT, HEIGHT);
        defaultFontSize = fontSize = Prefs.getInt(PREFS_FONT_SIZE, defaultFontSize);
        Dimension screen = IJ.getScreenSize();
        if (plotWidth>screen.width && plotHeight>screen.height) {
            plotWidth = WIDTH;
            plotHeight = HEIGHT;
        }
    }



    /**
     * @deprecated
     * replaced by the Plot class.
     */
    //public CustomPlotWindow(String title, String xLabel, String yLabel, float[] xValues, float[] yValues) {
    //    super(createImage(title, xLabel, yLabel, xValues, yValues));
    //    plot = staticPlot;
    //    ((PlotCanvas)getCanvas()).setPlot(plot);
    //}

    /**
     * @deprecated
     * replaced by the Plot class.
     */
    //public CustomPlotWindow(String title, String xLabel, String yLabel, double[] xValues, double[] yValues) {
    //    this(title, xLabel, yLabel, Tools.toFloat(xValues), Tools.toFloat(yValues));
    //}

    /** Creates a PlotWindow from a given ImagePlus with a Plot object.
     *  (called when reading an ImagePlus with an associated plot from a file) */

   // public CustomPlotWindow(ImagePlus imp, CustomPlot plot) {
   //     super(imp);
   //     ((PlotCanvas)getCanvas()).setPlot(plot);
   //     this.plot = plot;
   //     draw();
   // }


    CustomPlotWindow(CustomPlot plot){
        super(plot.getImagePlus(), plot);
        super.setPlot(plot);
        //this.plotCustom = plot;

    }

    /** Called by the constructor to generate the image the plot will be drawn on.
     This is a static method because constructors cannot call instance methods. */
    static ImagePlus createImage(String title, String xLabel, String yLabel, float[] xValues, float[] yValues) {
        staticPlot = new CustomPlot(title, xLabel, yLabel, xValues, yValues, null);
        return new ImagePlus(title, staticPlot.getBlankProcessor());
    }


    /** Sets the x-axis and y-axis range.
     *  @deprecated use the corresponding method of the Plot class */
    public void setLimits(double xMin, double xMax, double yMin, double yMax) {
        plotCustom.setLimits(xMin, xMax, yMin, yMax);
    }

    /** Adds a set of points to the plot or adds a curve if shape is set to LINE.
     *	Note that there are more options available by using the methods of the Plot class instead.
     *  @param x			the x-coodinates
     *  @param y			the y-coodinates
     *  @param shape		Plot.CIRCLE, X, BOX, TRIANGLE, CROSS, LINE etc.
     *  @deprecated use the corresponding method of the Plot class */
    public void addPoints(float[] x, float[] y, int shape) {
        plotCustom.addPoints(x, y, shape);
    }

    /** Adds a set of points to the plot using double arrays.
     *	Must be called before the plot is displayed.
     *	Note that there are more options available by using the methods of the Plot class instead.
     *  @deprecated use the corresponding method of the Plot class */
    public void addPoints(double[] x, double[] y, int shape) {
        addPoints(Tools.toFloat(x), Tools.toFloat(y), shape);
    }

    /** Adds vertical error bars to the plotCustom.
     *	Must be called before the plot is displayed.
     *	Note that there are more options available by using the methods of the Plot class instead.
     *  @deprecated use the corresponding method of the Plot class */
    public void addErrorBars(float[] errorBars) {
        plotCustom.addErrorBars(errorBars);
    }

    /** Draws a label.
     *	Note that there are more options available by using the methods of the Plot class instead.
     *  @deprecated use the corresponding method of the Plot class */
    public void addLabel(double x, double y, String label) {
        plotCustom.addLabel(x, y, label);
    }

    /** Changes the drawing color. The frame and labels are
     *	always drawn in black.
     *	Must be called before the plot is displayed.
     *	Note that there are more options available by using the methods of the Plot class instead.
     *  @deprecated use the corresponding method of the Plot class */
    public void setColor(Color c) {
        plotCustom.setColor(c);
    }

    /** Changes the line width.
     *	Must be called before the plot is displayed.
     *	Note that there are more options available by using the methods of the Plot class instead.
     *  @deprecated use the corresponding method of the Plot class */
    public void setLineWidth(int lineWidth) {
        plotCustom.setLineWidth(lineWidth);
    }

    /** Changes the font.
     *	Must be called before the plot is displayed.
     *	Note that there are more options available by using the methods of the Plot class instead.
     *  @deprecated use the corresponding method of the Plot class */
    public void changeFont(Font font) {
        plotCustom.changeFont(font);
    }

    /** Displays the plotCustom. */
    public void draw() {
        //super.draw();

        plotCustom = (CustomPlot) super.getPlot();
        Panel bottomPanel = new Panel();
        int hgap = IJ.isMacOSX()?1:5;

        list = new Button(" List ");
        list.addActionListener(this);
        bottomPanel.add(list);
        bottomPanel.setLayout(new FlowLayout(FlowLayout.RIGHT,hgap,0));
        data = new Button(dataButtonLabel);
        data.addActionListener(this);
        bottomPanel.add(data);
        more = new Button(moreButtonLabel);
        more.addActionListener(this);
        bottomPanel.add(more);
        if (plotCustom!=null && plotCustom.getPlotMaker()!=null) {
            live = new Button("Live");
            live.addActionListener(this);
            bottomPanel.add(live);
        }
        statusLabel = new Label();
        statusLabel.setFont(new Font("Monospaced", Font.PLAIN, 12));
        statusLabel.setBackground(new Color(220, 220, 220));
        bottomPanel.add(statusLabel);
        add(bottomPanel);
        data.add(getDataPopupMenu());
        more.add(getMorePopupMenu());

        plotCustom.menuItems = menuItems;
        plotCustom.draw();
        LayoutManager lm = getLayout();
        if (lm instanceof ImageLayout)
            ((ImageLayout)lm).ignoreNonImageWidths(true);  //don't expand size to make the panel fit
        GUI.scale(bottomPanel);
        maximizeCoordinatesLabelWidth();
        pack();

        ImageProcessor ip = plotCustom.getProcessor();
        boolean ipIsColor = ip instanceof ColorProcessor;
        boolean impIsColor = imp.getProcessor() instanceof ColorProcessor;
        if (ipIsColor != impIsColor)
            imp.setProcessor(null, ip);
        else
            imp.updateAndDraw();
        if (listValues)
            showList(/*useLabels=*/false);
        else
            ic.requestFocus();	//have focus on the canvas, not the button, so that pressing the space bar allows panning
    }

    /** Sets the Plot object shown in this PlotWindow. Does not update the window. */
    public void setPlot(CustomPlot plot) {
        this.plotCustom = plot;
        ((PlotCanvas)getCanvas()).setPlot(plot);
    }

    /** Releases the resources used by this PlotWindow */
    public void dispose() {
        if (plotCustom!=null)
            plotCustom.dispose();
        disableLivePlot();
        plotCustom = null;
        plotMaker = null;
        srcImp = null;
        super.dispose();
    }

    /** Called when the window is activated (WindowListener)
     *  Window layout is finished at latest a few millisec after windowActivated, then the
     *  'wasActivated' boolean is set to tell the ImageCanvas that resize events should
     *  lead to resizing the canvas (before, creating the layout can lead to resize events)*/
    public void windowActivated(WindowEvent e) {
        super.windowActivated(e);
        if (!wasActivated) {
            new Thread(new Runnable() {
                public void run() {
                    IJ.wait(50);  //sometimes, window layout is done only a few millisec after windowActivated
                    wasActivated = true;
                }
            }).start();
        }
    }

    /** Called when the canvas is resized */

    void canvasResized() {
        if (plotCustom == null) return;
		/*Dimension d1 = getExtraSize();
		Dimension d2 = plotCustom.getMinimumSize();
		setMinimumSize(new Dimension(d1.width + d2.width, d1.height + d2.height));*/
        maximizeCoordinatesLabelWidth();
    }

    /** Maximizes the width for the coordinate&status readout field and its parent bottomPanel */
    void maximizeCoordinatesLabelWidth() {
        Insets insets = getInsets();                    //by default, left & right insets are 0 anyhow
        Component parent = statusLabel.getParent();     //the bottomPanel, has insets of 0
        if (!parent.isValid()) parent.validate();
        int cWidth = getWidth() - 2*HGAP - statusLabel.getX() - insets.left - insets.right;
        int cHeight = statusLabel.getPreferredSize().height;
        statusLabel.setPreferredSize(new Dimension(cWidth, cHeight));
        parent.setSize(getWidth() - 2*HGAP, parent.getHeight());
    }

    /** Shows the text in the coordinate&status readout field at the bottom.
     *  This text may get temporarily replaced for 'tooltips' (mouse over range arrows etc.).
     *  Call with a null argument to enable coordinate readout again. */
    public void showStatus(String text) {
        userStatusText = text;
        if (statusLabel != null)
            statusLabel.setText(text == null ? "" : text);
    }

    /** Names for popupMenu items. Update NUM_MENU_ITEMS at the top when adding new ones! */
    private static int SAVE=0, COPY=1, COPY_ALL=2, LIST_SIMPLE=3, ADD_FROM_TABLE=4, ADD_FROM_PLOT=5, ADD_FIT=6, //data menu
            SET_RANGE=7, PREV_RANGE=8, RESET_RANGE=9, FIT_RANGE=10,  //the rest is in the more menu
            ZOOM_SELECTION=11, AXIS_OPTIONS=12, LEGEND=13, STYLE=14, TEMPLATE=15, RESET_PLOT=16,
            FREEZE=17, HI_RESOLUTION=18, PROFILE_PLOT_OPTIONS=19, RESCALE_LUT=20, AREA_HIST=21;
    //the following commands are disabled when the plot is frozen
    private static int[] DISABLED_WHEN_FROZEN = new int[]{ADD_FROM_TABLE, ADD_FROM_PLOT, ADD_FIT,
            SET_RANGE, PREV_RANGE, RESET_RANGE, FIT_RANGE, ZOOM_SELECTION, AXIS_OPTIONS, LEGEND, STYLE, RESET_PLOT};

    /** Prepares and returns the popupMenu of the Data>> button */
    PopupMenu getDataPopupMenu() {
        dataPopupMenu = new PopupMenu();
        GUI.scalePopupMenu(dataPopupMenu);
        MenuItem[] menuItems = new MenuItem[NUM_MENU_ITEMS];
        menuItems[SAVE] = addPopupItem(dataPopupMenu, "Save Data...");
        menuItems[COPY] = addPopupItem(dataPopupMenu, "Copy 1st Data Set");
        menuItems[COPY_ALL] = addPopupItem(dataPopupMenu, "Copy All Data");
        menuItems[LIST_SIMPLE] = addPopupItem(dataPopupMenu, "List (Simple Headings)");
        dataPopupMenu.addSeparator();
        menuItems[ADD_FROM_TABLE] = addPopupItem(dataPopupMenu, "Add from Table...");
        menuItems[ADD_FROM_PLOT] = addPopupItem(dataPopupMenu, "Add from Plot...");
        menuItems[ADD_FIT] = addPopupItem(dataPopupMenu, "Add Fit...");
        this.menuItems = menuItems;
        return dataPopupMenu;
    }

    /** Prepares and returns the popupMenu of the More>> button */
    PopupMenu getMorePopupMenu() {
        morePopupMenu = new PopupMenu();
        GUI.scalePopupMenu(morePopupMenu);
        menuItems[SET_RANGE] = addPopupItem(morePopupMenu, "Set Range...");
        menuItems[PREV_RANGE] = addPopupItem(morePopupMenu, "Previous Range");
        menuItems[RESET_RANGE] = addPopupItem(morePopupMenu, "Reset Range");
        menuItems[FIT_RANGE] = addPopupItem(morePopupMenu, "Set Range to Fit All");
        menuItems[ZOOM_SELECTION] = addPopupItem(morePopupMenu, "Zoom to Selection");
        morePopupMenu.addSeparator();
        menuItems[AXIS_OPTIONS] = addPopupItem(morePopupMenu, "Axis Options...");
        menuItems[LEGEND] = addPopupItem(morePopupMenu, "Legend...");
        menuItems[STYLE] = addPopupItem(morePopupMenu, "Contents Style...");
        menuItems[TEMPLATE] = addPopupItem(morePopupMenu, "Use Template...");
        menuItems[RESET_PLOT] = addPopupItem(morePopupMenu, "Reset Format");
        menuItems[FREEZE] = addPopupItem(morePopupMenu, "Freeze Plot", true);
        menuItems[HI_RESOLUTION] = addPopupItem(morePopupMenu, "High-Resolution Plot...");
        morePopupMenu.addSeparator();
        menuItems[PROFILE_PLOT_OPTIONS] = addPopupItem(morePopupMenu, "Plot Defaults...");
        morePopupMenu.addSeparator();
        menuItems[RESCALE_LUT] = addPopupItem(morePopupMenu, "Rescale LUT...");
        menuItems[AREA_HIST] = addPopupItem(morePopupMenu, "Create Hist from selection");
        return morePopupMenu;
    }

    MenuItem addPopupItem(PopupMenu popupMenu, String s) {
        return addPopupItem(popupMenu, s, false);
    }

    MenuItem addPopupItem(PopupMenu popupMenu, String s, boolean isCheckboxItem) {
        MenuItem mi = null;
        if (isCheckboxItem) {
            mi = new CheckboxMenuItem(s);
            ((CheckboxMenuItem)mi).addItemListener(this);
        } else {
            mi = new MenuItem(s);
            mi.addActionListener(this);
        }
        popupMenu.add(mi);
        return mi;
    }

    /** Called if user has activated a button or popup menu item */
    public void actionPerformed(ActionEvent e) {
        menuItems = plotCustom.menuItems;
        try {
            Object b = e.getSource();
            if (b==live)
                toggleLiveProfiling();
            else if (b==list)
                showList(/*useLabels=*/true);
            else if (b==data) {
                enableDisableMenuItems();
                dataPopupMenu.show((Component)b, 1, 1);
            } else if (b==more) {
                enableDisableMenuItems();
                morePopupMenu.show((Component)b, 1, 1);
            } else if (b==menuItems[SAVE])
                saveAsText();
            else if (b==menuItems[COPY])
                copyToClipboard(false);
            else if (b==menuItems[COPY_ALL])
                copyToClipboard(true);
            else if (b==menuItems[LIST_SIMPLE])
                showList(/*useLabels=*/false);
            else if (b==menuItems[ADD_FROM_TABLE])
                new PlotContentsDialog(plotCustom, PlotContentsDialog.ADD_FROM_TABLE).showDialog(this);
            else if (b==menuItems[ADD_FROM_PLOT])
                new PlotContentsDialog(plotCustom, PlotContentsDialog.ADD_FROM_PLOT).showDialog(this);
            else if (b==menuItems[ADD_FIT])
                new PlotContentsDialog(plotCustom, PlotContentsDialog.ADD_FIT).showDialog(this);
            else if (b==menuItems[ZOOM_SELECTION]) {
                if (imp!=null && imp.getRoi()!=null && imp.getRoi().isArea())
                    plotCustom.zoomToRect(imp.getRoi().getBounds());
            } else if (b==menuItems[SET_RANGE])
                new CustomPlotDialog(plotCustom, CustomPlotDialog.SET_RANGE).showDialog(this);
            else if (b==menuItems[PREV_RANGE])
                plotCustom.setPreviousMinMax();
            else if (b==menuItems[RESET_RANGE])
                plotCustom.setLimitsToDefaults(true);
            else if (b==menuItems[FIT_RANGE])
                plotCustom.setLimitsToFit(true);
            else if (b==menuItems[AXIS_OPTIONS])
                new CustomPlotDialog(plotCustom, CustomPlotDialog.AXIS_OPTIONS).showDialog(this);
            else if (b==menuItems[LEGEND])
                new CustomPlotDialog(plotCustom, CustomPlotDialog.LEGEND).showDialog(this);
            else if (b==menuItems[STYLE])
                new PlotContentsDialog(plotCustom, PlotContentsDialog.STYLE).showDialog(this);
            else if (b==menuItems[TEMPLATE])
                new CustomPlotDialog(plotCustom, CustomPlotDialog.TEMPLATE).showDialog(this);
            else if (b==menuItems[RESET_PLOT]) {
                plotCustom.setFont(Font.PLAIN, fontSize);
                plotCustom.setAxisLabelFont(Font.PLAIN, fontSize);
                plotCustom.setFormatFlags(Plot.getDefaultFlags());
                plotCustom.setFrameSize(plotWidth, plotHeight); //updates the image only when size changed
                plotCustom.updateImage();
            } else if (b==menuItems[HI_RESOLUTION])
                new CustomPlotDialog(plotCustom, CustomPlotDialog.HI_RESOLUTION).showDialog(this);
            else if (b==menuItems[PROFILE_PLOT_OPTIONS])
                IJ.doCommand("Plots...");
            else if (b==menuItems[RESCALE_LUT]) {
                new CustomPlotDialog(plotCustom, CustomPlotDialog.LUT_CHANGE, plotCustom.ownColorTable).showDialog(this);
            } else if (b==menuItems[AREA_HIST]) {
                //RoiManager roiManager = new RoiManager();
                Roi roi = imp.getRoi();
                if(roi == null || !roi.isArea()) {
                    IJ.showMessage("No region is selected");
                } else {
                   plotCustom.makeAreaHist(roi);
                }

            } else {
                System.out.println("dont see this");
            }

            ic.requestFocus();	//have focus on the canvas, not the button, so that pressing the space bar allows panning
        } catch (Exception ex) { IJ.handleException(ex); }
    }

    private void enableDisableMenuItems() {
        menuItems = plotCustom.menuItems;

        boolean frozen = plotCustom.isFrozen();	//prepare menu according to 'frozen' state of plot


        ((CheckboxMenuItem)menuItems[FREEZE]).setState(frozen);
        for (int i : DISABLED_WHEN_FROZEN)
            menuItems[i].setEnabled(!frozen);
        if (!PlotContentsDialog.tableWindowExists())
            menuItems[ADD_FROM_TABLE].setEnabled(false);
        if (plotCustom.getDataObjectDesignations().length == 0)
            menuItems[ADD_FIT].setEnabled(false);
    }

    /** Called if the user activates/deactivates a CheckboxMenuItem */
    public void itemStateChanged(ItemEvent e) {
        if (e.getSource()==menuItems[FREEZE]) {
            boolean frozen = ((CheckboxMenuItem)menuItems[FREEZE]).getState();
            plotCustom.setFrozen(frozen);
        }
    }

    /**
     * Updates the X and Y values when the mouse is moved and, if appropriate,
     * shows/hides the overlay with the triangular buttons for changing the axis
     * range limits.
     * Overrides mouseMoved() in ImageWindow.
     *
     * @see ij.gui.ImageWindow#mouseMoved
     */
    public void mouseMoved(int x, int y) {
        super.mouseMoved(x, y);
        if (plotCustom == null)
            return;
        String statusText = null; //coordinate readout, status or tooltip, will be shown in coordinate&status line

        //arrows and other symbols for modifying the plot range
        if (x < plotCustom.leftMargin || y > plotCustom.topMargin + plotCustom.frameHeight) {
            if (!rangeArrowsVisible && !plotCustom.isFrozen())
                showRangeArrows();
            if (activeRangeArrow < 0)       //mouse is not on one of the symbols, ignore (nothing to display)
            {}
            else if (activeRangeArrow < 8)  //mouse over an arrow: 0,3,4,7 for increase, 1,2,5,6 for decrease
                statusText = ((activeRangeArrow+1)&0x02) != 0 ? "Decrease Range" : "Increase Range";
            else if (activeRangeArrow == 8) //it's the 'R' icon
                statusText = "Reset Range";
            else if (activeRangeArrow == 9) //it's the 'F' icon
                statusText = "Full Range (Fit All)";
            else if (activeRangeArrow >= 10 &&
                    activeRangeArrow < 14)  //space between arrow-pairs for single number
                statusText = "Set limit...";
            else if (activeRangeArrow >= 14)
                statusText = "Axis Range & Options...";
            boolean repaint = false;
            if (activeRangeArrow >= 0 && !rangeArrowRois[activeRangeArrow].contains(x, y)) {
                rangeArrowRois[activeRangeArrow].setFillColor(
                        activeRangeArrow < 10 ? inactiveRangeArrowColor : inactiveRangeRectColor);
                repaint = true;             //de-highlight arrow where cursor has moved out
                activeRangeArrow = -1;
            }
            if (activeRangeArrow < 0) {     //no currently highlighted arrow, do we have a new one?
                int i = getRangeArrowIndex(x, y);
                if (i >= 0) {               //we have an arrow or symbol at cursor position
                    rangeArrowRois[i].setFillColor(
                            i < 14 ? activeRangeArrowColor : activeRangeRectColor);
                    activeRangeArrow = i;
                    repaint = true;
                }
            }
            if (repaint) ic.repaint();
        } else if (rangeArrowsVisible)
            hideRangeArrows();

        if (statusText == null)
            try{
            statusText = userStatusText != null ? userStatusText : plotCustom.getCoordinates(x, y);
            } catch (Exception ignored){}
        if (statusLabel != null)
            statusLabel.setText(statusText);
    }

    /** Called by PlotCanvas */
    void mouseExited(MouseEvent e) {
        if (rangeArrowsVisible)
            hideRangeArrows();
    }

    /** Mouse wheel: zooms when shift or ctrl is pressed, scrolls in x if space bar down, in y otherwise. */
    public synchronized void mouseWheelMoved(MouseWheelEvent e) {
        if (plotCustom.isFrozen() || !(ic instanceof PlotCanvas)) {	   //frozen plots are like normal images
            super.mouseWheelMoved(e);
            return;
        }

        int rotation = e.getWheelRotation();
        int amount = e.getScrollAmount();
        if (e.getX() < plotCustom.leftMargin || e.getX() > plotCustom.leftMargin + plotCustom.frameWidth)//n__
            return;
        if (e.getY() < plotCustom.topMargin || e.getY() > plotCustom.topMargin + plotCustom.frameHeight)
            return;
        boolean ctrl = (e.getModifiers()&Event.CTRL_MASK)!=0;
        if (amount<1) amount=1;
        if (rotation==0)
            return;
        if (ctrl||IJ.shiftKeyDown()) {
            double zoomFactor = rotation<0 ? Math.pow(2, 0.2) : Math.pow(0.5, 0.2);
            Point loc = ic.getCursorLoc();
            int x = ic.screenX(loc.x);
            int y = ic.screenY(loc.y);
            ((CustomPlot)((PlotCanvas)ic).getPlot()).zoom(x, y, zoomFactor);
        } else if (IJ.spaceBarDown())
            plotCustom.scroll(rotation*amount*Math.max(ic.getWidth()/50, 1), 0);
        else
            plotCustom.scroll(0, rotation*amount*Math.max(ic.getHeight()/50, 1));
    }

    /**
     * Creates an overlay with triangular buttons and othr symbols for changing the axis range
     * limits and shows it
     */
    void showRangeArrows() {
        if (imp == null)
            return;
        hideRangeArrows(); //in case we have old arrows from a different plot size or so
        rangeArrowRois = new Roi[4 * 2 + 2 + 4 + 2]; //4 arrows per axis, + 'Reset' and 'Fit All' icons, + 4 numerical input boxes + 2 axes
        int i = 0;
        int height = imp.getHeight();
        int arrowH = plotCustom.topMargin < 14 ? 6 : 8; //height of arrows and distance between them; base is twice that value
        float[] yP = new float[]{height - arrowH / 2, height - 3 * arrowH / 2, height - 5 * arrowH / 2 - 0.1f};

        for (float x : new float[]{plotCustom.leftMargin, plotCustom.leftMargin + plotCustom.frameWidth}) { //create arrows for x axis
            float[] x0 = new float[]{x - arrowH / 2, x - 3 * arrowH / 2 - 0.1f, x - arrowH / 2};
            rangeArrowRois[i++] = new PolygonRoi(x0, yP, 3, Roi.POLYGON);
            float[] x1 = new float[]{x + arrowH / 2, x + 3 * arrowH / 2 + 0.1f, x + arrowH / 2};
            rangeArrowRois[i++] = new PolygonRoi(x1, yP, 3, Roi.POLYGON);
        }
        float[] xP = new float[]{arrowH / 2 - 0.1f, 3 * arrowH / 2, 5 * arrowH / 2 + 0.1f};
        for (float y : new float[]{plotCustom.topMargin + plotCustom.frameHeight, plotCustom.topMargin}) { //create arrows for y axis
            float[] y0 = new float[]{y + arrowH / 2, y + 3 * arrowH / 2 + 0.1f, y + arrowH / 2};
            rangeArrowRois[i++] = new PolygonRoi(xP, y0, 3, Roi.POLYGON);
            float[] y1 = new float[]{y - arrowH / 2, y - 3 * arrowH / 2 - 0.1f, y - arrowH / 2};
            rangeArrowRois[i++] = new PolygonRoi(xP, y1, 3, Roi.POLYGON);
        }
        Font theFont = new Font("SansSerif", Font.BOLD, 13);

        TextRoi txtRoi = new TextRoi(1, height - 19, "\u2009R\u2009", theFont);  //thin spaces to make roi slightly wider
        rangeArrowRois[8] = txtRoi;
        TextRoi txtRoi2 = new TextRoi(20, height - 19, "\u2009F\u2009", theFont);
        rangeArrowRois[9] = txtRoi2;

        rangeArrowRois[10] = new Roi(plotCustom.leftMargin - arrowH/2 + 1, height - 5 * arrowH / 2, arrowH - 2, arrowH * 2);//numerical box left
        rangeArrowRois[11] = new Roi(plotCustom.leftMargin + plotCustom.frameWidth - arrowH/2 + 1, height - 5 * arrowH / 2, arrowH - 2, arrowH * 2);//numerical box right
        rangeArrowRois[12] = new Roi(arrowH / 2, plotCustom.topMargin + plotCustom.frameHeight - arrowH/2 + 1, arrowH * 2, arrowH -2);//numerical box bottom
        rangeArrowRois[13] = new Roi(arrowH / 2, plotCustom.topMargin - arrowH/2 + 1,  arrowH * 2, arrowH - 2   );//numerical box top

        int topMargin = plotCustom.topMargin;
        int bottomMargin = topMargin + plotCustom.frameHeight;
        int leftMargin = plotCustom.leftMargin;
        int rightMargin = plotCustom.leftMargin + plotCustom.frameWidth;
        rangeArrowRois[14] = new Roi(leftMargin, bottomMargin+2,        // area to click for x axis options
                rightMargin - leftMargin + 1, 2*arrowH);
        rangeArrowRois[15] = new Roi(leftMargin-2*arrowH-2, topMargin,  // area to click for y axis options
                2*arrowH, bottomMargin - topMargin + 1);

        Overlay ovly = imp.getOverlay();
        if (ovly == null)
            ovly = new Overlay();
        for (Roi roi : rangeArrowRois) {
            if (roi instanceof PolygonRoi)
                roi.setFillColor(inactiveRangeArrowColor);
            else if (roi instanceof TextRoi) {
                roi.setStrokeColor(Color.WHITE);
                roi.setFillColor(inactiveRangeArrowColor);
            } else
                roi.setFillColor(inactiveRangeRectColor); //transparent gray for single number boxes and axis range
            ovly.add(roi);
        }
        imp.setOverlay(ovly);
        ic.repaint();
        rangeArrowsVisible = true;
    }

    void hideRangeArrows() {
        if (imp == null || rangeArrowRois==null) return;
        Overlay ovly = imp.getOverlay();
        if (ovly == null) return;
        for (Roi roi : rangeArrowRois)
            ovly.remove(roi);
        ic.repaint();
        rangeArrowsVisible = false;
        activeRangeArrow = -1;
    }

    /** Returns the index of the range-modifying symbol or axis at the
     *  cursor position x,y, or -1 of none.
     *  Index numbers for arrows start with 0 at the 'down' arrow of the
     *  lower side of the x axis and end with 7 the up arrow at the upper
     *  side of the y axis. Numbers 8 & 9 are for "Reset Range" and "Fit All";
     *  numbers 10-13 for a dialog to set a single limit, and 14-15 for the axis options. */

    int getRangeArrowIndex(int x, int y) {
        if (!rangeArrowsVisible) return -1;
        for (int i=0; i<rangeArrowRois.length; i++)
            if (rangeArrowRois[i].getBounds().contains(x,y))
                return i;
        return -1;
    }


    /** Shows the data of the backing plot in a Textwindow with columns */
    void showList(boolean useLabels){
        ResultsTable rt = plotCustom.getResultsTable(saveXValues, useLabels);
        if (rt==null) return;
        rt.show("Plot Values");
        if (autoClose) {
            imp.changes=false;
            close();
        }
    }

    /** Returns the plot values with simple headings (X, Y, Y1 etc, not the labels) as a ResultsTable.
     *  Use plotCustom.getResultsTableWithLabels for a table with data set labels as column headings */
    public ResultsTable getResultsTable() {
        return plotCustom.getResultsTable(saveXValues);
    }

    /** creates the data that fills the showList() function values */
    private String getValuesAsString(){
        ResultsTable rt = getResultsTable();
        StringBuffer sb = new StringBuffer();
        for (int i=0; i<rt.size(); i++) {
            sb.append(rt.getRowAsString(i));
            sb.append("\n");
        }
        return sb.toString();
    }

    /** Saves the data of the plot in a text file */
    void saveAsText() {
        if (plotCustom.getXValues() == null) {
            IJ.error("Plot has no data");
            return;
        }
        SaveDialog sd = new SaveDialog("Save as Text", "Values", Prefs.defaultResultsExtension());
        String name = sd.getFileName();
        if (name==null) return;
        String directory = sd.getDirectory();
        IJ.wait(250);  // give system time to redraw ImageJ window
        IJ.showStatus("Saving plot values...");
        ResultsTable rt = plotCustom.getResultsTable(/*writeFirstXColumn=*/saveXValues, /*useLabels=*/true);
        try {
            rt.saveAs(directory+name);
        } catch (IOException e) {
            IJ.error("" + e);
            return;
        }
        if (autoClose)
        {imp.changes=false; close();}
    }

    /** Copy the first dataset or all values to the clipboard */
    void copyToClipboard(boolean writeAllColumns) {
        float[] xValues = plotCustom.getXValues();
        float[] yValues = plotCustom.getYValues();
        if (xValues == null) return;
        Clipboard systemClipboard = null;
        try {systemClipboard = getToolkit().getSystemClipboard();}
        catch (Exception e) {systemClipboard = null; }
        if (systemClipboard==null)
        {IJ.error("Unable to copy to Clipboard."); return;}
        IJ.showStatus("Copying plot values...");
        CharArrayWriter aw = new CharArrayWriter(10*xValues.length);
        PrintWriter pw = new PrintWriter(aw); //uses platform's line termination characters

        if (writeAllColumns) {
            ResultsTable rt = plotCustom.getResultsTableWithLabels();
            if (!Prefs.dontSaveHeaders) {
                String headings = rt.getColumnHeadings();
                pw.println(headings);
            }
            for (int i=0; i<rt.size(); i++)
                pw.println(rt.getRowAsString(i));
        } else {
            int xdigits = 0;
            if (saveXValues)
                xdigits = plotCustom.getPrecision(xValues);
            int ydigits = plotCustom.getPrecision(yValues);
            for (int i=0; i<Math.min(xValues.length, yValues.length); i++) {
                if (saveXValues)
                    pw.println(IJ.d2s(xValues[i],xdigits)+"\t"+IJ.d2s(yValues[i],ydigits));
                else
                    pw.println(IJ.d2s(yValues[i],ydigits));
            }
        }
        String text = aw.toString();
        pw.close();
        StringSelection contents = new StringSelection(text);
        systemClipboard.setContents(contents, this);
        IJ.showStatus(text.length() + " characters copied to Clipboard");
        if (autoClose)
        {imp.changes=false; close();}
    }

    public void lostOwnership(Clipboard clipboard, Transferable contents) {}

    public float[] getXValues() {
        return plotCustom.getXValues();
    }

    public float[] getYValues() {
        return plotCustom.getYValues();
    }

    /** Draws a new plot in this window. */
    public void drawPlot(CustomPlot plot) {
        this.plotCustom = plot;
        if (imp!=null) {
            if (ic instanceof PlotCanvas)
                ((PlotCanvas)ic).setPlot(plot);
            imp.setProcessor(null, plotCustom.getProcessor());
            plotCustom.setImagePlus(imp); //also adjusts the calibration of imp
        }
    }

    /** Called once when ImageJ quits. */
    public static void savePreferences(Properties prefs) {
        double min = ProfilePlot.getFixedMin();
        double max = ProfilePlot.getFixedMax();
        prefs.put(PREFS_WIDTH, Integer.toString(plotWidth));
        prefs.put(PREFS_HEIGHT, Integer.toString(plotHeight));
        prefs.put(PREFS_FONT_SIZE, Integer.toString(defaultFontSize));
        int options = 0;
        if (!interpolate) options |= INTERPOLATE; // true=0, false=1
        prefs.put(OPTIONS, Integer.toString(options));
    }

    private void toggleLiveProfiling() {
        boolean liveMode = bgThread != null;
        if (liveMode)
            disableLivePlot();
        else
            enableLivePlot();
    }

    /* Enable live plotting.
     * This requires that the PlotWindow has been initialized with a Plot having a PlotMaker */
    private void enableLivePlot() {
        if (plotMaker==null)
            plotMaker = plotCustom!=null?plotCustom.getPlotMaker():null;
        if (plotMaker==null) return;
        srcImp = plotMaker.getSourceImage();
        if (srcImp==null)
            return;
        if (bgThread==null) {
            bgThread = new Thread(this, "Live Plot");
            bgThread.setPriority(Math.max(bgThread.getPriority()-3, Thread.MIN_PRIORITY));
            doUpdate = true;
            bgThread.start();
        }
        if (IJ.debugMode) IJ.log("PlotWindow.createListeners");
        ImagePlus.addImageListener(this);
        Roi.addRoiListener(this);
        Font font = live.getFont();
        live.setFont(new Font(font.getName(), Font.BOLD, font.getSize()));
        live.setForeground(Color.red);
    }

    private void disableLivePlot() {
        if (IJ.debugMode) IJ.log("PlotWindow.disableLivePlot: "+srcImp);
        if (srcImp==null)
            return;
        if (bgThread!=null)
            bgThread.interrupt();
        bgThread = null;
        ImagePlus.removeImageListener(this);
        Roi.removeRoiListener(this);
        if (live != null) {
            Font font = live.getFont();
            live.setFont(new Font(font.getName(), Font.PLAIN, font.getSize()));
            live.setForeground(Color.black);
        }
    }


    /** For live plots, update the plot if the ROI of the source image changes */
    public synchronized void roiModified(ImagePlus img, int id) {
        if (IJ.debugMode) IJ.log("PlotWindow.roiModified: "+img+"  "+id);
        if (img==srcImp) {
            doUpdate=true;
            notify();
        }
    }

    // Unused
    public void imageOpened(ImagePlus imp) {
    }

    /** For live plots, this method is called if the source image content is changed. */
    public synchronized void imageUpdated(ImagePlus imp) {
        if (imp==srcImp) {
            doUpdate = true;
            notify();
        }
    }

    /** For live plots, if either the source image or this image are closed, exit live mode */
    public void imageClosed(ImagePlus imp) {
        if (imp==srcImp || imp==this.imp) {
            disableLivePlot();
            srcImp = null;
            plotMaker = null;
        }
    }

    // the background thread for live plotting.
    public void run() {
        while (true) {
            IJ.wait(50);	//delay to make sure the roi has been updated
            CustomPlot plot = plotMaker!=null? (CustomPlot)plotMaker.getPlot() :null;
            if (doUpdate && plot!=null && plotCustom.getNumPlotObjects()>0) {
                plotCustom.useTemplate(this.plotCustom, this.plotCustom.templateFlags | Plot.COPY_SIZE | Plot.COPY_LABELS | Plot.COPY_AXIS_STYLE |
                        Plot.COPY_CONTENTS_STYLE | Plot.COPY_LEGEND | Plot.COPY_EXTRA_OBJECTS);
                plotCustom.setPlotMaker(plotMaker);
                this.plotCustom = plot;
                ((PlotCanvas)ic).setPlot(plot);
                ImageProcessor ip = plotCustom.getProcessor();
                if (ip!=null && imp!=null) {
                    imp.setProcessor(null, ip);
                    plotCustom.setImagePlus(imp);
                }
            }
            synchronized(this) {
                if (doUpdate) {
                    doUpdate = false;		//and loop again
                } else {
                    try {wait();}	//notify wakes up the thread
                    catch(InterruptedException e) { //interrupted tells the thread to exit
                        return;
                    }
                }
            }
        }
    }

    /** Returns the Plot associated with this PlotWindow. */
    public Plot getPlot() {
        return plotCustom;
    }

    /** Freezes the active plot window, so the image does not get redrawn for zooming,
     *  setting the range, etc. */
    public static void freeze() {
        Window win = WindowManager.getActiveWindow();
        if (win!=null && (win instanceof ij.gui.PlotWindow))
            ((CustomPlotWindow)win).getPlot().setFrozen(true);
    }


    public static void setDefaultFontSize(int size) {
        if (size < 9) size = 9;
        if (size > 36) size = 36;
        defaultFontSize = size;
    }

    public static int getDefaultFontSize() {
        return defaultFontSize;
    }

}
