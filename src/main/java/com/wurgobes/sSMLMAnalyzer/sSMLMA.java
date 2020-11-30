package com.wurgobes.sSMLMAnalyzer;



import ij.*;

import ij.gui.HistogramWindow;
import ij.gui.Plot;
import ij.process.FloatProcessor;


import net.imagej.ImageJ;

import net.imagej.lut.LUTService;

import net.imglib2.type.numeric.IntegerType;

import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;


import java.io.IOException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import fiji.util.gui.GenericDialogPlus;

import org.jblas.FloatMatrix;
import org.jblas.exceptions.LapackException;
import org.scijava.plugin.Plugin;


import static com.wurgobes.sSMLMAnalyzer.Util.*;
import static com.wurgobes.sSMLMAnalyzer.levenshtein.getTheClosestMatch;
import static ij.util.ThreadUtil.*;



// test flip z for artifact

// put 1/2th order into seperate file
// make thunderSTORM z plot file

@Plugin(type = Command.class, menuPath = "Plugins>Spectral Analyzer>Analyze Pairs")
public class  sSMLMA < T extends IntegerType<T>> implements Command {


    @Parameter
    private LogService logService;

    @Parameter
    private LUTService lutService;

    //@Parameter
    //private OpService opService;

    //private final int[] unit_decades = {0, 0, 0, -2, -3, -6, -9, -10, -12, -15};

    private final String[] possible_options = {"id", "frame", "x", "y", "z", "intensity", "offset", "bkgstd", "sigma1", "sigma2", "uncertainty", "detections", "chi"};
    private final String[] unit_prefixes = {"null", "photons", "m", "cm", "mm", "um", "nm", "ang", "pm", "fm"};


    private String filePath = "F:\\ThesisData\\output\\output3_drift.csv";
    //private String filePath = "F:\\ThesisData\\output\\combined_drift.csv";
    //private String filePath = "F:\\ThesisData\\output\\4_grating_drift.csv";


    //private String CSV_FILE_NAME = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\grating_cleaned.csv";
    private String CSV_FILE_NAME = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\temp.csv";
    //private String CSV_FILE_NAME = "F:\\ThesisData\\output\\results1.csv";

    private boolean saveSCV = false;

    private float[] angRange = new float[] {0, 0};
    //private final float[] angRange = {(float) (-1 * Math.PI), (float) (-0.95 * Math.PI) }; //more than and less than


    private float[] distRange = new float[] {0, 0}; //default was {1500, 2500}
    //private final float[] distRange = {1940, 2600}; //1800 3000 (1940, 2240)


    private int orders = 4; //This its the number of orders, including the 0th, so 3 would be 0th + 1st + 2nd
    private final int orderCollumns = 6;
    private int totalCollumns = orders * orderCollumns;


    private float binwidth = 2.5f;


    private boolean toCleanup = false;
    private int neighbours = 7;
    private float cleanDistance = 100;

    private boolean processing = true;

    private static boolean debug = false;

    private OwnColorTable ownColorTable;

    private boolean flipAngles = false;

    private static boolean runningFromIDE = false;

    private String defaultLUT = "Spectrum.LUT";


    public int setup(){
        if (debug) return 1;

        GenericDialogPlus gd = new GenericDialogPlus("settings");

        gd.addFileField("CSV input", filePath, 25);

        gd.addCheckbox("Save to CSV?", saveSCV);
        gd.addToSameRow();
        gd.addFileField("CSV output", CSV_FILE_NAME, 25);

        gd.addMessage("------------------------------------------Angles and Distances------------------------------------------------------------------------------------------------------------------------------");

        gd.addMessage("If you want to override the calculated values change the values below.\nWill only overwrite values if they are NOT 0.");


        gd.addNumericField("Start Angle", angRange[0]);
        gd.addToSameRow();
        gd.addNumericField("End Angle", angRange[1]);

        gd.addNumericField("Start Distance", distRange[0]);
        gd.addToSameRow();
        gd.addNumericField("End Distance", distRange[1]);

        gd.addNumericField("Number of Orders", orders);


        gd.addCheckbox("Flip angle?", flipAngles);
        gd.addToSameRow();
        gd.addMessage("Sometimes the angle might be calculated 180 degrees off.\nIf all orders contain the same amount of points, check this box.");

        gd.addMessage("------------------------------------------Filtering------------------------------------------------------------------------------------------------------------------------------");

        gd.addCheckbox("Remove Lone Points", toCleanup);
        gd.addToSameRow();
        gd.addNumericField("Required Neightbours", neighbours);
        gd.addToSameRow();
        gd.addNumericField("Required Distance", cleanDistance);
        gd.addMessage("Removes points if there are not at least a number of neighbours in a certain distance.\nWarning: Extremely slow for large datasets");

        gd.addMessage("------------------------------------------Visualisation------------------------------------------------------------------------------------------------------------------------------");

        gd.addNumericField("Histogram binwidth (nm)", binwidth);
        String[] colors = ownColorTable.getLuts();
        Arrays.sort(colors);

        if(runningFromIDE) defaultLUT = "NCSA PalEdit/royal.lut";
        gd.addChoice("LUT", colors, defaultLUT);

        gd.showDialog();

        if (gd.wasCanceled())
            return 0;

        filePath = gd.getNextString();

        saveSCV = gd.getNextBoolean();
        CSV_FILE_NAME = gd.getNextString();

        angRange[0] = (float) gd.getNextNumber();
        angRange[1] = (float) gd.getNextNumber();

        distRange[0] = (float) gd.getNextNumber();
        distRange[1] = (float) gd.getNextNumber();

        flipAngles = gd.getNextBoolean();

        orders = (int) gd.getNextNumber();

        toCleanup = gd.getNextBoolean();
        neighbours = (int) gd.getNextNumber();
        cleanDistance = (float) gd.getNextNumber();

        binwidth = (float) gd.getNextNumber();

        try {
            ownColorTable.setLut(colors[gd.getNextChoiceIndex()]);
        } catch (Exception e){
            logService.info("Failed to set LUT.\nTrying Default: NCSA PalEdit/6_shades.lut");
            try {
                ownColorTable.setLut("NCSA PalEdit/6_shades.lut");
            } catch (Exception e2) {
                return 0;
            }

        }

        if(filePath.equals("")){
            logService.error("No input CSV was set");
            return 0;
        }

        if(saveSCV && CSV_FILE_NAME.equals("")){
            logService.error("Set saving to CSV but no filepath was provided.");
            return 0;
        }

        totalCollumns = orders * orderCollumns;

        return 1;
    }

    @Override
    public void run() {

        ownColorTable = new OwnColorTable(lutService);

        if(setup() != 0) {

            double csvTime = System.nanoTime();

            FloatMatrix floatMatrix = null;
            List<String> collumns = new ArrayList<>();

            OwnFloatMatrixLoader ownFloatMatrix = new OwnFloatMatrixLoader();

            if(debug) {
                filePath = CSV_FILE_NAME;
                processing = false;
                try {
                    System.out.println("Found " + ownColorTable.getLuts().length + " LUTs");
                    System.out.println(Arrays.toString(ownColorTable.getLuts()));
                    ownColorTable.setLut("NCSA PalEdit/royal.lut");
                } catch (Exception e2) {
                    System.out.println("Failed to set LUT");
                    System.exit(0);
                }
            }

            try {
                floatMatrix = ownFloatMatrix.loadCSVFile(filePath);
                collumns = ownFloatMatrix.collumns;
            } catch (IOException e) {
                System.out.println("File not found.");
            } catch (LapackException e) {
                System.out.println("Lapack error");
                e.printStackTrace();
            } catch (Exception e) {
                System.out.println("Does the csv start with a header?");
            }

            assert !collumns.isEmpty();
            assert floatMatrix != null;


            int[] revOptionsIndices = new int[possible_options.length];
            int[] unitsIndices = new int[collumns.size()];


            Pattern pattern = Pattern.compile("(\\w+)( [ (\\[](\\w+)[)\\] ])?");

            for (int i = 0; i < collumns.size(); i++) {
                String header = collumns.get(i);
                Matcher matcher = pattern.matcher(header);
                if (matcher.find()) {
                    revOptionsIndices[getTheClosestMatch(possible_options, matcher.group(1))] = i;
                    unitsIndices[i] = getTheClosestMatch(unit_prefixes, matcher.group(3));
                }
            }

            csvTime = System.nanoTime() - csvTime;
            System.out.println("Loading CSV took " + String.format("%.3f", csvTime / 1000000000) + " s");
            double processingTime = System.nanoTime();

            final FloatMatrix data;
            FloatMatrix finalPossibilities;

            if(processing){
                //frame, x, y, intensity
                data = floatMatrix.getColumns(new int[]{revOptionsIndices[1], revOptionsIndices[2], revOptionsIndices[3], revOptionsIndices[5]});

                if(angRange[0] * angRange[1] * distRange[0] * distRange[1] == 0){
                    AngleAnalyzer<T> angleAnalyzer = new AngleAnalyzer<>(data, flipAngles, logService);
                    angleAnalyzer.run();


                    float[] angResult = angleAnalyzer.getAngles();
                    float[] distResult = angleAnalyzer.getDistances();

                    angRange[0] = angRange[0] == 0f ? angResult[0] : angRange[0];
                    angRange[1] = angRange[1] == 0f ? angResult[1] : angRange[1];

                    distRange[0] = distRange[0] == 0f ? distResult[0] : distRange[0];
                    distRange[1] = distRange[1] == 0f ? distResult[1] : distRange[1];


                    if(distRange[0] > distRange[1]) {
                        logService.error("The distance had to be positive: " + distRange[0] + " is larger than " + distRange[1]);
                    }
                }



                final int frames = (int) data.getColumn(0).max(); //1 indexed
                // Implement FFT to guess angle and distance

                System.out.println("Total Frames: " + frames);
                System.out.println("Total Points: " + (int) floatMatrix.getColumn(revOptionsIndices[0]).max());

                final AtomicInteger ai = new AtomicInteger(0); //Atomic Integer is a thread safe incremental integer
                final Thread[] threads = createThreadArray(); //Get maximum of threads

                final FloatMatrix[] intermediateFinals = new FloatMatrix[getNbCpus()];

                //Set the run function for each thread
                for (int ithread = 0; ithread < threads.length; ithread++) {
                    final int finalIthread = ithread;
                    threads[ithread] = new Thread(() -> {

                        intermediateFinals[finalIthread] = new FloatMatrix(0, totalCollumns);

                        for (int frame = ai.getAndIncrement(); frame <= frames; frame = ai.getAndIncrement()) {
                            final int currentAI = ai.get();
                            if (currentAI % 1000 == 0) logService.info("\r" + currentAI + "/" + frames);
                            IJ.showProgress(currentAI, frames);
                            IJ.showStatus(currentAI + "/" + frames);

                            final int[] frameIndicices = data.getColumn(0).eq(frame).findIndices();

                            final FloatMatrix frameData = data.getRows(frameIndicices);

                            final FloatMatrix subtractedX = makeSubstractedMatrix(frameData.getColumn(1));
                            final FloatMatrix subtractedY = makeSubstractedMatrix(frameData.getColumn(2));

                            final FloatMatrix distances = Distance(subtractedX, subtractedY);
                            final FloatMatrix angles = atan2(subtractedX, subtractedY, distances);


                            final int[] correctAngleAndDistance = distances.gt(distRange[0]).and(distances.lt(distRange[1]))
                                    .and(angles.gt(angRange[0])).and(angles.lt(angRange[1])).findIndices();


                            if (correctAngleAndDistance.length > 1) {

                                final FloatMatrix possibilities = new FloatMatrix(correctAngleAndDistance.length, totalCollumns);
                                FloatMatrix intermediateFinalPossibilities = new FloatMatrix(0, totalCollumns);

                                for (int i = 0; i < correctAngleAndDistance.length; i++) {
                                    int index = correctAngleAndDistance[i];
                                    possibilities.putRow(i, extend(new FloatMatrix(1, orderCollumns * 2,
                                            0,                                                  //0
                                            frame,                                                          //1
                                            1 + Math.floorDiv(index, distances.rows),                       //2
                                            frameData.get(index / distances.rows, 1),    //3
                                            frameData.get(index / distances.rows, 2),    //4
                                            frameData.get(index / distances.rows, 3),    //5
                                            1 + Math.floorMod(index, distances.rows),                       //6
                                            frameData.get(index % distances.rows, 1),    //7
                                            frameData.get(index % distances.rows, 2),    //8
                                            frameData.get(index % distances.rows, 3),    //9
                                            distances.get(index),                                           //10
                                            angles.get(index)                                               //11
                                    ), 1, totalCollumns));
                                }
                                FloatMatrix connectedIds = possibilities.getColumn(6);
                                List<Integer> toSkip = new ArrayList<>();
                                for (int i = 0; i < possibilities.rows; i++) {
                                    if (!toSkip.contains(i)) {
                                        FloatMatrix identicalIds = connectedIds.eq(connectedIds.get(i));
                                        float identicalIdSum = identicalIds.sum();
                                        if (identicalIdSum > 1.0f) {

                                            float cumX = possibilities.get(i, 3);
                                            float cumY = possibilities.get(i, 4);

                                            for (int index2 : identicalIds.findIndices()) {
                                                if (index2 != i) {
                                                    cumX += possibilities.get(index2, 3);
                                                    cumY += possibilities.get(index2, 4);
                                                    toSkip.add(index2);
                                                }
                                            }

                                            FloatMatrix temp_arr = possibilities.getRow(i).dup();
                                            temp_arr.put(3, cumX / identicalIdSum);
                                            temp_arr.put(4, cumY / identicalIdSum);

                                            intermediateFinalPossibilities = FloatMatrix.concatVertically(intermediateFinalPossibilities, temp_arr);

                                        } else {
                                            intermediateFinalPossibilities = FloatMatrix.concatVertically(intermediateFinalPossibilities, possibilities.getRow(i));
                                        }
                                    }
                                }

                                intermediateFinals[finalIthread] = FloatMatrix.concatVertically(intermediateFinals[finalIthread], connectOrders(intermediateFinalPossibilities, orders, orderCollumns));

                            } else if (correctAngleAndDistance.length == 1) {
                                int index = correctAngleAndDistance[0];
                                intermediateFinals[finalIthread] = FloatMatrix.concatVertically(intermediateFinals[finalIthread], extend(new FloatMatrix(1, orderCollumns * 2,
                                        0,
                                        frame,
                                        1 + Math.floorDiv(index, distances.rows),
                                        frameData.get(index / distances.rows, 1),
                                        frameData.get(index / distances.rows, 2),
                                        frameData.get(index / distances.rows, 3),
                                        1 + Math.floorMod(index, distances.rows),
                                        frameData.get(index % distances.rows, 1),
                                        frameData.get(index % distances.rows, 2),
                                        frameData.get(index % distances.rows, 3),
                                        distances.get(index),
                                        angles.get(index)
                                ), 1, totalCollumns));
                            }
                        }
                    }); //end of thread creation
                }

//              Start actual processing
                startAndJoin(threads);

                finalPossibilities = new FloatMatrix(0, totalCollumns);
                for (FloatMatrix floatMatrix1 : intermediateFinals) {
                    finalPossibilities = FloatMatrix.concatVertically(finalPossibilities, floatMatrix1);
                }


                //Fix ids
                int id = 0;
                for (int i = 0; i < finalPossibilities.rows; i++) {
                    finalPossibilities.put(i, 0, id++);
                }

                if (toCleanup) {
                    IJ.showStatus("Cleaning Data");
                    logService.info("Cleaning up data");
                    finalPossibilities = cleanup(finalPossibilities, neighbours, cleanDistance);
                }

                if (saveSCV) {
                    //Create Header
                    List<String> Headers = new ArrayList<>();
                    Headers.add("id");
                    Headers.add("frame");
                    String distanceUnit = (unit_prefixes[unitsIndices[revOptionsIndices[2]]].equals(unit_prefixes[unitsIndices[revOptionsIndices[3]]]) ? unit_prefixes[unitsIndices[revOptionsIndices[2]]] : ("(" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "*" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + ")^½"));
                    for (int i = 0; i < orders; i++) {
                        Headers.add(i + " index");
                        Headers.add(i + " x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]");
                        Headers.add(i + " y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]");
                        Headers.add(i + " intensity [" + unit_prefixes[unitsIndices[revOptionsIndices[5]]] + "]");
                        if (i > 0) {
                            Headers.add((i - 1) + "-" + i + " distance" + distanceUnit);
                            Headers.add((i - 1) + "-" + i + "angle");
                        }
                    }

                    SaveCSV(finalPossibilities, Headers, CSV_FILE_NAME);
                }


                processingTime = System.nanoTime() - processingTime;
                System.out.println("\nProcessing data took " + String.format("%.3f", processingTime / 1000000000) + " s");


            } else {
                finalPossibilities = floatMatrix;

            }

            assert(finalPossibilities != null);

            System.gc();

            if(finalPossibilities.rows < 10){
                logService.info("No pairs detected. Try flipping the angle or manually adjusting the distance and angles.");
            } else {
                HistogramWindow[] histograms = new HistogramWindow[orders - 1];

                for (int i = 0; i < orders - 1; i++) {
                    FloatMatrix relevantData;
                    if (i == 0) {
                        relevantData = finalPossibilities.getColumn(10);
                    } else {
                        relevantData = finalPossibilities.getColumn(10 + (i * orderCollumns));
                        relevantData = relevantData.get(relevantData.ne(0.0f).findIndices());
                    }

                    System.out.println("Found " + relevantData.rows + " connections in the " + getTitleHist(i) + " order");
                    if (relevantData.rows < 50) {
                        orders = i + 1;
                        break;
                    } //No relevant amount of data above this point
                    ImagePlus dummy = new ImagePlus("", new FloatProcessor(relevantData.toArray2()));
                    histograms[i] = new HistogramWindow(getTitleHist(i), dummy, getBins(relevantData, binwidth), distRange[0], distRange[1]);
                    if (runningFromIDE) histograms[i].getImagePlus().show();
                }

                ImagePlus Angles = new ImagePlus("", new FloatProcessor(finalPossibilities.getColumn(11).toArray2()));
                new HistogramWindow("Angles", Angles, getBins(finalPossibilities.getColumn(11), 0.005f), angRange[0], angRange[1]).getImagePlus().show();


                Plot orderPlot = new Plot("Orders", "x", "y");

                String[] colors = {"blue", "red", "green", "black"};
                String[] shapes = {"circle", "cross", "box", "diamond"};//  "line", "connected circle", "filled", "bar", "separated bar", "circle", "box", "triangle", "diamond", "cross", "x", "dot", "error bars" or "xerror bars"


                orderPlot.setColor(colors[0]); //0th
                orderPlot.add(shapes[0], toDouble(finalPossibilities.getColumn(3)), toDouble(finalPossibilities.getColumn(4)));
                for (int i = 1; i < orders; i++) {
                    orderPlot.setColor(colors[i]);
                    orderPlot.add(shapes[i], toDouble(finalPossibilities.getColumn(1 + (i * orderCollumns))), toDouble(finalPossibilities.getColumn(2 + (i * orderCollumns))));
                }

                orderPlot.setLegend(getTitlePlot(orders), Plot.AUTO_POSITION);
                orderPlot.show();


                Plot distancePlot = new Plot("Distance", "x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]", "y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]");

                for (int i = 0; i < finalPossibilities.rows; i++) {
                    distancePlot.setColor(ownColorTable.getColor(finalPossibilities.get(i, 10), distRange[0] * 1.05, distRange[1] * 0.95));

                    distancePlot.add(shapes[0], toDouble(finalPossibilities.get(i, 3)), toDouble(finalPossibilities.get(i, 4))); // 3 4
                }

                distancePlot.setLimitsToFit(true);
                addLutLegend(distancePlot, ownColorTable, "Distance", 512, distRange[0] * 1.05, distRange[1] * 0.95);
                distancePlot.show();
                distancePlot.setLimits(Float.NaN, Float.NaN, Float.NaN, Float.NaN);

            }
        }
    }

    public static void main(String[] args) {
        debug = false;
        runningFromIDE = true; //this is really dumb

        net.imagej.ImageJ ij = new ImageJ();
        ij.ui().showUI();

        ij.command().run(sSMLMA.class, true);

    }

}

