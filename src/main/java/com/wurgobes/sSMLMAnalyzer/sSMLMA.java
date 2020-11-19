package com.wurgobes.sSMLMAnalyzer;



import ij.*;

import ij.gui.Plot;

import ij.gui.HistogramWindow;


import ij.process.FloatProcessor;


import net.imagej.ImageJ;

import net.imagej.lut.LUTService;
import net.imglib2.*;

import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;


import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

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


@Plugin(type = Command.class, menuPath = "Plugins>Spectral Analyzer>Analyze Pairs")
public class sSMLMA implements Command {


    @Parameter
    private LogService logService;

    @Parameter
    private LUTService lutService;


    private final int[] unit_decades = {0, 0, 0, -2, -3, -6, -9, -10, -12, -15};

    private final String[] possible_options = {"id", "frame", "x", "y", "z", "intensity", "offset", "bkgstd", "sigma1", "sigma2", "uncertainty", "detections", "chi"};
    private final String[] unit_prefixes = {"null", "photons", "m", "cm", "mm", "um", "nm", "ang", "pm", "fm"};


    //private String filePath = "F:\\ThesisData\\output\\output3_drift.csv";
    //private String filePath = "F:\\ThesisData\\output\\combined_drift.csv";
    private String filePath = "F:\\ThesisData\\output\\4_grating_drift.csv";


    private String CSV_FILE_NAME = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\grating_cleaned.csv";
    //private String CSV_FILE_NAME = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\combined_output.csv";


    private boolean saveSCV = false;

    //@Parameter
    //private float[] angRange = new float[] {(float) (-0.04 * Math.PI), (float)(0.04 * Math.PI)};
    private final float[] angRange = {(float) (-1 * Math.PI), (float) (-0.95 * Math.PI) }; //more than and less than

    //@Parameter
    //private float[] distRange = new float[] {1500, 2500}; //default was {2500, 4500}
    private final float[] distRange = {1940, 2600}; //1800 3000 (1940, 2240)


    private int orders = 4; //This its the number of orders, including the 0th, so 3 would be 0th + 1st + 2nd
    private final int orderCollumns = 6;
    private int totalCollumns = orders * orderCollumns;


    private float binwidth = 2.5f;


    private boolean toCleanup = false;

    private boolean processing = true;

    private static boolean debug = false;

    private OwnColorTable ownColorTable;


    public int setup(){
        if (debug) return 1;

        GenericDialogPlus gd = new GenericDialogPlus("settings");

        gd.addFileField("CSV input", filePath, 50);

        gd.addCheckbox("Save to CSV?", saveSCV);
        gd.addFileField("CSV output", CSV_FILE_NAME, 50);

        gd.addNumericField("Start Angle", angRange[0]);
        gd.addNumericField("End Angle", angRange[1]);

        gd.addNumericField("Start Distance", distRange[0]);
        gd.addNumericField("End Distance", distRange[1]);

        gd.addNumericField("Number of Orders", orders);

        gd.addCheckbox("Remove Lone Points", toCleanup);
        gd.addNumericField("Histogram binwidth (nm)", binwidth);
        String[] colors = ownColorTable.getLuts();
        Arrays.sort(colors);
        gd.addChoice("LUT", colors, "NCSA PalEdit/6_shades.lut");

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

        orders = (int) gd.getNextNumber();

        toCleanup = gd.getNextBoolean();

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

            OwnFloatMatrix ownFloatMatrix = new OwnFloatMatrix();

            if(debug) {
                filePath = CSV_FILE_NAME;
                processing = false;
                try {
                    System.out.println("Found " + ownColorTable.getLuts().length + " LUTs");
                    System.out.println(Arrays.toString(ownColorTable.getLuts()));
                    ownColorTable.setLut("NCSA PalEdit/6_shades.lut");
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


            FloatMatrix finalPossibilities;
            if (processing) {
                //frame, x, y
                final FloatMatrix data = floatMatrix.getColumns(new int[]{revOptionsIndices[1], revOptionsIndices[2], revOptionsIndices[3], revOptionsIndices[5]});


                final int frames = (int) data.getColumn(0).max(); //1 indexed

                csvTime = System.nanoTime() - csvTime;
                System.out.println("Loading CSV took " + String.format("%.3f", csvTime / 1000000000) + " s");
                double processingTime = System.nanoTime();

                int id = 1;

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

                        for (int frame = ai.getAndIncrement(); frame <= frames;frame = ai.getAndIncrement()) {
                            final int currentAI = ai.get();
                            if(currentAI%1000 == 0)System.out.println(currentAI + "/" + frames);
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
                                            1 + (index % distances.rows),                                   //6
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
                                        1 + (index % distances.rows),
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
                for(FloatMatrix floatMatrix1 : intermediateFinals){
                    finalPossibilities = FloatMatrix.concatVertically(finalPossibilities, floatMatrix1);
                }




                //Fix ids
                id = 0;
                List<String> Headers = new ArrayList<>();
                Headers.add("id");
                Headers.add("frame");
                String distanceUnit = (unit_prefixes[unitsIndices[revOptionsIndices[2]]].equals(unit_prefixes[unitsIndices[revOptionsIndices[3]]]) ? unit_prefixes[unitsIndices[revOptionsIndices[2]]] : ("(" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "*" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + ")^½"));
                for (int i = 0; i < orders; i++) {
                    Headers.add(i + " index");
                    Headers.add(i + " x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]");
                    Headers.add(i + " y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]");
                    if (i > 0) {
                        Headers.add((i - 1) + "-" + i + " distance" + distanceUnit);
                        Headers.add((i - 1) + "-" + i + "angle");
                    }
                }

                for (int i = 0; i < finalPossibilities.rows; i++) {
                    finalPossibilities.put(i, 0, id++);
                }

                if (toCleanup) {
                    IJ.showStatus("Cleaning Data");
                    finalPossibilities = cleanup(finalPossibilities, 7, 100);
                }

                if (saveSCV) SaveCSV(finalPossibilities, Headers, CSV_FILE_NAME);


                processingTime = System.nanoTime() - processingTime;
                System.out.println("\nProcessing data took " + String.format("%.3f", processingTime / 1000000000) + " s");

            } else {
                finalPossibilities = floatMatrix;
            }
            System.gc();

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
            }




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







            Plot distancePlot = new Plot("Distance", "x", "y");

            for (int i = 0; i < finalPossibilities.rows; i++) {
                float distance = finalPossibilities.get(i, 10);


                //if(distance < 2125 && distance > 2050) distancePlot.setColor("green");
                //else if(distance < 2425 && distance > 2300) distancePlot.setColor("red");
                //else continue;

                distancePlot.setColor(ownColorTable.getColor(finalPossibilities.get(i, 10), distRange[0] * 1.05, distRange[1] * 0.95));
                //distancePlot.setColor(ownColorTable.getColor(finalPossibilities.get(i, 10), 1800f, 2100f));
                distancePlot.add(shapes[0], toDouble(finalPossibilities.get(i, 3)), toDouble(finalPossibilities.get(i, 4)));
            }

            distancePlot.show();
            distancePlot.setLimits(Float.NaN, Float.NaN, Float.NaN, Float.NaN);




            //x: 35000 to 39000
            //y: 24000 to 30000

            System.gc();

            boolean visualisation = true;

            if (visualisation) {
                long total_resolution = (long) (finalPossibilities.getColumn(3).max() * finalPossibilities.getColumn(4).max() * ((finalPossibilities.getColumn(10).max() - finalPossibilities.getColumn(10).min()) + 1));
                //Arrayimg max size = 2^31-1

                float reduce = (float) Math.sqrt(total_resolution / 2147483647f) + 15;

                Interval interval = new FinalInterval(
                        (long) (finalPossibilities.getColumn(3).max() / reduce),
                        (long) (finalPossibilities.getColumn(4).max() / reduce),
                        (long) (finalPossibilities.getColumn(10).max() - finalPossibilities.getColumn(10).min()) + 1);

                System.out.println(interval.toString());
                RealPointSampleList<FloatType> realPointSampleList = createPointList(finalPossibilities, new int[]{3, 4, 10}, 5, reduce, finalPossibilities.getColumn(10).min());

                final ImgFactory<FloatType> imgFactory = new ArrayImgFactory<>(new FloatType());
                final Img<FloatType> img1 = imgFactory.create(interval);

                RandomAccess<FloatType> randomAccess = img1.randomAccess();

                RealCursor<FloatType> cursor = realPointSampleList.cursor();
                cursor.next();

                while (cursor.hasNext()) {
                    float[] pos = new float[3];
                    cursor.localize(pos);
                    randomAccess.setPosition(toLong(pos));
                    randomAccess.get().set(cursor.get());
                    cursor.next();
                }
                ImagePlus imagePlus = ImageJFunctions.show(img1);
                IJ.run("Hyperstack to Stack");
            }


            //------------------------------------------------------------------------------
        /*
        for(int i = 0; i < orders - 1; i++){
            HistogramWindow hist = histograms[i];
            for(int maximum : getMaxima(hist.getHistogram(), 5)) System.out.println(i + ": " + hist.getXValues()[maximum]);
        }


        int[] firstIndices = finalPossibilities.getColumn(8).lt(1833).findIndices();
        int[] secondIndices = finalPossibilities.getColumn(8).ge(1833).findIndices();

        FloatMatrix first = finalPossibilities.getRows(firstIndices).getColumn(13);
        FloatMatrix second = finalPossibilities.getRows(secondIndices).getColumn(13);
        float[][] firstValues = first.toArray2();
        float[][] secondValues = second.toArray2();

        new HistogramWindow("dummy 1", new ImagePlus("", new FloatProcessor(firstValues)), (int) ((first.max()-first.min())/20), distRange[0], distRange[1]);
        new HistogramWindow("dummy 2", new ImagePlus("", new FloatProcessor(secondValues)), (int) ((second.max()-second.min())/20), distRange[0], distRange[1]);
         */

            //CurveFitter curveFitter = new CurveFitter(hist.getXValues(), toDouble(hist.getHistogram()));
            //curveFitter.doFit(CurveFitter.GAUSSIAN);

            //System.out.println(curveFitter.getStatusString());
            //System.out.println(Arrays.toString(curveFitter.getParams()));


        }
    }

    public static void main(String[] args) {
        debug = false;

        net.imagej.ImageJ ij = new ImageJ();
        ij.ui().showUI();

        ij.command().run(sSMLMA.class, true);

    }

}

