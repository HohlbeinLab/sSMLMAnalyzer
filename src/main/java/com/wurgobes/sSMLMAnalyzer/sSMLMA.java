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

// add intensity checks (higher orders should always be lower than the 0th order)
// ZOLA integration

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


    private String filePath = "";
    private String csv_target_dir = "";

    private boolean saveSCV = true;

    private final float[] angRange = new float[] {0, 0};
    private final float[] distRange = new float[] {0, 0}; //default was {1500, 2500}


    private int orders = 4; //This its the number of orders, including the 0th, so 3 would be 0th + 1st + 2nd
    private final int orderColumns = 6;
    private int totalCollumns = orders * orderColumns;


    private float binwidth = 2.5f;

    private final float[] angInput = new float[] {0, 0};
    private final float[] distInput = new float[] {0, 0};


    private boolean toCleanup = false;
    private int neighbours = 7;
    private float cleanDistance = 100;

    private boolean processing = true;

    private static boolean debug = false;

    private OwnColorTable ownColorTable;


    private boolean searchAngle = true;
    private boolean deepSearchAngle = false;
    private boolean doingRetry = false;
    private boolean flipAngles = false;
    private boolean mirrorAngles = false;

    private final boolean[][] perm = new boolean[][]{{false, false}, {false, true}, {true, false}, {true, true}};
    private final boolean[][] permReference = new boolean[][]{{false, false}, {false, true}, {true, false}, {true, true}};
    private final FloatMatrix[] angleResults = new FloatMatrix[perm.length];
    private boolean retry = false;
    private boolean foundBestResult = false;
    private boolean displayInfo = true;
    private int runNumber = 1;

    private static boolean runningFromIDE = false;

    private String defaultLUT = "Spectrum.lut";
    private final float[] lutRange = new float[]{0,0};
    private boolean visualisation = true;

    private boolean checkForIntensity = false;
    private float ratioIntensity = 1.2f;

    private final int coreCount = Prefs.getThreads();

    int[] revOptionsIndices;
    int[] unitsIndices;
    FloatMatrix floatMatrix;
    OwnFloatMatrixLoader ownFloatMatrix = new OwnFloatMatrixLoader();

    private static String debug_arg_string = "";

    public int setup() {
        if (debug || doingRetry) return 1;

        String arg;

        if (debug_arg_string.equals("")) {
            arg = Macro.getOptions();
        } else {
            arg = debug_arg_string;
        }
        if (arg != null && !arg.equals("")) {
            saveSCV = false;
            visualisation = false;
            searchAngle = false;

            String[] arguments = arg.split(" ");
            String[] keywords = {
                    "csv_in", "csv_out","angle_start", "angle_end",
                    "distance_start", "distance_end",
                    "order_number", "check_order_intensity", "check_order_ratio",
                    "angle_flip", "angle_mirror", "angle_search", "angle_deep_search",
                    "lone_pair_remove", "lone_pair_neighbours", "lone_pair_distance",
                    "visualisation", "hist_binwidth", "LUT", "LUT_start", "LUT_end"
            };
            for(String a : arguments) {
                if (a.contains("=")) {
                    if (a.contains("=")) {
                        String[] keyword_val = a.split("=");
                        try {
                            switch (keyword_val[0]) {
                                case "csv_in":
                                    filePath = keyword_val[1];
                                    break;
                                case "csv_out":
                                    saveSCV = true;
                                    csv_target_dir = keyword_val[1];
                                    break;
                                case "angle_start":
                                    angRange[0] = Float.parseFloat(keyword_val[1]);
                                    break;
                                case "angle_end":
                                    angRange[1] = Float.parseFloat(keyword_val[1]);
                                    break;
                                case "distance_start":
                                    distRange[0] = Float.parseFloat(keyword_val[1]);
                                    break;
                                case "distance_end":
                                    distRange[1] = Float.parseFloat(keyword_val[1]);
                                    break;
                                case "order_number":
                                    orders = Integer.parseInt(keyword_val[1]);
                                    break;
                                case "check_order_intensity":
                                    checkForIntensity = Boolean.parseBoolean(keyword_val[1]);
                                    break;
                                case "check_order_ratio":
                                    ratioIntensity = Float.parseFloat(keyword_val[1]);
                                    break;
                                case "angle_flip":
                                    flipAngles = Boolean.parseBoolean(keyword_val[1]);
                                    break;
                                case "angle_mirror":
                                    mirrorAngles = Boolean.parseBoolean(keyword_val[1]);
                                    break;
                                case "angle_search":
                                    searchAngle = Boolean.parseBoolean(keyword_val[1]);
                                    break;
                                case "angle_deep_search":
                                    deepSearchAngle = Boolean.parseBoolean(keyword_val[1]);
                                    break;
                                case "lone_pair_remove":
                                    toCleanup = Boolean.parseBoolean(keyword_val[1]);
                                    break;
                                case "lone_pair_neighbours":
                                    neighbours = Integer.parseInt(keyword_val[1]);
                                    break;
                                case "lone_pair_distance":
                                    cleanDistance = Float.parseFloat(keyword_val[1]);
                                    break;
                                case "visualisation":
                                    visualisation = true;
                                    break;
                                case "hist_binwidth":
                                    visualisation = true;
                                    binwidth = Float.parseFloat(keyword_val[1]);
                                    break;
                                case "LUT":
                                    visualisation = true;
                                    defaultLUT = keyword_val[1];
                                    break;
                                case "LUT_start":
                                    visualisation = true;
                                    lutRange[0] = Float.parseFloat(keyword_val[1]);
                                    break;
                                case "LUT_end":
                                    visualisation = true;
                                    lutRange[1] = Float.parseFloat(keyword_val[1]);
                                    break;
                                default:
                                    logService.error("Keyword " + keyword_val[0] + " not found\nDid you mean: " + getTheClosestMatch(keywords, keyword_val[0]) + "?");
                                    return 0;
                            }
                        } catch (ArrayIndexOutOfBoundsException e) {
                            logService.error("Malformed token: " + a + ".\nDid you remember to format it as keyword=value?");
                            return 0;
                        } catch (Exception e){
                            logService.error("Failed to parse argument:" + keyword_val[1]);
                            return 0;
                        }
                    } else {
                        logService.error("Malformed token: " + a + ".\nDid you remember to format it as keyword=value?");
                        return 0;
                    }

                    try {
                        ownColorTable.setLut(defaultLUT);
                    } catch (Exception e) {
                        logService.info("Failed to set LUT.\nTrying Default: NCSA PalEdit/6_shades.lut");
                        try {
                            ownColorTable.setLut("NCSA PalEdit/6_shades.lut");
                        } catch (Exception e2) {
                            logService.info("Failed to set LUT again?.\n");
                            return 0;
                        }

                    }
                }
            }
        } else {
            GenericDialogPlus gd = new GenericDialogPlus("settings");

            gd.addFileField("CSV input", filePath, 25);

            gd.addCheckbox("Save to CSV?", saveSCV);
            gd.addToSameRow();
            gd.addDirectoryField("CSV output directory", csv_target_dir, 25);

            gd.addMessage("------------------------------------------Angles and Distances------------------------------------------------------------------------------------------------------------------------------");

            gd.addMessage("If you want to override the calculated values change the values below.\nWill only overwrite values if they are NOT 0.");


            gd.addNumericField("Start Angle", angRange[0]);
            gd.addToSameRow();
            gd.addNumericField("End Angle", angRange[1]);

            gd.addNumericField("Start Distance", distRange[0]);
            gd.addToSameRow();
            gd.addNumericField("End Distance", distRange[1]);

            gd.addNumericField("Number of Orders", orders);

            gd.addCheckbox("Intensity Order Required?", checkForIntensity);
            gd.addToSameRow();
            gd.addNumericField("Ratio between orders:", ratioIntensity);
            gd.addToSameRow();
            gd.addMessage("Require that each order has at least (ratio) intensity more than its next order");

            gd.addCheckbox("Flip angle?", flipAngles);
            gd.addToSameRow();
            gd.addCheckbox("Mirror angle?", mirrorAngles);
            gd.addToSameRow();
            gd.addMessage("Sometimes the angle might be calculated 180 degrees off, or mirrored.\n Change these boxes if you get very few or no points.");
            gd.addCheckbox("Search for angles?", searchAngle);
            gd.addToSameRow();
            gd.addMessage("Sometimes the angle is not correctly calculated the first time. With this setting the plugin will search for the correct angle and tell you what settings it used.");
            gd.addCheckbox("Search for the angle with the most pairs?", deepSearchAngle);
            gd.addToSameRow();
            gd.addMessage("Try and find the permutation of above options that results in the most pairs.");

            gd.addMessage("------------------------------------------Filtering------------------------------------------------------------------------------------------------------------------------------");

            gd.addCheckbox("Remove Lone Points", toCleanup);
            gd.addToSameRow();
            gd.addNumericField("Required Neightbours", neighbours);
            gd.addToSameRow();
            gd.addNumericField("Required Distance", cleanDistance);
            gd.addMessage("Removes points if there are not at least a number of neighbours in a certain distance.\nWarning: Extremely slow for large datasets");

            gd.addMessage("------------------------------------------Visualisation------------------------------------------------------------------------------------------------------------------------------");

            gd.addCheckbox("Visualise results", visualisation);

            gd.addNumericField("Histogram binwidth", binwidth);
            String[] colors = ownColorTable.getLuts();
            Arrays.sort(colors);


            if (runningFromIDE) defaultLUT = "NCSA PalEdit/royal.lut";
            gd.addChoice("LUT", colors, defaultLUT);

            gd.addMessage("Modifies the LUT range of the distances. Will use the calculated distance ranges if not set.");
            gd.addNumericField("Start LUT", lutRange[0]);
            gd.addToSameRow();
            gd.addNumericField("End LUT", lutRange[1]);

            gd.showDialog();

            if (gd.wasCanceled())
                return 0;

            filePath = gd.getNextString();

            saveSCV = gd.getNextBoolean();
            csv_target_dir = gd.getNextString();

            angInput[0] = (float) gd.getNextNumber();
            angInput[1] = (float) gd.getNextNumber();

            distInput[0] = (float) gd.getNextNumber();
            distInput[1] = (float) gd.getNextNumber();

            orders = (int) gd.getNextNumber();

            checkForIntensity = gd.getNextBoolean();
            ratioIntensity = (float) gd.getNextNumber();

            flipAngles = gd.getNextBoolean();
            mirrorAngles = gd.getNextBoolean();

            searchAngle = gd.getNextBoolean();
            deepSearchAngle = gd.getNextBoolean();

            toCleanup = gd.getNextBoolean();
            neighbours = (int) gd.getNextNumber();
            cleanDistance = (float) gd.getNextNumber();

            visualisation = gd.getNextBoolean();
            binwidth = (float) gd.getNextNumber();

            try {
                defaultLUT = colors[gd.getNextChoiceIndex()];
                ownColorTable.setLut(defaultLUT);
            } catch (Exception e) {
                logService.info("Failed to set LUT.\nTrying Default: NCSA PalEdit/6_shades.lut");
                try {
                    ownColorTable.setLut("NCSA PalEdit/6_shades.lut");
                } catch (Exception e2) {
                    return 0;
                }

            }

            lutRange[0] = (float) gd.getNextNumber();
            lutRange[1] = (float) gd.getNextNumber();
        }

        if(deepSearchAngle) searchAngle = true;

        if(filePath.equals("")){
            logService.error("No input CSV was set");
            return 0;
        }

        if(saveSCV && csv_target_dir.equals("")){
            logService.error("Set saving to CSV but no filepath was provided.");
            return 0;
        }

        if (!(visualisation || saveSCV)) {
            logService.error("No output method of any sorts is selected.\nSelect either Visualisation or Save to SCV.");
            return 0;
        }

        totalCollumns = orders * orderColumns;

        return 1;
    }

    @Override
    public void run() {

        ownColorTable = new OwnColorTable(lutService);
        if(doingRetry) {
            try {
                ownColorTable.setLut(defaultLUT);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if(setup() != 0) {
            double csvTime = System.nanoTime();

            List<String> collumns = new ArrayList<>();



            if (debug) {
                filePath = csv_target_dir + "\\all_orders.csv";
                processing = false;
                try {
                    System.out.println("Found " + ownColorTable.getLuts().length + " LUTs");
                    System.out.println(Arrays.toString(ownColorTable.getLuts()));
                    ownColorTable.setLut("NCSA PalEdit/royal.lut");
                } catch (Exception e2) {
                    System.out.println("Failed to set LUT");
                    System.exit(0);
                }
                //flipAngles = true;
                //mirrorAngles = true;
            }

            if(!doingRetry){
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

                revOptionsIndices = new int[possible_options.length];
                unitsIndices = new int[collumns.size()];

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

            }

            if(searchAngle){
                System.out.println("Run: " + runNumber + ". Determining Angle with settings: Flip Angle: " + flipAngles  + ", Mirror Angle: " + mirrorAngles);
            }


            double processingTime = System.nanoTime();
            final FloatMatrix data;
            FloatMatrix finalPossibilities;


            if(processing){
                //frame, x, y, intensity
                data = floatMatrix.getColumns(new int[]{revOptionsIndices[1], revOptionsIndices[2], revOptionsIndices[3], revOptionsIndices[5]});

                if(angRange[0] * angRange[1] * distRange[0] * distRange[1] == 0){
                    AngleAnalyzer<T> angleAnalyzer = new AngleAnalyzer<>(data, flipAngles, mirrorAngles,logService, debug);
                    angleAnalyzer.run();


                    float[] angResult = angleAnalyzer.getAngles();
                    float[] distResult = angleAnalyzer.getDistances();

                    angRange[0] = angInput[0] == 0f ? angResult[0] : angInput[0];
                    angRange[1] = angInput[1] == 0f ? angResult[1] : angInput[1];

                    distRange[0] = distInput[0] == 0f ? distResult[0] : distInput[0];
                    distRange[1] = distInput[1] == 0f ? distResult[1] : distInput[1];

                    lutRange[0] = distRange[0] == 0f ? distResult[0] : lutRange[0];
                    lutRange[1] = distRange[1] == 0f ? distResult[1] : lutRange[1];

                    if(distRange[0] > distRange[1]) {
                        logService.error("The distance had to be positive: " + distRange[0] + " is larger than " + distRange[1]);
                    }
                }



                final int frames = (int) data.getColumn(0).max(); //1 indexed
                // Implement FFT to guess angle and distance

                System.out.println("Total Frames: " + frames);
                System.out.println("Total Points: " + (int) floatMatrix.getColumn(revOptionsIndices[0]).max());

                final AtomicInteger ai = new AtomicInteger(0); //Atomic Integer is a thread safe incremental integer
                final Thread[] threads = createThreadArray(coreCount); //Get maximum of threads

                final FloatMatrix[] intermediateFinals = new FloatMatrix[coreCount];

                //Set the run function for each thread
                for (int ithread = 0; ithread < threads.length; ithread++) {
                    final int finalIthread = ithread;
                    final boolean finalIntensityCheck = checkForIntensity;
                    final float finalRatioIntensity = ratioIntensity;
                    threads[ithread] = new Thread(() -> {

                        intermediateFinals[finalIthread] = new FloatMatrix(0, totalCollumns);

                        for (int frame = ai.getAndIncrement(); frame <= frames; frame = ai.getAndIncrement()) {
                            if (frame % 1000 == 0) logService.info("\r" + frame + "/" + frames);
                            IJ.showProgress(frame, frames);
                            IJ.showStatus(frame + "/" + frames);

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
                                    if(finalIntensityCheck &&
                                            (frameData.get(index / distances.rows, 3)/
                                            frameData.get(index % distances.rows, 3)) > finalRatioIntensity
                                            ) continue;
                                    possibilities.putRow(i, extend(new FloatMatrix(1, orderColumns * 2,
                                            0,                                                  //0
                                            frame,                                                 //1
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

                                intermediateFinals[finalIthread] = FloatMatrix.concatVertically(intermediateFinals[finalIthread], connectOrders(intermediateFinalPossibilities, orders, orderColumns));

                            } else if (correctAngleAndDistance.length == 1) {
                                int index = correctAngleAndDistance[0];
                                intermediateFinals[finalIthread] = FloatMatrix.concatVertically(intermediateFinals[finalIthread], extend(new FloatMatrix(1, orderColumns * 2,
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

                processingTime = System.nanoTime() - processingTime;
                System.out.println("\nProcessing data took " + String.format("%.3f", processingTime / 1000000000) + " s");


            } else {
                data = floatMatrix;
                finalPossibilities = floatMatrix;

            }

            assert(finalPossibilities != null);

            System.gc();


            System.out.println("Pairs in the 0th-1st order found: " + finalPossibilities.rows);



            if(finalPossibilities.rows < 10) {
                logService.info("No pairs detected. Try flipping the angle or manually adjusting the distance and angles.");
                retry = true;
            }
            if(finalPossibilities.rows < 0.1 * data.rows){
                if(searchAngle)
                    logService.info("Very few pairs were detected.");
                else
                    logService.info("Very few pairs were detected. A search for a better angle could improve results.");
                retry = true;
            }

            boolean[][] checks = checkForRetry(finalPossibilities.getColumn(11));

            if(sum(checks[0]) < 2){
                logService.info("Nothing resembling a guassian was found for the angles.");
                retry = true;
            }

            if(checks[1][0]) {
                if(checks[1][1]) logService.info("The left tail of the angle histogram seems to be partially cut-off");
                else logService.info("The right tail of the angle histogram seems to be partially cut-off");
            }

            if((retry && searchAngle) | deepSearchAngle) {
                boolean[] curr_perm = new boolean[]{flipAngles, mirrorAngles};

                for(int i = 0; i < perm.length ; i++){
                    if(Arrays.equals(perm[i], curr_perm)){
                        angleResults[i] = new FloatMatrix();
                        angleResults[i].copy(finalPossibilities);
                        perm[i] = null;

                        if(perm[(i+1)%perm.length] != null){
                            flipAngles = perm[(i+1)%perm.length][0];
                            mirrorAngles = perm[(i+1)%perm.length][1];
                        }
                        break;
                    }
                }




                boolean all_filled = true;
                for(FloatMatrix result : angleResults){
                    if (result == null) {
                        all_filled = false;
                        break;
                    }
                }
                if(!all_filled){
                    doingRetry = true;

                    angRange[0] = angInput[0] == 0f ? angInput[0] : angRange[0];
                    angRange[1] = angInput[1] == 0f ? angInput[1] : angRange[1];

                    distRange[0] = distInput[0] == 0f ? distInput[0] : distRange[0];
                    distRange[1] = distInput[1] == 0f ? distInput[1] : distRange[1];

                    runNumber++;
                    run();

                } else {
                    int max = angleResults[0].rows;
                    int ind = 0;
                    for(int i = 1; i < angleResults.length; i++){
                        if(angleResults[i].rows > max) {
                            max = angleResults[i].rows;
                            ind = i;
                        }

                    }
                    finalPossibilities = angleResults[ind];

                    angRange[0] = finalPossibilities.getColumn(11).min();
                    angRange[1] = finalPossibilities.getColumn(11).max();

                    distRange[0] = finalPossibilities.getColumn(10).min();
                    distRange[1] = finalPossibilities.getColumn(10).max();

                    foundBestResult = true;

                    String message = "\nBest result was " + max + " pairs found.\n" +
                        "The following values were found:\n" +
                        "\tAngle(rad): " +  angRange[0] + " to " + angRange[1] + "\n" +
                        "\tDistance: " +  distRange[0] + " to " + distRange[1]+ "\n" +
                        "With the settings:\n" +
                        "\tFlip Angle: " + permReference[ind][0] + "\n" +
                        "\tMirror Angle: " + permReference[ind][1] + "\n";
                    if(saveSCV) message += ("CSV files were saved to the folder: " + csv_target_dir + "\n");

                    if(runningFromIDE) System.out.println(message);
                    IJ.showMessage(message);
                }

            }
            if((!(retry && searchAngle) | (foundBestResult && deepSearchAngle)) && displayInfo) {
                displayInfo = false;

                if (processing && toCleanup) {
                    IJ.showStatus("Cleaning Data");
                    logService.info("Cleaning up data");

                    try {
                        finalPossibilities = cleanup(finalPossibilities, neighbours, cleanDistance, coreCount);
                    } catch(Exception e){
                        logService.info("Checking for neighbours failed. A retry could work.");
                    }
                }

                //id, frame, x, y, intensity, distance
                final FloatMatrix halfOrderMatrix = new FloatMatrix(finalPossibilities.rows, orderColumns);
                final FloatMatrix allOrdersCombined = new FloatMatrix(finalPossibilities.rows, orderColumns);

                halfOrderMatrix.putColumn(0, finalPossibilities.getColumn(0)); //id
                halfOrderMatrix.putColumn(1, finalPossibilities.getColumn(1)); //frame
                halfOrderMatrix.putColumn(4, finalPossibilities.getColumn(5)); //intensity
                halfOrderMatrix.putColumn(5, finalPossibilities.getColumn(10)); //distance

                allOrdersCombined.putColumn(0, finalPossibilities.getColumn(0)); //id
                allOrdersCombined.putColumn(1, finalPossibilities.getColumn(1)); //frame
                allOrdersCombined.putColumn(4, finalPossibilities.getColumn(5)); //intensity
                allOrdersCombined.putColumn(5, finalPossibilities.getColumn(10)); //distance

                halfOrderMatrix.putColumn(2, finalPossibilities.getColumn(3).add(finalPossibilities.getColumn(7)).divi(2.0f)); //x
                halfOrderMatrix.putColumn(3, finalPossibilities.getColumn(4).add(finalPossibilities.getColumn(8)).divi(2.0f)); //y

                FloatMatrix offsets = new FloatMatrix(finalPossibilities.rows, 2);
                offsets.putColumn(0, finalPossibilities.getColumn(7).sub(finalPossibilities.getColumn(3)).divi(2.0f));
                offsets.putColumn(1, finalPossibilities.getColumn(8).sub(finalPossibilities.getColumn(4)).divi(2.0f));


                for(int i = 0; i < orders - 2; i++){
                    FloatMatrix relevantRows = finalPossibilities.getColumn(12 + (i * orderColumns)).ne(0.0f);
                    offsets.getColumn(0).addi(finalPossibilities.getColumn(13 + (i * orderColumns)).subi(finalPossibilities.getColumn(7 + (i * orderColumns))).divi(2.0f).muli(relevantRows));
                    offsets.getColumn(1).addi(finalPossibilities.getColumn(14 + (i * orderColumns)).subi(finalPossibilities.getColumn(8 + (i * orderColumns))).divi(2.0f).muli(relevantRows));
                }


                allOrdersCombined.putColumn(2, finalPossibilities.getColumn(3));
                allOrdersCombined.putColumn(3, finalPossibilities.getColumn(4));

                FloatMatrix relevantRows = finalPossibilities.getColumn((orders-1) * orderColumns).ne(0.0f);
                for(int i = orders; i > 1; i--){
                    allOrdersCombined.getColumn(2).addi(offsets.getColumn(0).divi((float) i).muli(relevantRows));
                    allOrdersCombined.getColumn(3).addi(offsets.getColumn(1).divi((float) i).muli(relevantRows));

                    relevantRows.xori(finalPossibilities.getColumn((i-2) * orderColumns).ne(0.0f));
                }


                if(visualisation) {
                    HistogramWindow[] histograms = new HistogramWindow[orders - 1];

                    for (int i = 0; i < orders - 1; i++) {
                        FloatMatrix relevantData;
                        if (i == 0) {
                            relevantData = finalPossibilities.getColumn(10);
                        } else {
                            relevantData = finalPossibilities.getColumn(10 + (i * orderColumns));
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
                    /////////////////
                    ImagePlus Angles = new ImagePlus("", new FloatProcessor(finalPossibilities.getColumn(11).toArray2()));
                    HistogramWindow angleHist = new HistogramWindow("Angles", Angles, getBins(finalPossibilities.getColumn(11), 0.005f), angRange[0], angRange[1]);
                    if(runningFromIDE) angleHist.getImagePlus().show();

                    ///////////////////////////////////////////////////////////// PLOTS
                    String[] colors = {"blue", "red", "green", "black"};
                    String[] shapes = {"cross", "circle", "box", "diamond"};//  "line", "connected circle", "filled", "bar", "separated bar", "circle", "box", "triangle", "diamond", "cross", "x", "dot", "error bars" or "xerror bars"


                    //////////////////////////////////////////////
                    Plot orderPlot = new Plot("Orders", "x", "y");

                    orderPlot.setColor(colors[0]); //0th
                    orderPlot.add(shapes[0], toDouble(finalPossibilities.getColumn(3)), toDouble(finalPossibilities.getColumn(4)));
                    for (int i = 1; i < orders; i++) {
                        orderPlot.setColor(colors[i]);
                        orderPlot.add(shapes[i], toDouble(finalPossibilities.getColumn(1 + (i * orderColumns))), toDouble(finalPossibilities.getColumn(2 + (i * orderColumns))));
                    }

                    orderPlot.setLegend(getTitlePlot(orders), Plot.AUTO_POSITION);
                    orderPlot.show();
                    //////////////////////////////////////////////

                    Plot distancePlot = new Plot("Distance", "x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]", "y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]");

                    for (int i = 0; i < finalPossibilities.rows; i++) {
                        distancePlot.setColor(ownColorTable.getColor(allOrdersCombined.get(i, 5), distRange[0], distRange[1]));

                        //distancePlot.add(shapes[0], toDouble(finalPossibilities.get(i, 3)), toDouble(finalPossibilities.get(i, 4))); // 3 4
                        distancePlot.add(shapes[0], toDouble(halfOrderMatrix.get(i, 2)), toDouble(halfOrderMatrix.get(i, 3))); // 3 4
                    }

                    distancePlot.setLimitsToFit(true);
                    addLutLegend(distancePlot, ownColorTable, "Distance", 512, distRange[0], distRange[1]);
                    distancePlot.show();
                    distancePlot.setLimits(Float.NaN, Float.NaN, Float.NaN, Float.NaN);

                }

                if (saveSCV) {
                    //Create Header
                    List<String> LongHeader = new ArrayList<>();
                    List<String> ShortHeader = new ArrayList<>();
                    LongHeader.add("id");
                    ShortHeader.add("id");
                    LongHeader.add("frame");
                    ShortHeader.add("frame");

                    String distanceUnit = (unit_prefixes[unitsIndices[revOptionsIndices[2]]].equals(unit_prefixes[unitsIndices[revOptionsIndices[3]]]) ? unit_prefixes[unitsIndices[revOptionsIndices[2]]] : ("(" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "*" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + ")^½"));
                    for (int i = 0; i < orders; i++) {
                        LongHeader.add("index " + i);
                        LongHeader.add("x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "] " + i);
                        LongHeader.add("y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "] " + i);
                        LongHeader.add("intensity [" + unit_prefixes[unitsIndices[revOptionsIndices[5]]] + "] " + i);
                        if (i > 0) {
                            LongHeader.add((i - 1) + "-" + i + " distance [" + distanceUnit + "]");
                            LongHeader.add((i - 1) + "-" + i + "angle");
                        }
                    }
                    ShortHeader.add("x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]");
                    ShortHeader.add("y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]");
                    ShortHeader.add("intensity [" + unit_prefixes[unitsIndices[revOptionsIndices[5]]] + "]");
                    ShortHeader.add("z [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]");

                    SaveCSV(finalPossibilities, LongHeader, csv_target_dir + "\\all_orders.csv");
                    SaveCSV(halfOrderMatrix, ShortHeader, csv_target_dir + "\\accurate_positions.csv");
                    SaveCSV(finalPossibilities.getColumns(new int[]{0, 1, 3, 4, 5, 10}), ShortHeader, csv_target_dir + "\\thunderSTORM.csv");
                }
            }
        }
    }

    public static void main(String[] args) {
        debug = false;
        runningFromIDE = true; //this is really dumb

        /*
        "csv_in", "csv_out","angle_start", "angle_end",
        "distance_start", "distance_end",
        "order_number", "check_order_intensity", "check_order_ratio",
        "angle_flip", "angle_mirror", "angle_search", "angle_deep_search",
        "lone_pair_remove", "lone_pair_neighbours", "lone_pair_distance",
        "visualisation", "hist_binwidth", "LUT", "LUT_start", "LUT_end"
        */

        //private String csv_target_dir = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\results";

        //private String filePath = "F:\\ThesisData\\output\\output3_drift.csv";
        //private String filePath = "F:\\ThesisData\\output\\combined_drift.csv";
        //private float[] angRange = new float[] {(float) (-0.03 * Math.PI), (float) (0.07 * Math.PI)};
        //private float[] distRange = new float[] {1500, 2200};

        //private String filePath = "F:\\ThesisData\\output\\4_grating_drift.csv";
        //private final float[] angRange = {(float) (-1 * Math.PI), (float) (-0.95 * Math.PI) }; //more than and less than
        //private final float[] distRange = {1940, 2600}; //1800 3000 (1940, 2240)

        debug_arg_string = "csv_in=F:\\ThesisData\\output\\combined_drift.csv angle_start=-0.094 angle_end=0.22 distance_start=1500 distance_end=2200 visualisation=true";

        net.imagej.ImageJ ij = new ImageJ();
        ij.ui().showUI();

        ij.command().run(sSMLMA.class, true);

    }

}

