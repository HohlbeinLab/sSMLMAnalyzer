package com.wurgobes.sSMLMAnalyzer;

/*
Spectral Super Resolution Pair Finder
(c) 2021 Martijn Gobes, Wageningen University.
Based on the work by Koen Martens

This plugin works by either going to the menu path: Plugins>Spectral Analyzer>Analyze Pairs
or by using a macro (example provided in readme.md

The class is instantiated by ImageJ and run(String arg) is called.
If any angle/distance value is not set, AngleAnalyzer is called on the csv data to calculate these.

Pairs are created by filtering the options based on the angle and distance.
If multiple options are found, the best one is selected.

Pairs are then chained together to create multi-order connections (1-2, 2-3, 3-4 -> 1-2-3-4)

This software is released under the GPL v3. You may copy, distribute and modify
the software as long as you track changes/dates in source files. Any
modifications to or software including (via compiler) GPL-licensed code
must also be made available under the GPL along with build & install instructions.
https://www.gnu.org/licenses/gpl-3.0.en.html

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */


import ij.*;
import static ij.util.ThreadUtil.*;
import ij.gui.HistogramWindow;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import net.imagej.ImageJ;
import net.imagej.lut.LUTService;

import net.imglib2.type.numeric.IntegerType;

import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.*;

import java.io.*;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.charset.StandardCharsets;
import java.security.CodeSource;
import java.security.ProtectionDomain;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import fiji.util.gui.GenericDialogPlus;

import org.jblas.FloatMatrix;
import org.jblas.exceptions.LapackException;

import static com.wurgobes.sSMLMAnalyzer.Util.*;
import static com.wurgobes.sSMLMAnalyzer.levenshtein.getTheClosestMatch;


@Plugin(type = Command.class, menuPath = "Plugins>Spectral Analyzer>Analyze Pairs")

// TODO
// custom backup LUT incase imagej cant provide one
// refactoring:
//   FloatMatrix to DoubleMatrix
//   No concating in calculations
// test angle settings 0 to 360?
// include ZOLA rendering
// care about z?

// T is only used when calling AngleAnalyzer and denotes the type of image created
public class sSMLMA <T extends IntegerType<T>> implements Command {

    //The services are passed through from ImageJ automatically
    @Parameter
    private LogService logService;

    @Parameter
    private LUTService lutService;

    // Variables related to loading the csv (which collumn is what variable)
    private final String[] possible_options = {"id", "frame", "x", "y", "z", "intensity", "offset", "bkgstd", "sigma1", "sigma2", "uncertainty", "detections", "chi"};
    private final String[] unit_prefixes = {"null", "photons", "m", "cm", "mm", "um", "nm", "ang", "pm", "fm"};
    int[] revOptionsIndices;
    int[] unitsIndices;

    // Variable in which the csv is loaded into (and the class used for loading it)
    OwnFloatMatrixLoader ownFloatMatrixLoader = new OwnFloatMatrixLoader();
    FloatMatrix floatMatrix;

    //Filepaths for the input csv and result csv directory
    private String filePath = "";
    private String csv_target_dir = "";

    private boolean saveSCV = true;

    //Ranges used in actual code
    private final float[] angRange = new float[] {0, 0};
    private final float[] distRange = new float[] {0, 0};

    //Ranges where the users input is stored
    private final float[] angInput = new float[] {0, 0};
    private final float[] distInput = new float[] {0, 0};

    private int orders = 4; //This is the number or orders (including 0th)
    private final int orderColumns = 6; //Collums in final result per order
    private int totalCollumns = orders * orderColumns;

    // Variables related to checking for intensities
    private boolean checkForIntensity = false;
    private float ratioIntensity = 1.2f;

    // Variables and bool to remove pairs with few neighbours around in the same order
    private boolean toCleanup = false;
    private int neighbours = 7;
    private float cleanDistance = 100;

    private OwnColorTable ownColorTable; // Class to load LUT's

    //variables related to searching for the best angle
    // permutations to try next, storing results, what retry we are in
    // Uses tail recursion, sort of
    // (if its stupid and it works, its slightly less stupid)
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

    // Variables related to visualisation
    private boolean visualisation = true;
    private boolean visualiseZOLA = false;
    private String defaultLUT = "Spectrum.lut";
    private final float[] lutRange = new float[]{0,0};
    private float binwidth = 2.5f; //Binwidth in units used in distance for histograms
    //private final int[] unit_decades = {0, 0, 0, -2, -3, -6, -9, -10, -12, -15};

    // Core count as set by ImageJ preferences
    private final int coreCount = Prefs.getThreads();

    //Debug variables (only ever set or used when running from the IDE)
    private boolean processing = true;
    private static boolean debug = false;
    private static String debug_arg_string = "";
    private static boolean runningFromIDE = false;

    private static String GetInputStream(InputStream is) {
        StringBuilder result = new StringBuilder();
        try (InputStreamReader streamReader =
                     new InputStreamReader(is, StandardCharsets.UTF_8);
             BufferedReader reader = new BufferedReader(streamReader)) {

            String line;
            while ((line = reader.readLine()) != null) {
                result.append(line);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        return result.toString();
    }

    // Setup() ensures all variables get filled that need filling
    // and throws errors if it doesn't work out
    private boolean setup() {
        // During debugging we dont need to set any
        // When doing another check we have already set them all
        if (debug || doingRetry) return true;

        String arg; // Holds the macro arguments

        // Loading the debug string if its set
        if (debug_arg_string.equals("")) {
            arg = Macro.getOptions();
        } else {
            arg = debug_arg_string;
        }

        // If macro arguments are found, do not run any UI
        if (arg != null && !arg.equals("")) {
            // Remove some defaults
            saveSCV = false;
            visualisation = false;
            searchAngle = false;


            String[] arguments = arg.split(" ");

            // All accepted keywords
            String[] keywords = {
                    "csv_in", "csv_out","angle_start", "angle_end",
                    "distance_start", "distance_end",
                    "order_number", "check_order_intensity", "check_order_ratio",
                    "angle_flip", "angle_mirror", "angle_search", "angle_deep_search",
                    "lone_pair_remove", "lone_pair_neighbours", "lone_pair_distance",
                    "visualisation", "visualisationZOLA", "hist_binwidth", "LUT", "LUT_start", "LUT_end"
            };
            // For each keyword find the right variable and set it
            // if not found, or the value is malformed, throw an error
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
                                    if(deepSearchAngle) searchAngle = true;
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
                                    visualisation = Boolean.parseBoolean(keyword_val[1]);
                                    break;
                                case "visualisationZOLA":
                                    visualiseZOLA = Boolean.parseBoolean(keyword_val[1]);
                                    break;
                                case "hist_binwidth":
                                    binwidth = Float.parseFloat(keyword_val[1]);
                                    break;
                                case "LUT":
                                    defaultLUT = keyword_val[1];
                                    break;
                                case "LUT_start":
                                    lutRange[0] = Float.parseFloat(keyword_val[1]);
                                    break;
                                case "LUT_end":
                                    lutRange[1] = Float.parseFloat(keyword_val[1]);
                                    break;
                                default:
                                    logService.error("Keyword " + keyword_val[0] + " not found\nDid you mean: " + getTheClosestMatch(keywords, keyword_val[0]) + "?");
                                    return false;
                            }
                        } catch (ArrayIndexOutOfBoundsException e) {
                            logService.error("Malformed token: " + a + ".\nDid you remember to format it as keyword=value?");
                            return false;
                        } catch (Exception e){
                            logService.error("Failed to parse argument:" + a);
                            return false;
                        }
                    } else {
                        logService.error("Malformed token: " + a + ".\nDid you remember to format it as keyword=value?");
                        return false;
                    }
                }
            }

            // Check if the LUT provided works, otherwise try the default
            if(visualisation) {
                try {
                    ownColorTable.setLut(defaultLUT);
                } catch (Exception e) {
                    logService.info("Failed to set LUT.\nTrying Default: NCSA PalEdit/6_shades.lut");
                    try {
                        ownColorTable.setLut("NCSA PalEdit/6_shades.lut");
                    } catch (Exception e2) {
                        logService.info("Failed to set LUT again?.\n");
                        return false;
                    }
                }
            }
        } else {
            InputStream is = sSMLMA.class.getResourceAsStream("/Readme.txt"); // OK
            String content = GetInputStream(is);

            // Declare and fill the UI with all options
            GenericDialogPlus gd = new GenericDialogPlus("settings");
            String[] colors = ownColorTable.getLuts();
            {
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
                gd.addToSameRow();
                gd.addCheckbox("Visualise using ZOLA-3D?", visualiseZOLA);

                gd.addNumericField("Histogram binwidth", binwidth);

                Arrays.sort(colors);


                if (runningFromIDE) defaultLUT = "NCSA PalEdit/royal.lut";
                gd.addChoice("LUT", colors, defaultLUT);

                gd.addMessage("Modifies the LUT range of the distances. Will use the calculated distance ranges if not set.");
                gd.addNumericField("Start LUT", lutRange[0]);
                gd.addToSameRow();
                gd.addNumericField("End LUT", lutRange[1]);

                gd.addHelp(content);
            }

            gd.showDialog();

            if (gd.wasCanceled())
                return false;

            // Read all the resulting values into the variables
            // Type checking is done by the UI (doesnt allow alphanumerical in number field)
            {
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
                visualiseZOLA = gd.getNextBoolean();
                binwidth = (float) gd.getNextNumber();



                lutRange[0] = (float) gd.getNextNumber();
                lutRange[1] = (float) gd.getNextNumber();
            }

            // Check if the chosen LUT can be set
            if(visualisation) {
                try {
                    defaultLUT = colors[gd.getNextChoiceIndex()];
                    ownColorTable.setLut(defaultLUT);
                } catch (Exception e) {
                    logService.info("Failed to set LUT.\nTrying Default: NCSA PalEdit/6_shades.lut");
                    try {
                        ownColorTable.setLut("NCSA PalEdit/6_shades.lut");
                    } catch (Exception e2) {
                        return false;
                    }

                }
            }
        }


        if(deepSearchAngle) searchAngle = true; // If we do a deep search, this includes the normal search

        // Require input CSV
        if(filePath.equals("")){
            logService.error("No input CSV was set");
            return false;
        }

        // Require output directory if you want to save
        if(saveSCV && csv_target_dir.equals("")){
            logService.error("Set saving to CSV but no filepath was provided.");
            return false;
        }

        // Require user to either save or show results
        if (!(visualisation || saveSCV || visualiseZOLA)) {
            logService.error("No output method of any sorts is selected.\nSelect either Visualisation or Save to SCV.");
            return false;
        }

        return true;
    }

    @Override
    public void run() {

        ownColorTable = new OwnColorTable(lutService); // Instantiate colortable with our lutservice

        // If we are doing a retry, we need to reset the colorTable
        // Should (should) never fail since it succeeded before.
        if (doingRetry) {
            try {
                ownColorTable.setLut(defaultLUT);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        // Setup is done here
        // If it returns a failure (false) we execute nothing else
        if (setup()) {
            double csvTime = System.nanoTime(); // Timing loading csv

            List<String> collumns = new ArrayList<>(); // will store the collumns as found in csv
            totalCollumns = orders * orderColumns; // total number of collumns in result table

            // Debug stuff for skipping calculation and testing visualisation
            if (debug) {
                filePath = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\results\\thunderSTORM.csv";
                processing = false;
                visualisation = false;
                visualiseZOLA = true;
                saveSCV = false;
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

            // If we are not doing a retry we need to parse the CSV
            if (!doingRetry) {

                // Load our file into a matrix and retrieve the collumheaders
                // also catch any errors that might arise
                try {
                    floatMatrix = ownFloatMatrixLoader.loadCSVFile(filePath);
                    collumns = ownFloatMatrixLoader.getColumns();
                } catch (IOException e) {
                    logService.info("File not found.");
                } catch (LapackException e) {
                    logService.info("Lapack error");
                    e.printStackTrace();
                } catch (Exception e) {
                    logService.info("Error reading file. Does the csv start with a header?");
                    e.printStackTrace();
                }


                // Initialise the units and header arrays
                revOptionsIndices = new int[possible_options.length];
                unitsIndices = new int[collumns.size()];

                // Regex pattern that extracts the header name and any possible units
                Pattern pattern = Pattern.compile("(\\w+)( [ (\\[](\\w+)[)\\] ])?");

                // For each collumn get the header and unit it is using
                for (int i = 0; i < collumns.size(); i++) {
                    String header = collumns.get(i);
                    Matcher matcher = pattern.matcher(header);
                    if (matcher.find()) {
                        revOptionsIndices[getTheClosestMatch(possible_options, matcher.group(1))] = i;
                        unitsIndices[i] = getTheClosestMatch(unit_prefixes, matcher.group(3));
                    }
                }

                // Feedback about CSV loading time
                csvTime = System.nanoTime() - csvTime;
                System.out.println("Loading CSV took " + String.format("%.3f", csvTime / 1000000000) + " s");

            }

            double processingTime = System.nanoTime();
            final FloatMatrix data; // Holds all data (including datapoints that are not pairs)
            FloatMatrix finalPossibilities; // Holds all pairs

            // This boolean indicates if the angle finding succeeded or failed, along side some other checks
            boolean succes = false;


            if (processing && floatMatrix != null) {
                // Echo back settings used when searching for the best angle
                if (searchAngle) {
                    System.out.println("Run: " + runNumber + ". Determining Angle with settings: Flip Angle: " + flipAngles + ", Mirror Angle: " + mirrorAngles);
                }

                // frame, x, y, intensity
                // Load the relevant data into
                data = floatMatrix.getColumns(new int[]{revOptionsIndices[1], revOptionsIndices[2], revOptionsIndices[3], revOptionsIndices[5]});

                // If any var is not set, we have to calculate them all
                if (angRange[0] * angRange[1] * distRange[0] * distRange[1] == 0) {

                    // Instantiate the analyzer with the proper data and run it
                    AngleAnalyzer<T> angleAnalyzer = new AngleAnalyzer<>(data, flipAngles, mirrorAngles, logService, debug);
                    angleAnalyzer.run();

                    //get the results
                    float[] angResult = angleAnalyzer.getAngles();
                    float[] distResult = angleAnalyzer.getDistances();


                    // Parse all the results, only overwriting values not set by the user
                    angRange[0] = angInput[0] == 0f ? angResult[0] : angInput[0];
                    angRange[1] = angInput[1] == 0f ? angResult[1] : angInput[1];

                    distRange[0] = distInput[0] == 0f ? distResult[0] : distInput[0];
                    distRange[1] = distInput[1] == 0f ? distResult[1] : distInput[1];

                    lutRange[0] = distRange[0] == 0f ? distResult[0] : lutRange[0];
                    lutRange[1] = distRange[1] == 0f ? distResult[1] : lutRange[1];

                    // check if all calculations went right
                    succes = angleAnalyzer.getSucces();
                    if (distRange[0] > distRange[1])  succes = false;

                } else {
                    // if all values are filled in this is the only check we need to do
                    succes = !(distRange[0] > distRange[1]);

                }
            } else {
                // If we dont process (during debug) we directly load the data
                data = floatMatrix;
            }

            // If we failed (and are processing) we report the error and stop execution
            if ((!succes && processing) | floatMatrix == null){
                if (distRange[0] > distRange[1]) {
                    logService.error("The distance had to be positive: " + distRange[0] + " is larger than " + distRange[1]);
                } else {
                    logService.error("No features were detected. Are there pairs in this sample?");
                }

            } else {
                // Declaring a few vars due to scope juggling
                FloatMatrix halfOrderMatrix = null;
                FloatMatrix allOrdersCombined = null;
                if (processing) {

                    final int frames = (int) data.getColumn(0).max(); //Get end frame number
                    final int startFrame = (int) data.getColumn(0).min(); //Get start frame number

                    // Echo back amount of frames and points
                    System.out.println("Total Frames: " + (frames - startFrame));
                    System.out.println("Total Points: " + (int) floatMatrix.getColumn(revOptionsIndices[0]).max());

                    final AtomicInteger ai = new AtomicInteger(startFrame); //Atomic Integer is a thread safe incremental integer
                    final Thread[] threads = createThreadArray(coreCount); //Get array of threads

                    final FloatMatrix[] intermediateFinals = new FloatMatrix[coreCount]; // All intermediate results to be merged later
                    // prepare final verions of variables
                    final boolean finalIntensityCheck = checkForIntensity;
                    final float finalRatioIntensity = ratioIntensity;
                    //Set the run function for each thread
                    for (int ithread = 0; ithread < threads.length; ithread++) {

                        final int finalIthread = ithread;
                        threads[finalIthread] = new Thread(() -> {

                            // Will hold the final values for this thread
                            intermediateFinals[finalIthread] = new FloatMatrix(0, totalCollumns);

                            // Process each frame
                            for (int frame = ai.getAndIncrement(); frame <= frames; frame = ai.getAndIncrement()) {

                                // Showing process to the user
                                if (runningFromIDE && frame % 1000 == 0) logService.info("\r" + frame + "/" + frames);
                                IJ.showProgress(frame, frames);
                                IJ.showStatus(frame + "/" + frames);


                                final int[] frameIndicices = data.getColumn(0).eq(frame).findIndices(); // All rows indices for current frame

                                final FloatMatrix frameData = data.getRows(frameIndicices); //All rows for current frame

                                // Distance Matrices that show distance from one point to each other point in X and Y
                                final FloatMatrix subtractedX = makeSubstractedMatrix(frameData.getColumn(1));
                                final FloatMatrix subtractedY = makeSubstractedMatrix(frameData.getColumn(2));

                                // Matrices that have distances and angles from any point to another point in this frame
                                final FloatMatrix distances = Distance(subtractedX, subtractedY);
                                final FloatMatrix angles = atan2(subtractedX, subtractedY, distances);

                                // Filter to only points that match the distance and angle restrictions;
                                // The value indicates both the start and end point
                                final int[] correctAngleAndDistance = distances.gt(distRange[0]).and(distances.lt(distRange[1]))
                                        .and(angles.gt(angRange[0])).and(angles.lt(angRange[1])).findIndices();

                                // If points were found, we must process them
                                if (correctAngleAndDistance.length > 1) {

                                    // Matrices to hold the possibilities and combined ones
                                    final FloatMatrix possibilities = new FloatMatrix(correctAngleAndDistance.length, totalCollumns);
                                    FloatMatrix intermediateFinalPossibilities = new FloatMatrix(0, totalCollumns);

                                    // For each pair, record it to out temporary matrix
                                    for (int i = 0; i < correctAngleAndDistance.length; i++) {
                                        int index = correctAngleAndDistance[i];

                                        // Check if the intensity ratio checks out, if that check is enabled
                                        if (finalIntensityCheck &&
                                                (frameData.get(index / distances.rows, 3) /
                                                        frameData.get(index % distances.rows, 3)) > finalRatioIntensity
                                        ) continue;

                                        // Put all info about the two points into our matrix
                                        possibilities.putRow(i, extend(new FloatMatrix(1, orderColumns * 2,
                                                0,                                                     //0 (global index goes here later)
                                                frame,                                                          //1 frame
                                                1 + Math.floorDiv(index, distances.rows),                       //2 start index in frame
                                                frameData.get(index / distances.rows, 1),    //3 start x
                                                frameData.get(index / distances.rows, 2),    //4 start y
                                                frameData.get(index / distances.rows, 3),    //5 start intensity
                                                1 + Math.floorMod(index, distances.rows),                       //6 end index in frame
                                                frameData.get(index % distances.rows, 1),    //7 end x
                                                frameData.get(index % distances.rows, 2),    //8 end y
                                                frameData.get(index % distances.rows, 3),    //9 end intensity
                                                distances.get(index),                                           //10 distance
                                                angles.get(index)                                               //11 angle
                                        ), 1, totalCollumns));
                                    }

                                    FloatMatrix connectedIds = possibilities.getColumn(6); // all frame indices to what a point is connected
                                    List<Integer> toSkip = new ArrayList<>(); // holds indices to skip over

                                    // For each pair in possibilities, ensure there are none that have overlapping starts or ends
                                    // (1-3, 2-3) ->  4-3 where 4 is average position of 2 and 3
                                    for (int i = 0; i < possibilities.rows; i++) {
                                        if (!toSkip.contains(i)) { // Skip it if we processed it already

                                            //If more ids match the target we average the position/intensity and skip any of the matches
                                            //If there is only one, we add it to our final array
                                            FloatMatrix identicalIds = connectedIds.eq(connectedIds.get(i));
                                            float identicalIdSum = identicalIds.sum();
                                            if (identicalIdSum > 1.0f) {

                                                float cumX = possibilities.get(i, 3);
                                                float cumY = possibilities.get(i, 4);
                                                float cumInt = possibilities.get(i, 5);

                                                for (int index2 : identicalIds.findIndices()) {
                                                    if (index2 != i) {
                                                        cumX += possibilities.get(index2, 3);
                                                        cumY += possibilities.get(index2, 4);
                                                        cumInt += possibilities.get(index2, 5);
                                                        toSkip.add(index2);
                                                    }
                                                }

                                                FloatMatrix temp_arr = possibilities.getRow(i).dup();
                                                temp_arr.put(3, cumX / identicalIdSum);
                                                temp_arr.put(4, cumY / identicalIdSum);
                                                temp_arr.put(5, cumInt / identicalIdSum);

                                                intermediateFinalPossibilities = FloatMatrix.concatVertically(intermediateFinalPossibilities, temp_arr);

                                            } else {
                                                intermediateFinalPossibilities = FloatMatrix.concatVertically(intermediateFinalPossibilities, possibilities.getRow(i));
                                            }
                                        }
                                    }
                                    // Add the points calculated in this frame to this threads buffer
                                    // We also connect the orders here: (1-2, 2-3, 3-4 -> 1-2-3-4)
                                    intermediateFinals[finalIthread] = FloatMatrix.concatVertically(intermediateFinals[finalIthread], connectOrders(intermediateFinalPossibilities, orders, orderColumns));

                                } else if (correctAngleAndDistance.length == 1) { // We found only a single point here, so we just add it and do no connecting
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

                    // Start the threads
                    // This actually starts the calculations
                    startAndJoin(threads);

                    //Combine all Matrices together
                    finalPossibilities = new FloatMatrix(0, totalCollumns);
                    for (FloatMatrix floatMatrix1 : intermediateFinals) {
                        finalPossibilities = FloatMatrix.concatVertically(finalPossibilities, floatMatrix1);
                    }

                    // Add global index for each row
                    int id = 0;
                    for (int i = 0; i < finalPossibilities.rows; i++) {
                        finalPossibilities.put(i, 0, id++);
                    }

                    // Echo back time it took
                    processingTime = System.nanoTime() - processingTime;
                    System.out.println("\nProcessing data took " + String.format("%.3f", processingTime / 1000000000) + " s");


                    // Ensure nothing went wrong and echo back how many points we found
                    // Also clean up some garbage since we are done processing and there are many things we no longer need
                    System.out.println("Pairs in the 0th-1st order found: " + finalPossibilities.rows);
                    System.gc();

                    // If we find very few or no pairs, we echo this back and set the retry flag indicating we should do another attempt
                    // Or tell the user they could enabling searching
                    if (finalPossibilities.rows < 10) {
                        logService.info("No pairs detected. Try flipping the angle or manually adjusting the distance and angles.");
                        retry = true;
                    }
                    if (finalPossibilities.rows < 0.1 * data.rows) {
                        if (searchAngle)
                            logService.info("Very few pairs were detected.");
                        else
                            logService.info("Very few pairs were detected. A search for a better angle could improve results.");
                        retry = true;
                    }

                    // This tells us some information about the angles, like if it has a guassian shape, or if one of the tails is cutoff
                    boolean[][] checks = checkForRetry(finalPossibilities.getColumn(11));

                    if (sum(checks[0]) < 2) {
                        logService.info("Nothing resembling a guassian was found for the angles.");
                        retry = true;
                    }

                    if (checks[1][0]) {
                        if (checks[1][1])
                            logService.info("The left tail of the angle histogram seems to be partially cut-off");
                        else logService.info("The right tail of the angle histogram seems to be partially cut-off");
                        logService.info("You can take the angle value given and adjust it manually.");
                    }

                    // If we are searching for a better angle, save the results and determine if another retry is needed
                    if ((retry && searchAngle) | deepSearchAngle) {
                        boolean[] curr_perm = new boolean[]{flipAngles, mirrorAngles};

                        // Save our results in the correct position and load the next permutation of variables
                        for (int i = 0; i < perm.length; i++) {
                            if (Arrays.equals(perm[i], curr_perm)) {
                                angleResults[i] = new FloatMatrix();
                                angleResults[i].copy(finalPossibilities);
                                perm[i] = null;

                                if (perm[(i + 1) % perm.length] != null) {
                                    flipAngles = perm[(i + 1) % perm.length][0];
                                    mirrorAngles = perm[(i + 1) % perm.length][1];
                                }
                                break;
                            }
                        }

                        // Check if any permutation had not been attempted yet
                        boolean all_filled = true;
                        for (FloatMatrix result : angleResults) {
                            if (result == null) {
                                all_filled = false;
                                break;
                            }
                        }
                        // If we haven't done all permutations of settings we arent done yet
                        if (!all_filled) {
                            // We reload any user input into the varialbes and do a retry
                            // This runs run() again
                            // this is not good, but it was the best solution
                            // for now
                            doingRetry = true;

                            angRange[0] = angInput[0] == 0f ? angInput[0] : angRange[0];
                            angRange[1] = angInput[1] == 0f ? angInput[1] : angRange[1];

                            distRange[0] = distInput[0] == 0f ? distInput[0] : distRange[0];
                            distRange[1] = distInput[1] == 0f ? distInput[1] : distRange[1];

                            runNumber++;
                            run();

                        } else {
                            // We are done and load the best results back into out matrix, as well as the best settings
                            int max = angleResults[0].rows;
                            int ind = 0;
                            for (int i = 1; i < angleResults.length; i++) {
                                if (angleResults[i].rows > max) {
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

                            // Echo back the settings found for this best permutations
                            String message = "\nBest result was " + max + " pairs found.\n" +
                                    "The following values were found:\n" +
                                    "\tAngle(rad): " + angRange[0] + " to " + angRange[1] + "\n" +
                                    "\tDistance: " + distRange[0] + " to " + distRange[1] + "\n" +
                                    "With the settings:\n" +
                                    "\tFlip Angle: " + permReference[ind][0] + "\n" +
                                    "\tMirror Angle: " + permReference[ind][1] + "\n";
                            if (saveSCV) message += ("CSV files were saved to the folder: " + csv_target_dir + "\n");

                            if (runningFromIDE) System.out.println(message);
                            IJ.showMessage(message);
                        }

                    }


                    // Anything after this point is skipped if we are not in the final run
                    // So this point is only reached with the best (hopefully) results
                    if ((!(retry && searchAngle) | (foundBestResult && deepSearchAngle)) && displayInfo) {
                        displayInfo = false; // ensures this path is only ran once

                        // We remove any points with too few neighbours here, if enabled
                        if (processing && toCleanup) {
                            IJ.showStatus("Cleaning Data");
                            logService.info("Cleaning up data");

                            FloatMatrix backup = finalPossibilities.dup();

                            // This sometimes fails and i have not been able to determine why
                            // A second attempt seems to always works, somehow
                            try {
                                finalPossibilities = cleanup(finalPossibilities, neighbours, cleanDistance, coreCount);
                            } catch (Exception e) {
                                try {
                                    finalPossibilities = cleanup(finalPossibilities, neighbours, cleanDistance, coreCount);
                                } catch (Exception e2) {
                                    logService.error("Cleaning up points failed.");
                                    logService.info("Continuing without cleanup.");

                                    finalPossibilities = backup; // In case the matrix got corrupted above

                                }
                            }
                        }


                        // We have our final data
                        // but some other versions could give improved information by combining position data from multiple orders into one
                        // We do this into two ways here
                        // Once with a only one pair (0th-1st) order
                        // Another one with all orders combined, as many as there are for one row or data
                        // It is a lot of code doing some simple combining
                        //id, frame, x, y, intensity, distance
                        halfOrderMatrix = new FloatMatrix(finalPossibilities.rows, orderColumns);
                        allOrdersCombined = new FloatMatrix(finalPossibilities.rows, orderColumns);
                        {


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


                            for (int i = 0; i < orders - 2; i++) {
                                FloatMatrix relevantRows = finalPossibilities.getColumn(12 + (i * orderColumns)).ne(0.0f);
                                offsets.getColumn(0).addi(finalPossibilities.getColumn(13 + (i * orderColumns)).subi(finalPossibilities.getColumn(7 + (i * orderColumns))).divi(2.0f).muli(relevantRows));
                                offsets.getColumn(1).addi(finalPossibilities.getColumn(14 + (i * orderColumns)).subi(finalPossibilities.getColumn(8 + (i * orderColumns))).divi(2.0f).muli(relevantRows));
                            }


                            allOrdersCombined.putColumn(2, finalPossibilities.getColumn(3));
                            allOrdersCombined.putColumn(3, finalPossibilities.getColumn(4));

                            FloatMatrix relevantRows = finalPossibilities.getColumn((orders - 1) * orderColumns).ne(0.0f);
                            for (int i = orders; i > 1; i--) {
                                allOrdersCombined.getColumn(2).addi(offsets.getColumn(0).divi((float) i).muli(relevantRows));
                                allOrdersCombined.getColumn(3).addi(offsets.getColumn(1).divi((float) i).muli(relevantRows));

                                relevantRows.xori(finalPossibilities.getColumn((i - 2) * orderColumns).ne(0.0f));
                            }
                        }
                    }
                } else { // If we skip all processing in debug mode we are done
                    finalPossibilities = floatMatrix;
                    halfOrderMatrix = floatMatrix;
                    allOrdersCombined = floatMatrix;
                }
                if (visualisation) {
                    assert allOrdersCombined != null;

                    HistogramWindow[] histograms = new HistogramWindow[orders - 1]; // We create some histograms for each distance order we want to visualise

                    // For each order we fill the distance into each matrix
                    // this allows one to compare the distances for each order
                    // they should be the same for each order to order, but its still useful
                    for (int i = 0; i < orders - 1; i++) {
                        FloatMatrix relevantData;

                        // Visualise either every row or only the rows that have data in it
                        // Since not each row might have data for each order
                        if (i == 0) {
                            relevantData = finalPossibilities.getColumn(10);
                        } else {
                            relevantData = finalPossibilities.getColumn(10 + (i * orderColumns));
                            relevantData = relevantData.get(relevantData.ne(0.0f).findIndices());
                        }

                        System.out.println("Found " + relevantData.rows + " connections in the " + getTitleHist(i) + " order");

                        // If there are few points, the graph is useless and we do not display it and do not visualise more orders than this
                        if (relevantData.rows < 50) {
                            orders = i + 1;
                            break;
                        }
                        // Create an image from the matrix, and then use that to make a histogram
                        // Thanks ImageJ
                        ImagePlus dummy = new ImagePlus("", new FloatProcessor(relevantData.toArray2()));
                        histograms[i] = new HistogramWindow(getTitleHist(i), dummy, getBins(relevantData, binwidth), distRange[0], distRange[1]);
                        if (runningFromIDE) histograms[i].getImagePlus().show();
                    }
                    ///////////////// Angles

                    // Take all the angles for the first pair, create an image and display it
                    // width of bins is hardcoded to 0.005 rad
                    ImagePlus Angles = new ImagePlus("", new FloatProcessor(finalPossibilities.getColumn(11).toArray2()));
                    HistogramWindow angleHist = new HistogramWindow("Angles", Angles, getBins(finalPossibilities.getColumn(11), 0.005f), angRange[0], angRange[1]);
                    if (runningFromIDE) angleHist.getImagePlus().show();

                    ///////////////////////////////////////////////////////////// PLOTS
                    // Arrays for the colors and shapes used for each order
                    // Quick reference to all possible shapes
                    // "line", "connected circle", "filled", "bar", "separated bar", "circle", "box", "triangle", "diamond", "cross", "x", "dot", "error bars" or "xerror bars"
                    String[] colors = {"blue", "red", "green", "black"};
                    String[] shapes = {"dot", "dot", "dot", "dot"};

                    //////////////////////////////////////////////

                    // This plot displays the x and y for all found points
                    Plot orderPlot = new Plot("Orders", "x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]", "y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]");


                    // Set the color and add the first order (and then all other orders in the for loop)
                    orderPlot.setColor(colors[0]);
                    orderPlot.add(shapes[0], toDouble(finalPossibilities.getColumn(3)), toDouble(finalPossibilities.getColumn(4)));
                    for (int i = 1; i < orders; i++) {
                        orderPlot.setColor(colors[i]);
                        orderPlot.add(shapes[i], toDouble(finalPossibilities.getColumn(1 + (i * orderColumns))), toDouble(finalPossibilities.getColumn(2 + (i * orderColumns))));
                    }

                    orderPlot.setLegend(getTitlePlot(orders), Plot.AUTO_POSITION); // Set the legend and its position
                    orderPlot.show();
                    //////////////////////////////////////////////

                    // Show a Plot with only the combined position for each row, reducing the amounts of points
                    // The color of each point is the distance between the 0th and 1st order, based on the selected LUT

                    Plot distancePlot = new Plot("Distance", "x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]", "y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]");

                    for (int i = 0; i < finalPossibilities.rows; i++) {

                        distancePlot.setColor(ownColorTable.getColor(allOrdersCombined.get(i, 5), distRange[0], distRange[1]));

                        distancePlot.add(shapes[0], toDouble(halfOrderMatrix.get(i, 2)), toDouble(halfOrderMatrix.get(i, 3))); // 3 4
                    }

                    distancePlot.setLimitsToFit(true); // Ensure all points are visible
                    addLutLegend(distancePlot, ownColorTable, "Distance", 512, distRange[0], distRange[1]); // Add the LUT as a legend
                    distancePlot.show();
                    distancePlot.setLimits(Float.NaN, Float.NaN, Float.NaN, Float.NaN); // Ensure all points are visible, again



                }

                if (saveSCV) {
                    assert halfOrderMatrix != null;
                    //Create Header
                    // Longheader for all the data
                    // Shortheader when tis just one position
                    List<String> LongHeader = new ArrayList<>();
                    List<String> ShortHeader = new ArrayList<>();
                    LongHeader.add("id");
                    ShortHeader.add("id");
                    LongHeader.add("frame");
                    ShortHeader.add("frame");

                    // Pre-create the distance unit for the distance
                    String distanceUnit = (unit_prefixes[unitsIndices[revOptionsIndices[2]]].equals(unit_prefixes[unitsIndices[revOptionsIndices[3]]]) ? unit_prefixes[unitsIndices[revOptionsIndices[2]]] : ("(" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "*" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + ")^½"));
                    // Add the headers for each order
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
                    // Only add one header level to the short ones
                    ShortHeader.add("x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]");
                    ShortHeader.add("y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]");
                    ShortHeader.add("intensity [" + unit_prefixes[unitsIndices[revOptionsIndices[5]]] + "]");
                    ShortHeader.add("z [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]");

                    // Save all data using the proper header, including one that easily is loaded into ThunderSTORM again for visualisation etc
                    SaveCSV(finalPossibilities, LongHeader, csv_target_dir + "\\all_orders.csv");

                    SaveCSV(halfOrderMatrix, ShortHeader, csv_target_dir + "\\accurate_positions.csv");
                    SaveCSV(finalPossibilities.getColumns(new int[]{0, 1, 3, 4, 5, 10}), ShortHeader, csv_target_dir + "\\thunderSTORM.csv");
                }

                if(visualiseZOLA) {
                    logService.info("ZOLA Visualisation");

                    try {
                        String tmpfile;

                        if (runningFromIDE) {
                            System.out.println("Running from IDE does not work for ZOLA integration");
                        } else {
                            if (saveSCV) {
                                tmpfile = csv_target_dir + "\\thunderSTORM.csv";
                            } else {
                                tmpfile = IJ.getDirectory("temp") + "tmp.csv";
                                saveThunderSTORM(tmpfile, finalPossibilities.getColumns(new int[]{0, 1, 3, 4, 5, 10}));
                            }


                            Prefs.set("Zola.showLUT", true);
                            Prefs.set("Zola.pathlocalization", tmpfile);
                            Prefs.set("Zola.is3Drendering", true);

                            IJ.run("Import table");
                            IJ.run("2D/3D histogram");
                        }
                    }  catch (Exception e) {
                        logService.info("ZOLA integration failed");
                        e.printStackTrace();
                    }
                }


            }
        }
    }


    // Only run from the IDE
    public static void main(String[] args) {
        debug = true;
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

        //debug_arg_string = "csv_in=F:\\ThesisData\\output\\output3_drift.csv  visualisation=true";
        debug_arg_string = "csv_in=F:\\ThesisData\\output\\output3_drift.csv angle_start=-0.094 angle_end=0.22 distance_start=1500 distance_end=2200 visualisation=true";
        //debug_arg_string = "csv_in=F:\\ThesisData\\Test3D\\localisations_drift.csv  visualisation=true";

        net.imagej.ImageJ ij = new ImageJ();
        ij.ui().showUI();

        ij.command().run(sSMLMA.class, true);

    }
}
