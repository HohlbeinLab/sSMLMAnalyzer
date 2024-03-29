package com.wurgobes.sSMLMAnalyzer;

/*
Spectral Super Resolution Pair Finder
(c) 2021 Martijn Gobes, Wageningen University.
Based on the work by Koen Martens

This plugin works by either going to the menu path: Plugins>Spectral Analyzer>Analyze Pairs
or by using a macro (example provided in readme.md

The class is instantiated by ImageJ and run(String arg) is called.
If any angle/distance value is not set, AngleAnalyzer is called on the csv data to calculate these.

All angles are given in rad, and assume a 0 value of (->), in a range of [-Pi, Pi),
The atan2 algorithm in Util.java is used for this.

Pairs are created by filtering the options based on the angle and distance.
Pairs are then chained together to create multi-order connections (1-2, 2-3, 3-4 -> 1-2-3-4)

If multiple options for pairs are found, the best one is selected.
Ties are broken by angle; closest angle to the original one is chosen = pairs closest to straight line.

The pair chains are written to the output in 3 versions.
One contains all data (all_orders.csv) for each pair-chain, with position (x, y, z), intensity, distance and angle (rad).
The file ThunderSTORM.csv contains only the starting point, where the z position is substituted by distance,
this allows one to load this into various plugins and visualise the z dimension as distance.

The last file, accurate_positions.csv, contains combined positional data.
It contains only a single point per pair-chain, but one where the positional data of all points in it's pair chain is combined.
Sometimes aberrations were noticed where the first points where right shifted, and the second points where left shifted.
This file takes the delta x and delta y between each pair and takes the average.
This value is halved (to get the midway point) and added to the first point in the pair-chain.
The distance and angle is also averaged.

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


import com.wurgobes.sSMLMAnalyzer.CustomPlot.CustomPlot;
import ij.*;
import static ij.util.ThreadUtil.*;
import ij.gui.HistogramWindow;
import ij.gui.Plot;
import ij.process.FloatProcessor;

import net.imagej.ImageJ;
import net.imagej.lut.LUTService;

import net.imglib2.type.numeric.IntegerType;

import org.jblas.util.Random;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.*;

import java.io.*;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;


import fiji.util.gui.GenericDialogPlus;

import org.jblas.FloatMatrix;
import org.jblas.exceptions.LapackException;

import static com.wurgobes.sSMLMAnalyzer.Util.*;
import static com.wurgobes.sSMLMAnalyzer.levenshtein.getTheClosestMatch;

import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.dom.*;
import javax.xml.transform.stream.*;
import org.w3c.dom.*;

@Plugin(type = Command.class, menuPath = "Plugins>Spectral Analyzer>Analyze Pairs")

// TODO
// refactoring:
//   No concating in calculations
// create a proper final FloatMatrix type, but performance benefit is unknown

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
    private final int orderColumns = 7; //Columns in final result per order (messing with this will break many many things)
    private int totalColumns = orders * orderColumns;

    // Variables related to checking for intensities and distances
    private boolean checkforZ = false;
    private float zMargin = 10f;
    private boolean hasZ = false;
    private boolean hasIntensity = false;
    private boolean checkForIntensity = false;
    private float ratioIntensity = 1.2f;
    private boolean checkDistanceOrderDelta = false;
    private float distanceDelta = 100f;

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
    private String defaultLUT = "physics.lut";
    private final float[] lutRange = new float[]{0,0};
    private float binwidth = 2.5f; //Binwidth in units used in distance for histograms
    //private final int[] unit_decades = {0, 0, 0, -2, -3, -6, -9, -10, -12, -15};

    // Core count as set by ImageJ preferences
    private final int coreCount = Prefs.getThreads();

    //Debug variables (only ever set or used when running from the IDE)
    private boolean processing = true;
    private final static boolean debug = false;
    private static String debug_arg_string = "";
    private static boolean runningFromMacro = false;
    private static boolean runningFromIDE = false;
    private final static boolean debugFlag = false;

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

    private void WriteXML(Path FilePath) {
        Document dom;
        Element e1;
        Element e2;
        Element e3;

        // instance of a DocumentBuilderFactory
        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        try {
            // use factory to get an instance of document builder
            DocumentBuilder db = dbf.newDocumentBuilder();
            // create instance of DOM
            dom = db.newDocument();

            // create the root element
            Element rootEle = dom.createElement("Spectral_SMLM_Analzyer");

            e1 = dom.createElement("info");

            e2 = dom.createElement("version");
            e2.appendChild(dom.createTextNode("1.0")); //should make this dynamic
            e1.appendChild(e2);

            e2 = dom.createElement("time");
            e2.appendChild(dom.createTextNode(String.valueOf(java.time.LocalDate.now()))); //yyyy-MM-dd
            e1.appendChild(e2);

            rootEle.appendChild(e1);

            e1 = dom.createElement("settings");

            e2 = dom.createElement("csv_in");
            e2.appendChild(dom.createTextNode(filePath));
            e1.appendChild(e2);

            if(saveSCV) {
                e2 = dom.createElement("csv_out");
                e2.appendChild(dom.createTextNode(csv_target_dir));
                e1.appendChild(e2);
            }

            e2 = dom.createElement("angles_input");

                e3 = dom.createElement("start");
                e3.appendChild(dom.createTextNode(String.valueOf(angInput[0])));
            e2.appendChild(e3);

                e3 = dom.createElement("end");
                e3.appendChild(dom.createTextNode(String.valueOf(angInput[1])));
            e2.appendChild(e3);

            e1.appendChild(e2);

            e2 = dom.createElement("distances_input");

                e3 = dom.createElement("start");
                e3.appendChild(dom.createTextNode(String.valueOf(distInput[0])));
            e2.appendChild(e3);

                e3 = dom.createElement("end");
                e3.appendChild(dom.createTextNode(String.valueOf(distInput[1])));
            e2.appendChild(e3);

            e1.appendChild(e2);

            e2 = dom.createElement("orders");
            e2.appendChild(dom.createTextNode(String.valueOf(orders)));
            e1.appendChild(e2);

            if(checkForIntensity) {
                e2 = dom.createElement("ratio_intensity");
                e2.appendChild(dom.createTextNode(String.valueOf(ratioIntensity)));
                e1.appendChild(e2);
            }

            e2 = dom.createElement("angle_flip");
            e2.appendChild(dom.createTextNode(String.valueOf(flipAngles)));
            e1.appendChild(e2);

            e2 = dom.createElement("angle_mirror");
            e2.appendChild(dom.createTextNode(String.valueOf(mirrorAngles)));
            e1.appendChild(e2);

            e2 = dom.createElement("angle_search");
            e2.appendChild(dom.createTextNode(String.valueOf(searchAngle)));
            e1.appendChild(e2);

            e2 = dom.createElement("angle_deep_search");
            e2.appendChild(dom.createTextNode(String.valueOf(deepSearchAngle)));
            e1.appendChild(e2);

            if(toCleanup) {
                e2 = dom.createElement("lone_pair_neighbours");
                e2.appendChild(dom.createTextNode(String.valueOf(neighbours)));
                e1.appendChild(e2);

                e2 = dom.createElement("lone_pair_distance");
                e2.appendChild(dom.createTextNode(String.valueOf(cleanDistance)));
                e1.appendChild(e2);
            }

            if(visualisation) {
                e2 = dom.createElement("ZOLA");
                e2.appendChild(dom.createTextNode(String.valueOf(visualiseZOLA)));
                e1.appendChild(e2);

                e2 = dom.createElement("histogram_binwidth");
                e2.appendChild(dom.createTextNode(String.valueOf(binwidth)));
                e1.appendChild(e2);

                e2 = dom.createElement("LUT");
                e2.appendChild(dom.createTextNode(defaultLUT));
                e1.appendChild(e2);

                e2 = dom.createElement("LUT_range");

                    e3 = dom.createElement("start");
                    e3.appendChild(dom.createTextNode(String.valueOf(lutRange[0])));
                e2.appendChild(e3);

                    e3 = dom.createElement("end");
                    e3.appendChild(dom.createTextNode(String.valueOf(lutRange[1])));
                e2.appendChild(e3);

                e1.appendChild(e2);
            }

            if(checkforZ) {
                e2 = dom.createElement("check_z_margin");
                e2.appendChild(dom.createTextNode(String.valueOf(zMargin)));
                e1.appendChild(e2);
            }

            if(checkDistanceOrderDelta) {
                e2 = dom.createElement("distance_delta");
                e2.appendChild(dom.createTextNode(String.valueOf(zMargin)));
                e1.appendChild(e2);
            }

            rootEle.appendChild(e1);

            e1 = dom.createElement("calculated");

            e2 = dom.createElement("angles_calculated");

                e3 = dom.createElement("start");
                e3.appendChild(dom.createTextNode(String.valueOf(angRange[0])));
            e2.appendChild(e3);

                e3 = dom.createElement("end");
                e3.appendChild(dom.createTextNode(String.valueOf(angRange[1])));
            e2.appendChild(e3);

            e1.appendChild(e2);

            e2 = dom.createElement("distance_calculated");

                e3 = dom.createElement("start");
                e3.appendChild(dom.createTextNode(String.valueOf(distRange[0])));
            e2.appendChild(e3);

                e3 = dom.createElement("end");
                e3.appendChild(dom.createTextNode(String.valueOf(distRange[1])));
            e2.appendChild(e3);

            e1.appendChild(e2);

            e2 = dom.createElement("LUT_calculated");

            e3 = dom.createElement("start");
            e3.appendChild(dom.createTextNode(String.valueOf(lutRange[0])));
            e2.appendChild(e3);

            e3 = dom.createElement("end");
            e3.appendChild(dom.createTextNode(String.valueOf(lutRange[1])));
            e2.appendChild(e3);

            e1.appendChild(e2);

            rootEle.appendChild(e1);


            dom.appendChild(rootEle);

            try {
                Transformer tr = TransformerFactory.newInstance().newTransformer();
                tr.setOutputProperty(OutputKeys.INDENT, "yes");
                tr.setOutputProperty(OutputKeys.METHOD, "xml");
                tr.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
                tr.setOutputProperty(OutputKeys.DOCTYPE_SYSTEM, "settings.dtd");
                tr.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");
                // send DOM to file
                tr.transform(new DOMSource(dom),
                        new StreamResult(Files.newOutputStream(FilePath)));
            } catch (TransformerException | IOException te) {
                logService.info(te.getMessage());
            }
        } catch (ParserConfigurationException pce) {
            logService.info("UsersXML: Error trying to instantiate DocumentBuilder " + pce);
        }
    }

    // Setup() ensures all variables get filled that need filling
    // and throws errors if it doesn't work out
    private boolean setup() {
        // During debugging we dont need to set any
        // When doing another check we have already set them all
        if (debug || doingRetry) return true;

        String arg; // Holds the macro arguments
        ownColorTable = new OwnColorTable(lutService);

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
            runningFromMacro = true;

            final Pattern pattern = Pattern.compile("(\\w+)(=('[^']+'|\\S+))?");
            Matcher m = pattern.matcher(arg);

            // All accepted keywords
            String[] keywords = {
                    "csv_in", "csv_out","angle_start", "angle_end",
                    "distance_start", "distance_end",
                    "order_number", "check_order_intensity", "check_order_ratio",
                    "angle_flip", "angle_mirror", "angle_search", "angle_deep_search",
                    "lone_pair_remove", "lone_pair_neighbours", "lone_pair_distance",
                    "visualisation", "visualisationZOLA", "hist_binwidth", "LUT", "LUT_start", "LUT_end",
                    "check_z", "check_z_margin", "check_distance_delta", "distance_delta",
                    //These are for macro recording mode. Don't think about it
                    "browse", "csv", "save", "csv_0", "start", "end", "start_0", "end_0",
                    "number", "restrict", "max", "intensity", "ratio", "flip", "mirror",
                    "search", "search_0", "remove", "required_0", "remove_0", "maximum", "visualise",
                    "visualise_0", "histogram", "start_1", "end_1", "file_path"
            };
            String[] macroRecordingKeywords = {
                    "browse", "csv", "save", "csv_0", "start", "end", "start_0", "end_0",
                    "number", "restrict", "max", "intensity", "ratio", "flip", "mirror",
                    "search", "search_0", "remove", "required_0", "remove_0", "maximum", "visualise",
                    "visualise_0", "histogram", "start_1", "end_1", "file_path"
            };
            // For each keyword find the right variable and set it
            // if not found, or the value is malformed, throw an error
            logService.info(arg);
            while (m.find()) {
                if (m.groupCount() == 3 || Arrays.asList(macroRecordingKeywords).contains(m.group(1))) {
                    String[] keyword_val = {m.group(1), m.group(3) != null ? m.group(3).replace("'", "") : String.valueOf(false)};
                    try {
                        switch (keyword_val[0]) {
                            case "csv_in":
                            case "browse":
                            case "csv":
                                filePath = keyword_val[1].replace("\\\\", "\\").replace("\\","\\\\");
                                break;
                            case "csv_0":
                            case "csv_out":
                                if (!keyword_val[1].equals("[]")) {
                                    saveSCV = true;
                                    csv_target_dir = keyword_val[1].replace("\\\\", "\\").replace("\\","\\\\");
                                }
                                break;
                            case "start":
                            case "angle_start":
                                angInput[0] = (float) (Float.parseFloat(keyword_val[1]) * (Math.PI / 180f));
                                break;
                            case "end":
                            case "angle_end":
                                angInput[1] = (float) (Float.parseFloat(keyword_val[1]) * (Math.PI / 180f));
                                break;
                            case "start_0":
                            case "distance_start":
                                distInput[0] = Float.parseFloat(keyword_val[1]);
                                break;
                            case "end_0":
                            case "distance_end":
                                distInput[1] = Float.parseFloat(keyword_val[1]);
                                break;
                            case "number":
                            case "order_number":
                                orders = Integer.parseInt(keyword_val[1]);
                                break;
                            case "intensity":
                                checkForIntensity = true;
                                break;
                            case "check_order_intensity":
                                checkForIntensity = Boolean.parseBoolean(keyword_val[1]);
                                break;
                            case "ratio":
                            case "check_order_ratio":
                                ratioIntensity = Float.parseFloat(keyword_val[1]);
                                break;
                            case "flip":
                                flipAngles = true;
                                break;
                            case "angle_flip":
                                flipAngles = Boolean.parseBoolean(keyword_val[1]);
                                break;
                            case "mirror":
                                mirrorAngles = true;
                                break;
                            case "angle_mirror":
                                mirrorAngles = Boolean.parseBoolean(keyword_val[1]);
                                break;
                            case "search":
                                searchAngle = true;
                                break;
                            case "angle_search":
                                searchAngle = Boolean.parseBoolean(keyword_val[1]);
                                break;
                            case "search_0":
                                deepSearchAngle = true;
                                searchAngle = true;
                                break;
                            case "angle_deep_search":
                                deepSearchAngle = Boolean.parseBoolean(keyword_val[1]);
                                if (deepSearchAngle) searchAngle = true;
                                break;
                            case "remove":
                                toCleanup = true;
                                break;
                            case "lone_pair_remove":
                                toCleanup = Boolean.parseBoolean(keyword_val[1]);
                                break;
                            case "required":
                            case "lone_pair_neighbours":
                                neighbours = Integer.parseInt(keyword_val[1]);
                                break;
                            case "required_0":
                            case "lone_pair_distance":
                                cleanDistance = Float.parseFloat(keyword_val[1]);
                                break;
                            case "visualise":
                                visualisation = true;
                                break;
                            case "visualisation":
                                visualisation = Boolean.parseBoolean(keyword_val[1]);
                                break;
                            case "visualise_0":
                                visualiseZOLA = true;
                                break;
                            case "visualisationZOLA":
                                visualiseZOLA = Boolean.parseBoolean(keyword_val[1]);
                                break;
                            case "histogram":
                            case "hist_binwidth":
                                binwidth = Float.parseFloat(keyword_val[1]);
                                break;
                            case "LUT":
                                defaultLUT = keyword_val[1];
                                break;
                            case "start_1":
                            case "LUT_start":
                                lutRange[0] = Float.parseFloat(keyword_val[1]);
                                break;
                            case "end_1":
                            case "LUT_end":
                                lutRange[1] = Float.parseFloat(keyword_val[1]);
                                break;
                            case "restrict":
                                checkforZ = true;
                                break;
                            case "check_z":
                                checkforZ = Boolean.parseBoolean(keyword_val[1]);
                                break;
                            case "max":
                            case "check_z_margin":
                                zMargin = Float.parseFloat(keyword_val[1]);
                                break;
                            case "remove_0":
                                checkDistanceOrderDelta = true;
                                break;
                            case "check_distance_delta":
                                checkDistanceOrderDelta = Boolean.parseBoolean(keyword_val[1]);
                                break;
                            case "maximum":
                            case "distance_delta":
                                distanceDelta = Float.parseFloat(keyword_val[1]);
                                break;
                            case "file_path":
                                break;
                            case "save":
                                saveSCV = true;
                                break;
                            default:
                                logService.error("Keyword " + keyword_val[0] + " not found\nDid you mean: " + getTheClosestMatch(keywords, keyword_val[0]) + "?");
                                return false;
                        }
                    } catch (Exception e) {
                        logService.error("Failed to parse argument:" + keyword_val[0]);
                        return false;
                    }
                } else {
                    logService.error("Malformed argument String. Did you remember to format it as keyword=value and wrap filepaths with quotes? Entire argument string was: " + arg);
                    return false;
                }

            }

        } else {
            InputStream is = sSMLMA.class.getResourceAsStream("/Readme.txt"); // OK
            String content = GetInputStream(is);

            // Declare and fill the UI with all options
            GenericDialogPlus gd = new GenericDialogPlusPlus("settings");
            String[] colors = ownColorTable.getLuts();
            {
                gd.addFileField("CSV input", filePath, 25);

                gd.addCheckbox("Save to CSV?", saveSCV);
                gd.addToSameRow();
                gd.addDirectoryField("CSV output directory", csv_target_dir, 25);

                gd.addMessage("------------------------------------------Angles and Distances-----------------------------------------------------------------------------------------------------------------");

                gd.addMessage("If you want to override the calculated values change the values below.\nWill only overwrite values if they are NOT 0.");


                gd.addNumericField("Start Angle (deg)", angRange[0]);
                gd.addToSameRow();
                gd.addNumericField("End Angle (deg)", angRange[1]);

                gd.addNumericField("Start Distance", distRange[0]);
                gd.addToSameRow();
                gd.addNumericField("End Distance", distRange[1]);

                gd.addNumericField("Number of Orders", orders);

                gd.addCheckbox("Restrict delta z", checkforZ);
                gd.addToSameRow();
                gd.addNumericField("max delta z [nm]", zMargin);
                gd.addToSameRow();
                gd.addMessage("If a Z column is found, this will assume two connected points will not differ in their z value by more than delta z.");

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

                gd.addMessage("------------------------------------------Filtering----------------------------------------------------------------------------------------------------------------------------");

                gd.addCheckbox("Remove Lone Points", toCleanup);
                gd.addToSameRow();
                gd.addNumericField("Required Neighbours", neighbours);
                gd.addToSameRow();
                gd.addNumericField("Required Distance", cleanDistance);
                gd.addMessage("Removes points if there are not at least a number of neighbours in a certain distance.\nWarning: Extremely slow for large datasets");
                gd.addCheckbox("Remove points with high distance delta", checkDistanceOrderDelta);
                gd.addToSameRow();
                gd.addNumericField("Maximum delta", distanceDelta);
                gd.addMessage("Remove points with a higher delta between the distance from the first to second point and the second to third point.");

                gd.addMessage("------------------------------------------Visualisation------------------------------------------------------------------------------------------------------------------------");

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

                gd.addMessage("The graphs created have two new functions under the More... button:\n" +
                        "The 'rescale LUT...' button allows you to rescale the LUT to different distance values or change the LUT altogether\n" +
                        "The 'create Hist from selection' button allows you to create a detailed histogram of all points that lie within the selection\n" +
                        "    this selection can be made using the box selection, freehand or polygon, as g as the area is enclosed.\n" +
                        "These graphs are built on top of of ImageJ's graphs, but feature some untied ends, and thus not all buttons will work " +
                        "on these graphs, and they will occasionally error, but these will not be fatal and can be ignored.");

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

                angInput[0] = (float) (gd.getNextNumber() * (Math.PI / 180f));
                angInput[1] = (float) (gd.getNextNumber() * (Math.PI / 180f));

                distInput[0] = (float) gd.getNextNumber();
                distInput[1] = (float) gd.getNextNumber();

                orders = (int) gd.getNextNumber();

                checkforZ = gd.getNextBoolean();
                zMargin = (float) gd.getNextNumber();

                checkForIntensity = gd.getNextBoolean();
                ratioIntensity = (float) gd.getNextNumber();

                flipAngles = gd.getNextBoolean();
                mirrorAngles = gd.getNextBoolean();

                searchAngle = gd.getNextBoolean();
                deepSearchAngle = gd.getNextBoolean();

                toCleanup = gd.getNextBoolean();
                neighbours = (int) gd.getNextNumber();
                cleanDistance = (float) gd.getNextNumber();

                checkDistanceOrderDelta = gd.getNextBoolean();
                distanceDelta = (float) gd.getNextNumber();

                visualisation = gd.getNextBoolean();
                visualiseZOLA = gd.getNextBoolean();
                binwidth = (float) gd.getNextNumber();

                lutRange[0] = (float) gd.getNextNumber();
                lutRange[1] = (float) gd.getNextNumber();
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
        if(!(visualisation || saveSCV || visualiseZOLA)) {
            logService.error("No output method of any sorts is selected.\nSelect either Visualisation or Save to CSV.");
            return false;
        }

        if(saveSCV && !new File(csv_target_dir).exists()) {
            if(!new File(csv_target_dir).mkdir()) {
                logService.error("Failed to create CSV target directory: " + csv_target_dir);
                return false;
            }
        }

        return true;
    }

    @Override
    public void run() {
        boolean fileError = false;

        // Setup is done here
        // If it returns a failure (false) we execute nothing else
        if (setup()) {
            double csvTime = System.nanoTime(); // Timing loading csv

            List<String> collumns = new ArrayList<>(); // will store the collumns as found in csv
            totalColumns = orders * orderColumns; // total number of collumns in result table

            // Debug stuff for skipping calculation and testing visualisation
            if (debug) {
                //filePath = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\results\\thunderSTORM.csv";
                filePath = "C:\\Users\\Martijn\\Downloads\\Thunderstorm results no drift no filter.csv";
                processing = true;
                visualisation = true;
                visualiseZOLA = false;
                saveSCV = false;
                searchAngle = true;

                //angInput[0] = -0.1f;
                //angInput[1] = 0.1f;
                //distInput[0] = 3500;
                //distInput[1] = 4500;
                try {
                    logService.info("Found " + ownColorTable.getLuts().length + " LUTs");
                    logService.info(Arrays.toString(ownColorTable.getLuts()));
                    ownColorTable.setLut("NCSA PalEdit/royal.lut");
                } catch (Exception e2) {
                    logService.error("Failed to set LUT");
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
                    fileError = true;
                } catch (LapackException e) {
                    logService.info("Lapack error");
                    fileError = true;
                } catch (Exception e) {
                    logService.info("Error reading file. Does the csv start with a header?");
                    fileError = true;
                    e.printStackTrace();
                }



                // Initialise the units and header arrays
                revOptionsIndices = new int[possible_options.length];
                unitsIndices = new int[collumns.size()];

                //Load not_found (-1) as default
                for(int i = 0; i < possible_options.length; i++){
                    revOptionsIndices[i] = -1;
                }

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
                logService.info("Loading CSV took " + String.format("%.3f", csvTime / 1000000000) + " s");

            }

            double processingTime = System.nanoTime();
            final FloatMatrix data; // Holds all data (including datapoints that are not pairs)
            FloatMatrix finalPossibilities; // Holds all pairs

            // This boolean indicates if the angle finding succeeded or failed, along side some other checks
            boolean succes = false;

            if (angInput[0] * angInput[1] * distInput[0] * distInput[1] != 0) {
                searchAngle = false;
                deepSearchAngle = false;
                flipAngles = false;
                mirrorAngles = false;
            }

            if (processing && floatMatrix != null) {
                // Echo back settings used when searching for the best angle
                if (searchAngle) {
                    logService.info("Run: " + runNumber + ". Determining Angle with settings: Flip Angle: " + flipAngles + ", Mirror Angle: " + mirrorAngles);
                }

                // frame, x, y, z, intensity
                // Load the relevant data, if no z is found, load a column of 0's
                if(revOptionsIndices[4] == -1 || revOptionsIndices[5] == -1){
                    data = new FloatMatrix(floatMatrix.rows, 5);
                    data.putColumn(0, floatMatrix.getColumn(revOptionsIndices[1])); // frame
                    data.putColumn(1, floatMatrix.getColumn(revOptionsIndices[2])); // x
                    data.putColumn(2, floatMatrix.getColumn(revOptionsIndices[3])); // y
                    if(revOptionsIndices[4] != -1) {
                        hasZ = true;
                        data.putColumn(3, floatMatrix.getColumn(revOptionsIndices[4])); // z
                    }

                    if(revOptionsIndices[5] != -1) {
                        data.putColumn(4, floatMatrix.getColumn(revOptionsIndices[5])); // intensity
                        hasIntensity = true;
                    }
                    else {
                        fillCollumn(data, 4, 1);
                    }

                } else {
                    hasZ = true;
                    hasIntensity = true;
                    data = floatMatrix.getColumns(new int[]{revOptionsIndices[1], revOptionsIndices[2], revOptionsIndices[3], revOptionsIndices[4],revOptionsIndices[5]});
                }

                // If any var is not set, we have to calculate them all
                if (angInput[0] * angInput[1] * distInput[0] * distInput[1] == 0) {

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
                    succes = angleAnalyzer.getSuccess();
                    if (distRange[0] > distRange[1])  succes = false;

                } else {
                    // if all values are filled in this is the only check we need to do
                    succes = !(distRange[0] > distRange[1]);

                    angRange[0] = angInput[0];
                    angRange[1] = angInput[1];
                    distRange[0] = distInput[0];
                    distRange[1] = distInput[1];

                }
            } else {
                // If we dont process (during debug) we directly load the data
                data = floatMatrix;
            }

            // If we failed (and are processing) we report the error and stop execution
            if ((!succes && processing) | floatMatrix == null) {
                if (distRange[0] > distRange[1]) {
                    logService.error("The distance had to be positive: " + distRange[0] + " is larger than " + distRange[1]);
                } else if(fileError){
                    logService.error("Could not load the file. Is the path (" + filePath + ") correct?");
                } else {
                    logService.error("No features were detected. Are there pairs in this sample?");
                }

            } else {
                // Declaring a few vars due to scope juggling
                FloatMatrix halfOrderMatrix;
                FloatMatrix allOrdersCombined;
                if (processing) {

                    final FloatMatrix frameCollumn = data.getColumn(0);
                    final int[] frames = getFrameNumbers(data, 0);
                    final int numFrames = frames.length;

                    // Echo back amount of frames and points
                    logService.info("Total Frames: " + numFrames);
                    logService.info("Total Points: " + floatMatrix.rows);

                    final AtomicInteger ai = new AtomicInteger(0); //Atomic Integer is a thread safe incremental integer
                    final AtomicBoolean reportOrders = new AtomicBoolean(true); //To only report more orders once
                    final Thread[] threads = createThreadArray(coreCount); //Get array of threads

                    final FloatMatrix[] intermediateFinals = new FloatMatrix[coreCount]; // All intermediate results to be merged later
                    // prepare final verions of variables
                    final boolean finalIntensityCheck = checkForIntensity;
                    final float finalRatioIntensity = ratioIntensity;
                    final boolean finalcheckforZ = checkforZ;
                    final float finalzMargin = zMargin;
                    final boolean finalhasZ = hasZ;
                    //Set the run function for each thread
                    for (int ithread = 0; ithread < threads.length; ithread++) {

                        final int finalIthread = ithread;
                        threads[finalIthread] = new Thread(() -> {

                            // Will hold the final values for this thread
                            intermediateFinals[finalIthread] = new FloatMatrix(0, totalColumns);

                            // Process each frame
                            for (int frameIndex = ai.getAndIncrement(); frameIndex < numFrames; frameIndex = ai.getAndIncrement()) {
                                final int frame = frames[frameIndex];
                                // Showing process to the user
                                if (runningFromIDE && frame % 1000 == 0) logService.info("\r" + frame + "/" + numFrames);
                                IJ.showProgress(frame, numFrames);
                                IJ.showStatus(frame + "/" + numFrames);

                                final int[] frameIndicices = frameCollumn.eq(frame).findIndices(); // All rows indices for current frame
                                final FloatMatrix frameData = data.getRows(frameIndicices); //All rows for current frame


                                // Distance Matrices that show distance from one point to each other point in X and Y
                                final FloatMatrix subtractedX = makeSubstractedMatrix(frameData.getColumn(1));
                                final FloatMatrix subtractedY = makeSubstractedMatrix(frameData.getColumn(2));

                                // Matrices that have distances and angles from any point to another point in this frame
                                final FloatMatrix distances = Distance(subtractedX, subtractedY);
                                final FloatMatrix angles = atan2(subtractedX, subtractedY, distances);

                                // Filter to only points that match the distance and angle restrictions;
                                // The value indicates both the start and end point
                                final int[] correctAngleAndDistance;

                                if(finalhasZ && finalcheckforZ){
                                    final FloatMatrix zDistances = makeSubstractedMatrix(frameData.getColumn(3));
                                    correctAngleAndDistance = distances.gt(distRange[0]).and(distances.lt(distRange[1]))
                                            .and(angles.gt(angRange[0])).and(angles.lt(angRange[1])).and(zDistances.lt(finalzMargin)).findIndices();
                                } else {
                                    correctAngleAndDistance = distances.gt(distRange[0]).and(distances.lt(distRange[1]))
                                            .and(angles.gt(angRange[0])).and(angles.lt(angRange[1])).findIndices();
                                }

                                // If points were found, we must process them
                                if (correctAngleAndDistance.length > 1) {

                                    // Matrices to hold the possibilities and combined ones
                                    final FloatMatrix possibilities = new FloatMatrix(correctAngleAndDistance.length, totalColumns);
                                    FloatMatrix intermediateFinalPossibilities = new FloatMatrix(0, totalColumns);

                                    // For each pair, record it to out temporary matrix
                                    for (int i = 0; i < correctAngleAndDistance.length; i++) {
                                        int index = correctAngleAndDistance[i];

                                        // Check if the intensity ratio checks out, if that check is enabled
                                        if (finalIntensityCheck &&
                                                (frameData.get(index / distances.rows, 4) /
                                                        frameData.get(index % distances.rows, 4)) > finalRatioIntensity
                                        ) continue;

                                        // Put all info about the two points into our matrix
                                        possibilities.putRow(i, extend(new FloatMatrix(1, orderColumns * 2,
                                                0,                                                     //0 (global index goes here later)
                                                frame,                                                          //1 frame
                                                1 + Math.floorDiv(index, distances.rows),                       //2 start index in frame
                                                frameData.get(index / distances.rows, 1),    //3 start x
                                                frameData.get(index / distances.rows, 2),    //4 start y
                                                frameData.get(index / distances.rows, 3),    //5 start z
                                                frameData.get(index / distances.rows, 4),    //6 start intensity
                                                1 + Math.floorMod(index, distances.rows),                       //7 end index in frame
                                                frameData.get(index % distances.rows, 1),    //8 end x
                                                frameData.get(index % distances.rows, 2),    //9 end y
                                                frameData.get(index % distances.rows, 3),    //10 end z
                                                frameData.get(index % distances.rows, 4),    //11 end intensity
                                                distances.get(index),                                           //12 distance
                                                angles.get(index)                                               //13 angle
                                        ), 1, totalColumns));
                                    }

                                    FloatMatrix connectedIds = possibilities.getColumn(7); // all frame indices to what a point is connected
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
                                                float cumZ = possibilities.get(i, 5);
                                                float cumInt = possibilities.get(i, 6);

                                                for (int index2 : identicalIds.findIndices()) {
                                                    if (index2 != i) {
                                                        cumX += possibilities.get(index2, 3);
                                                        cumY += possibilities.get(index2, 4);
                                                        cumZ += possibilities.get(index2, 5);
                                                        cumInt += possibilities.get(index2, 6);
                                                        toSkip.add(index2);
                                                    }
                                                }

                                                FloatMatrix temp_arr = possibilities.getRow(i).dup();
                                                temp_arr.put(3, cumX / identicalIdSum);
                                                temp_arr.put(4, cumY / identicalIdSum);
                                                temp_arr.put(5, cumZ / identicalIdSum);
                                                temp_arr.put(6, cumInt / identicalIdSum);

                                                intermediateFinalPossibilities = FloatMatrix.concatVertically(intermediateFinalPossibilities, temp_arr);

                                            } else {
                                                intermediateFinalPossibilities = FloatMatrix.concatVertically(intermediateFinalPossibilities, possibilities.getRow(i));
                                            }
                                        }
                                    }
                                    // Add the points calculated in this frame to this threads buffer
                                    // We also connect the orders here: (1-2, 2-3, 3-4 -> 1-2-3-4)
                                    intermediateFinals[finalIthread] = FloatMatrix.concatVertically(intermediateFinals[finalIthread], connectOrders(intermediateFinalPossibilities, orders, orderColumns, reportOrders));

                                } else if (correctAngleAndDistance.length == 1) { // We found only a single point here, so we just add it and do no connecting
                                    int index = correctAngleAndDistance[0];
                                    intermediateFinals[finalIthread] = FloatMatrix.concatVertically(intermediateFinals[finalIthread], extend(new FloatMatrix(1, orderColumns * 2,
                                            0,
                                            frame,
                                            1 + Math.floorDiv(index, distances.rows),
                                            frameData.get(index / distances.rows, 1),
                                            frameData.get(index / distances.rows, 2),
                                            frameData.get(index / distances.rows, 3),
                                            frameData.get(index / distances.rows, 4),
                                            1 + Math.floorMod(index, distances.rows),
                                            frameData.get(index % distances.rows, 1),
                                            frameData.get(index % distances.rows, 2),
                                            frameData.get(index % distances.rows, 3),
                                            frameData.get(index % distances.rows, 4),
                                            distances.get(index),
                                            angles.get(index)
                                    ), 1, totalColumns));
                                }
                            }
                        }); //end of thread creation
                    }

                    // Start the threads
                    // This actually starts the calculations
                    startAndJoin(threads);

                    //Combine all Matrices together
                    finalPossibilities = new FloatMatrix(0, totalColumns);
                    for (FloatMatrix floatMatrix1 : intermediateFinals) {
                        finalPossibilities = FloatMatrix.concatVertically(finalPossibilities, floatMatrix1);
                    }


                    // Echo back time it took
                    processingTime = System.nanoTime() - processingTime;
                    logService.info("Processing data took " + String.format("%.3f", processingTime / 1000000000) + " s");

                    // Ensure nothing went wrong and echo back how many points we found
                    // Also clean up some garbage since we are done processing and there are many things we no longer need
                    logService.info("Pairs in the 0th-1st order found: " + finalPossibilities.rows);
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
                    try{
                        boolean[][] checks = checkForRetry(finalPossibilities.getColumn(13));

                        if (sum(checks[0]) < 2) {
                            logService.info("Nothing resembling a guassian was found for the angles.");
                            retry = true;
                        }


                        if (checks[1][0]) {
                            if (checks[1][1])
                                logService.info("The left tail of the angle histogram seems to be partially cut-off");
                            else
                                logService.info("The right tail of the angle histogram seems to be partially cut-off");
                            logService.info("You can take the angle value given and adjust it manually.");
                        }
                    } catch (Exception e){
                        logService.error("Failed to analyze angles. No features assumed");
                        retry = true;
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

                            angRange[0] = finalPossibilities.getColumn(13).min();
                            angRange[1] = finalPossibilities.getColumn(13).max();

                            distRange[0] = finalPossibilities.getColumn(12).min();
                            distRange[1] = finalPossibilities.getColumn(12).max();

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

                            if (runningFromIDE) logService.info(message);
                            if (!runningFromMacro)IJ.showMessage(message);
                        }

                    }

                    // Sort by frame id
                    finalPossibilities = sort(finalPossibilities, 1);

                    // Add global index for each row
                    int id = 0;
                    for (int i = 0; i < finalPossibilities.rows; i++) {
                        finalPossibilities.put(i, 0, id++);
                    }


                    // Anything after this point is skipped if we are not in the final run
                    // So this point is only reached with the best (hopefully) results
                    if ((!(retry && searchAngle) | (foundBestResult && deepSearchAngle)) && displayInfo && finalPossibilities.rows > 0) {
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
                            if(finalPossibilities.rows == 0) {
                                logService.info("Filtering resulted in no points. Restoring old data");
                                IJ.showMessage("Filtering resulted in no points. Restoring old data");
                                finalPossibilities = backup;
                            }
                        }

                        // Apply checks for distance delta
                        if(checkDistanceOrderDelta){
                            FloatMatrix deltaFirstSecondDistance = finalPossibilities.getColumn(12).sub(finalPossibilities.getColumn(12 + orderColumns));
                            FloatMatrix comparisonFirstSecondDistance = deltaFirstSecondDistance.get(finalPossibilities.getColumn(12 + orderColumns).ne(0.0f).findIndices());
                            ImagePlus comparisonDistanceImagePlus = new ImagePlus("", new FloatProcessor(comparisonFirstSecondDistance.toArray2()));
                            HistogramWindow comparisonDistanceHist = new HistogramWindow("delta distance (no filter)", comparisonDistanceImagePlus, getBins(comparisonFirstSecondDistance, binwidth), comparisonFirstSecondDistance.min(), comparisonFirstSecondDistance.max());

                            if (runningFromIDE) comparisonDistanceHist.getImagePlus().show();


                            int[] distanceDeltaFilterIndices = deltaFirstSecondDistance.le(distanceDelta * 0.5f).and(deltaFirstSecondDistance.ge(distanceDelta * -0.5f)).or(finalPossibilities.getColumn(12 + orderColumns).eq(0.0f)).findIndices();
                            logService.info("Removed " + (finalPossibilities.rows - distanceDeltaFilterIndices.length) + " points due to too high difference between 0-1 distance and 1-2 distance.");
                            finalPossibilities = finalPossibilities.getRows(distanceDeltaFilterIndices);
                        }

                        // We have our final data
                        // some other versions could give improved information by combining position data from multiple orders into one
                        // We do this into two ways here
                        // Once with a only one pair (0th-1st) order
                        // Another one with all orders combined, as many as there are for one row or data
                        // It is a lot of code doing some simple combining
                        // 0 id, 1 frame, 2 x, 3 y, 4 z, 5 intensity, 6 distance, 7 angle
                        halfOrderMatrix = new FloatMatrix(finalPossibilities.rows, orderColumns + 1);
                        allOrdersCombined = new FloatMatrix(finalPossibilities.rows, orderColumns + 1);
                        {
                            halfOrderMatrix.putColumn(0, finalPossibilities.getColumn(0)); //id
                            halfOrderMatrix.putColumn(1, finalPossibilities.getColumn(1)); //frame
                            halfOrderMatrix.putColumn(2, finalPossibilities.getColumn(3).add(finalPossibilities.getColumn(8)).divi(2.0f)); //x
                            halfOrderMatrix.putColumn(3, finalPossibilities.getColumn(4).add(finalPossibilities.getColumn(9)).divi(2.0f)); //y
                            if(hasZ) halfOrderMatrix.putColumn(4, finalPossibilities.getColumn(5).add(finalPossibilities.getColumn(10)).divi(2.0f)); // z
                            halfOrderMatrix.putColumn(5, finalPossibilities.getColumn(6)); //intensity
                            halfOrderMatrix.putColumn(6, finalPossibilities.getColumn(12)); //distance
                            halfOrderMatrix.putColumn(7, finalPossibilities.getColumn(13)); //angle

                            allOrdersCombined.putColumn(0, finalPossibilities.getColumn(0)); //id
                            allOrdersCombined.putColumn(1, finalPossibilities.getColumn(1)); //frame
                            allOrdersCombined.putColumn(5, finalPossibilities.getColumn(6)); //intensity


                            FloatMatrix offsets = new FloatMatrix(finalPossibilities.rows, 5);
                            offsets.putColumn(0, finalPossibilities.getColumn(8).sub(finalPossibilities.getColumn(3)).divi(2.0f)); // x
                            offsets.putColumn(1, finalPossibilities.getColumn(9).sub(finalPossibilities.getColumn(4)).divi(2.0f)); // y
                            offsets.putColumn(2, finalPossibilities.getColumn(10).sub(finalPossibilities.getColumn(5)).divi(2.0f)); // z
                            offsets.putColumn(3, finalPossibilities.getColumn(12)); // distance
                            offsets.putColumn(4, finalPossibilities.getColumn(13)); // angle

                            for (int i = 0; i < orders - 2; i++) {
                                FloatMatrix relevantRows = finalPossibilities.getColumn(14 + (i * orderColumns)).ne(0.0f);
                                offsets.getColumn(0).addi(finalPossibilities.getColumn(15 + (i * orderColumns)).subi(finalPossibilities.getColumn(8 + (i * orderColumns))).divi(2.0f).muli(relevantRows));
                                offsets.getColumn(1).addi(finalPossibilities.getColumn(16 + (i * orderColumns)).subi(finalPossibilities.getColumn(9 + (i * orderColumns))).divi(2.0f).muli(relevantRows));
                                offsets.getColumn(2).addi(finalPossibilities.getColumn(17 + (i * orderColumns)).subi(finalPossibilities.getColumn(10 + (i * orderColumns))).divi(2.0f).muli(relevantRows));
                                offsets.getColumn(3).addi(finalPossibilities.getColumn(19 + (i * orderColumns)));
                                offsets.getColumn(4).addi(finalPossibilities.getColumn(20 + (i * orderColumns)));
                            }


                            allOrdersCombined.putColumn(2, finalPossibilities.getColumn(3));
                            allOrdersCombined.putColumn(3, finalPossibilities.getColumn(4));
                            allOrdersCombined.putColumn(4, finalPossibilities.getColumn(5));
                            allOrdersCombined.putColumn(6, offsets.getColumn(3));
                            allOrdersCombined.putColumn(7, offsets.getColumn(4));



                            FloatMatrix relevantRows = finalPossibilities.getColumn((orders - 1) * orderColumns).ne(0.0f);
                            for (int i = orders; i > 1; i--) {
                                allOrdersCombined.getColumn(2).addi(offsets.getColumn(0).divi((float) i).muli(relevantRows)); // x
                                allOrdersCombined.getColumn(3).addi(offsets.getColumn(1).divi((float) i).muli(relevantRows)); // y
                                allOrdersCombined.getColumn(4).addi(offsets.getColumn(2).divi((float) i).muli(relevantRows)); // z
                                allOrdersCombined.getColumn(6).addi(offsets.getColumn(3).divi((float) i).muli(relevantRows)); // distance
                                allOrdersCombined.getColumn(7).addi(offsets.getColumn(4).divi((float) i).muli(relevantRows)); // angle

                                relevantRows.xori(finalPossibilities.getColumn((i - 2) * orderColumns).ne(0.0f));
                            }
                        }
                        if (visualisation) {

                            HistogramWindow[] histograms = new HistogramWindow[orders - 1]; // We create some histograms for each distance order we want to visualise

                            // For each order we fill the distance into each matrix
                            // this allows one to compare the distances for each order
                            // they should be the same for each order to order, but its still useful
                            for (int i = 0; i < orders - 1; i++) {
                                FloatMatrix relevantData;

                                // Visualise either every row or only the rows that have data in it
                                // Since not each row might have data for each order
                                if (i == 0) {
                                    relevantData = finalPossibilities.getColumn(12);
                                } else {
                                    relevantData = finalPossibilities.getColumn(12 + (i * orderColumns));
                                    relevantData = relevantData.get(relevantData.ne(0.0f).findIndices());
                                }

                                logService.info("Found " + relevantData.rows + " connections in the " + getTitleHist(i) + " order");

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
                            FloatMatrix AnglesRad = finalPossibilities.getColumn(13);
                            AnglesRad.muli((float) (180/Math.PI));
                            ImagePlus Angles = new ImagePlus("", new FloatProcessor(AnglesRad.toArray2()));
                            HistogramWindow angleHist = new HistogramWindow("Angles (degrees)", Angles, getBins(AnglesRad, (float) (0.005f * (180/Math.PI))), angRange[0] * (180/Math.PI), angRange[1]* (180/Math.PI));
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


                            CustomPlot distancePlot = new CustomPlot("Distance", "x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]", "y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]", lutService, defaultLUT);

                            distancePlot.add(shapes[0], toDouble(halfOrderMatrix.getColumn(2)), toDouble(halfOrderMatrix.getColumn(3)), toDouble(halfOrderMatrix.getColumn(6)));

                            distancePlot.setLimitsToFit(true); // Ensure all points are visible
                            distancePlot.setLUTLegend("Distance", 512, distRange[0], distRange[1]); // Add the LUT as a legend
                            distancePlot.setLimits(Float.NaN, Float.NaN, Float.NaN, Float.NaN); // Ensure all points are visible, again
                            distancePlot.show();


                            /* This can be used to demonstrate the positions the first two points will have
                            {
                                CustomPlot distancePlot2 = new CustomPlot("Distance2", "x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]", "y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]", lutService, defaultLUT);

                                distancePlot2.add(shapes[0], toDouble(finalPossibilities.getColumn(3)), toDouble(finalPossibilities.getColumn(4)), toDouble(finalPossibilities.getColumn(12)));
                                distancePlot2.add(shapes[0], toDouble(finalPossibilities.getColumn(8)), toDouble(finalPossibilities.getColumn(9)), toDouble(finalPossibilities.getColumn(12)));
                                distancePlot2.setLimitsToFit(true); // Ensure all points are visible
                                distancePlot2.setLUTLegend("Distance", 512, distRange[0], distRange[1]); // Add the LUT as a legend
                                distancePlot2.setLimits(Float.NaN, Float.NaN, Float.NaN, Float.NaN); // Ensure all points are visible, again
                                distancePlot2.show();
                            }
                             */

                        }

                        // Removed at runtime
                        if(debugFlag) {
                            float resolutionBinWidth = 250; // 5nm


                            FloatMatrix localMedianMatrix = finalPossibilities.getColumns(new int[]{3, 4});
                            localMedianMatrix.divi(resolutionBinWidth);

                            int[] xValues = localMedianMatrix.getColumn(0).toIntArray();
                            int[] yValues = localMedianMatrix.getColumn(1).toIntArray();
                            float[] distanceValues = finalPossibilities.getColumn(12).toArray();

                            IntSummaryStatistics Xstat = Arrays.stream(xValues).summaryStatistics();
                            IntSummaryStatistics Ystat = Arrays.stream(xValues).summaryStatistics();

                            ArrayList<Point3D> pointList = new ArrayList<>();

                            for (int i = 0; i < localMedianMatrix.rows; i++) {
                                pointList.add(new Point3D(xValues[i], yValues[i], distanceValues[i]));
                            }

                            ArrayList<Float> localMedianList = new ArrayList<>();
                            ArrayList<Double> localAverageList = new ArrayList<>();
                            ArrayList<Integer> localModeList = new ArrayList<>();

                            for (int x = Xstat.getMin(); x <= Xstat.getMax(); x++) {
                                final int finalX = x;

                                List<Point3D> sublist = pointList.stream()
                                        .filter(p -> p.getX() == finalX).collect(Collectors.toList());

                                if(x%100 == 0) System.out.println(x);

                                for (int y = Ystat.getMin(); y <= Ystat.getMax(); y++) {

                                    final int finalY = y;

                                    List<Float> currentMedianList = sublist.stream()
                                            .filter(p -> p.getY() == finalY).map(Point3D::getZ).collect(Collectors.toList());


                                    localModeList.add(mode(currentMedianList));
                                    double average = getAverage(currentMedianList);
                                    if(average != 0.0) localAverageList.add(average);

                                    float currentMedian = getMedian(currentMedianList);
                                    if (currentMedian != Float.MIN_VALUE)
                                        localMedianList.add(currentMedian);

                                    if(currentMedianList.size() > 200 && Random.nextDouble() < 0.1 && getAverage(currentMedianList) > 1800) {
                                        // float[] minMax = getFloatMinMax(currentMedianList);
                                        createHist(toFloat(currentMedianList.toArray(new Float[0])), getFBins(currentMedianList, 2), 1900,  2100, x + ", " + y, runningFromIDE);
                                    }
                                }

                                pointList.removeIf(p -> p.getX() == finalX);
                            }


                            float[] medianMinMax = getFloatMinMax(localMedianList);


                            //createHist(toFloat(localMedianList.toArray(new Float[0])), getFBins(localMedianList, 2), 1500,  1800, "medians1", runningFromIDE);
                            createHist(toFloat(localMedianList.toArray(new Float[0])), getFBins(localMedianList, 2), medianMinMax[0],  medianMinMax[1], "medians2", runningFromIDE);

                            float[] averageMinMax = getDoubleMinMax(localAverageList);
                            //createHist(toFloat(localAverageList.toArray(new Double[0])), getDBins(localAverageList, 2), 1500,  1800, "mean1", runningFromIDE);
                            createHist(toFloat(localAverageList.toArray(new Double[0])), getDBins(localAverageList, 2), averageMinMax[0],  averageMinMax[1], "mean2", runningFromIDE);

                            double[] modes = toDouble(localModeList);
                            //int[] ranges = new int[]{1500, 1800};
                            int[] ranges = new int[]{(int) averageMinMax[0], (int) averageMinMax[1]};
                            createHist(toFloat(modes), (int) ((ranges[1]-ranges[0])/2.5), ranges[0], ranges[1], "modes1", runningFromIDE);
                            //ranges = new int[]{1900, 2100};
                            //createHist(toFloat(modes), (int) ((ranges[1]-ranges[0])/2.5), ranges[0], ranges[1], "modes2", runningFromIDE);

                        }

                        if (saveSCV) {

                            // Remove unneeded Z column before saving
                            //System.out.println(finalPossibilities.getRow(0));
                            logService.info("Writing files to " + csv_target_dir);
                            FloatMatrix noZMatrix = null;
                            if(!hasZ){
                                noZMatrix = new FloatMatrix(finalPossibilities.rows, orders * (orderColumns - 1));

                                // No range copy, so this is slow
                                noZMatrix.putColumn(0, finalPossibilities.getColumn(0)); // id
                                noZMatrix.putColumn(1, finalPossibilities.getColumn(1)); // frame
                                noZMatrix.putColumn(2, finalPossibilities.getColumn(2)); // id in frame
                                noZMatrix.putColumn(3, finalPossibilities.getColumn(3)); // x
                                noZMatrix.putColumn(4, finalPossibilities.getColumn(4)); // y
                                noZMatrix.putColumn(5, finalPossibilities.getColumn(6)); // intensity

                                //i.e. 14: id, 15: x, 16: y, 17: z, 18: intensity, 19: distance, 20: angle
                                for(int i = 1; i < orders; i++){
                                    noZMatrix.putColumn((orderColumns - 1) * i, finalPossibilities.getColumn(orderColumns * i)); // id in frame
                                    noZMatrix.putColumn((orderColumns - 1) * i + 1, finalPossibilities.getColumn(orderColumns * i + 1)); // x
                                    noZMatrix.putColumn((orderColumns - 1) * i + 2, finalPossibilities.getColumn(orderColumns * i + 2)); // y
                                    noZMatrix.putColumn((orderColumns - 1) * i + 3, finalPossibilities.getColumn(orderColumns * i + 4)); // intensity
                                    noZMatrix.putColumn((orderColumns - 1) * i + 4, finalPossibilities.getColumn(orderColumns * i + 5)); // distance
                                    noZMatrix.putColumn((orderColumns - 1) * i + 5, finalPossibilities.getColumn(orderColumns * i + 6)); // angle
                                }
                                //System.out.println(noZMatrix.getRow(0));
                            }




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
                            String distanceUnit = (unit_prefixes[unitsIndices[revOptionsIndices[2]]].equals(unit_prefixes[unitsIndices[revOptionsIndices[3]]]) ? unit_prefixes[unitsIndices[revOptionsIndices[2]]] : ("(" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "*" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + ")^0.5"));
                            // Add the headers for each order
                            for (int i = 0; i <= orders; i++) {
                                LongHeader.add("index " + i);
                                LongHeader.add("x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "] " + i);
                                LongHeader.add("y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "] " + i);
                                if(hasZ) LongHeader.add("z [" + unit_prefixes[unitsIndices[revOptionsIndices[4]]] + "] " + i);
                                if(hasIntensity) LongHeader.add("intensity [" + unit_prefixes[unitsIndices[revOptionsIndices[5]]] + "] " + i);
                                else LongHeader.add("intensity [photons] " + i);
                                if (i > 0) {
                                    LongHeader.add((i - 1) + "-" + i + " distance [" + distanceUnit + "]");
                                    LongHeader.add((i - 1) + "-" + i + "angle");
                                }
                            }
                            // Only add one header level to the short ones
                            ShortHeader.add("x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]");
                            ShortHeader.add("y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]");
                            if(hasZ) ShortHeader.add("z [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]");
                            if(hasIntensity) ShortHeader.add("intensity [" + unit_prefixes[unitsIndices[revOptionsIndices[5]]] + "]");
                            else ShortHeader.add("intensity [photons] ");
                            ShortHeader.add("distance [" + distanceUnit + "]");
                            ShortHeader.add("angle");

                            try {
                                logService.info(("Writing all_orders.csv"));
                                // Save all data using the proper header, including one that easily is loaded into ThunderSTORM again for visualisation etc
                                if (hasZ)
                                    SaveCSV(finalPossibilities, LongHeader, Paths.get(csv_target_dir, "all_orders.csv"));
                                else
                                    SaveCSV(noZMatrix, LongHeader, Paths.get(csv_target_dir, "all_orders.csv"));
                            } catch(Exception e){
                                logService.error("Could not create file: all_orders.csv. Is the file opened anywhere?");
                            }

                            try {
                                logService.info(("Writing two_orders_combined_positions.csv"));
                                SaveCSV(halfOrderMatrix, ShortHeader, Paths.get(csv_target_dir, "two_orders_combined_positions.csv"));
                            } catch(Exception e){
                                logService.error("Could not create file: two_orders_combined_positions.csv. Is the file opened anywhere?");
                            }
                            try {
                                logService.info(("Writing all_orders_combined_positions.csv"));
                                SaveCSV(allOrdersCombined, ShortHeader, Paths.get(csv_target_dir, "all_orders_combined_positions.csv"));
                            } catch(Exception e){
                                logService.error("Could not create file: all_orders_combined_positions.csv. Is the file opened anywhere?");
                            }
                            try {
                                logService.info(("Writing thunderSTORM.csv"));
                                saveThunderSTORM(Paths.get(csv_target_dir, "thunderSTORM.csv"), finalPossibilities.getColumns(new int[]{0, 1, 3, 4, 6, 12}));
                            } catch(Exception e){
                                logService.error("Could not create file: thunderSTORM.csv. Is the file opened anywhere?");
                            }
                            logService.info(("Finished writing all csv files."));
                        }

                        if(visualiseZOLA) {
                            logService.info("ZOLA Visualisation");

                            try {
                                Path tmpfile; // We need a tmpfile if no saving is done

                                if (runningFromIDE) { // Doesn;t work from IDE because of the isolated ImageJ isntance
                                    logService.info("Running from IDE does not work for ZOLA integration");
                                } else {
                                    // If we already saved, great, otherwise save a tmp thunderstorm file to use
                                    if (saveSCV) {
                                        tmpfile = Paths.get(csv_target_dir, "thunderSTORM.csv");
                                    } else {
                                        tmpfile = Paths.get(IJ.getDirectory("temp"), "tmp.csv");
                                        saveThunderSTORM(tmpfile, finalPossibilities.getColumns(new int[]{0, 1, 3, 4, 6, 12}));
                                    }

                                    Prefs.set("Zola.showLUT", true); // Show lut on image
                                    Prefs.set("Zola.pathlocalization", tmpfile.toString()); // Load the file
                                    Prefs.set("Zola.is3Drendering", true); // Set 3D rendering

                                    IJ.run("Import table"); // Import our table (also shows 2D histogram for some reason)
                                    IJ.run("2D/3D histogram"); // Show the 3D histogram (does not work from macro, and shows 2D)
                                }
                            }  catch (Exception e) {
                                logService.info("ZOLA integration failed");
                                e.printStackTrace();
                            }
                        }

                        WriteXML(Paths.get(csv_target_dir, "info.xml"));
                        logService.info("Finished writing XML file");

                    }
                }
            }
        }
    }


    // Only run from the IDE
    public static void main(String[] args) {
        runningFromIDE = true; //this is really dumb

        /*
        "csv_in", "csv_out","angle_start", "angle_end",
        "distance_start", "distance_end",
        "order_number", "check_order_intensity", "check_order_ratio",
        "angle_flip", "angle_mirror", "angle_search", "angle_deep_search",
        "lone_pair_remove", "lone_pair_neighbours", "lone_pair_distance",
        "visualisation", "visualisationZOLA", "hist_binwidth", "LUT", "LUT_start", "LUT_end",
        "check_z", "check_z_margin", "check_distance_delta", "distance_delta"
        */

        // private String csv_target_dir = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\results";

        // private String filePath = "F:\\ThesisData\\output\\output3_drift.csv";
        // private String filePath = "F:\\ThesisData\\output\\combined_drift.csv";
        // private float[] angRange = new float[] {(float) (-0.03 * Math.PI), (float) (0.07 * Math.PI)};
        // private float[] distRange = new float[] {1500, 2200};

        // private String filePath = "F:\\ThesisData\\output\\4_grating_drift.csv";
        // private final float[] angRange = {(float) (-1 * Math.PI), (float) (-0.95 * Math.PI) }; //more than and less than
        // private final float[] distRange = {1940, 2600}; //1800 3000 (1940, 2240)

        //debug_arg_string = "angle_start=-0.5 angle_end=4 distance_start=3000 distance_end=3500 csv_in='\\\\WURNET.NL\\Homes\\gobes001\\AppData\\FolderRedirection\\Desktop\\PhD\\Projects\\SpectralSMLM\\Koen Martens Work\\RawData\\Lorenzo_pSMLM_wavelet15.csv' visualisation=true";
        //debug_arg_string = "angle_start=0 angle_end=4 distance_start=2250 distance_end=3000 csv_in='\\\\WURNET.NL\\Homes\\gobes001\\AppData\\FolderRedirection\\Desktop\\PhD\\Projects\\SpectralSMLM\\Koen Martens Work\\RawData\\SMAP_1B_20200204_combined.csv' visualisation=true";
        //debug_arg_string = "csv_in='C:\\\\Users\\\\gobes001\\\\PhD\\\\Projects\\\\SpectralSMLM\\\\Koen Martens Work\\\\RawData\\\\20200120 part Fig3\\\\Movie1_pSMLM3_wavelet15.csv' csv_out='C:\\\\Users\\\\gobes001\\\\PhD\\\\Projects\\\\SpectralSMLM\\\\Koen Martens Work\\\\RawData\\\\20200120 part Fig3\\\\Movie1_15_2' angle_start=-0.7 angle_end=4 distance_start=2500 distance_end=3500 visualisation=true lone_pair_remove=true lone_pair_neighbours=250 lone_pair_distance=400";
        debug_arg_string = "distance_start=1520 distance_end=2200 angle_start=-0.65 angle_end=0.5 csv_out='I:\\ThesisData\\input\\testouput\\' csv_in='I:\\ThesisData\\input\\localisations_driftcorr_first1000.csv' visualisation=true";

        // v1: 35.759 v_all: 70.142
        //debug_arg_string = "";
        net.imagej.ImageJ ij = new ImageJ();
        ij.ui().showUI();

        ij.command().run(sSMLMA.class, true);
    }
}

class Point3D {
    float x, y, z;

    Point3D(float x, float y, float z){
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public float getX() { return x;}
    public float getY() { return y;}
    public float getZ() { return z;}
}

