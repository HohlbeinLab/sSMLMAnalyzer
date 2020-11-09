package com.wurgobes.sSMLMAnalyzer;



import ij.*;

import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.gui.HistogramWindow;
import ij.process.FloatProcessor;


import org.scijava.plugin.Parameter;

import java.io.IOException;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import fiji.util.gui.GenericDialogPlus;

import org.jblas.FloatMatrix;
import org.jblas.exceptions.LapackException;

import static com.wurgobes.sSMLMAnalyzer.Util.*;
import static com.wurgobes.sSMLMAnalyzer.levenshtein.getTheClosestMatch;


public class sSMLMA implements PlugIn {


    private final int[] unit_decades = {0, 0, 0, -2, -3, -6, -9, -10, -12, -15};

    private final String[] possible_options = {"id", "frame", "x", "y", "z", "intensity", "offset", "bkgstd", "sigma1", "sigma2", "uncertainty", "detections", "chi"};
    private final String[] unit_prefixes = {"null", "photons", "m", "cm", "mm", "um", "nm", "ang", "pm", "fm"};

    @Parameter
    private String filePath = "F:\\ThesisData\\output\\output3_drift.csv";

    @Parameter
    private String CSV_FILE_NAME = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\small.csv";

    @Parameter
    private final boolean saveSCV = true;

    @Parameter
    private final float[] angRange = {(float) (-0.04 * Math.PI), (float)(0.04 * Math.PI)};

    @Parameter
    private final float[] distRange = {1500, 2500}; //default was {2500, 4500}

    @Parameter
    private final int orders = 4; //This its the number of orders, including the 0th, so 3 would be 0th + 1st + 2nd

    public void setup(){
        GenericDialogPlus gd = new GenericDialogPlus("settings");

        gd.addFileField("CSV input", filePath, 50);
        gd.addFileField("CSV output", CSV_FILE_NAME, 50);

        gd.showDialog();

        filePath = gd.getNextString();
        CSV_FILE_NAME = gd.getNextString();
    }

    @Override
    public void run(String arg) {

        //setup();

        double csvTime = System.nanoTime();

        //Primitive64Matrix.Factory matrixFactory = Primitive64Matrix.FACTORY;

        //LineSplittingParser lineSplittingParser = new LineSplittingParser();

        //lineSplittingParser.parse(filePath, false, );

        FloatMatrix floatMatrix = null;
        List<String> collumns = new ArrayList<>();

        OwnFloatMatrix ownFloatMatrix = new OwnFloatMatrix();

        try {
            floatMatrix =  ownFloatMatrix.loadCSVFile(filePath);
            collumns = ownFloatMatrix.collumns;
        } catch (IOException e) {
            System.out.println("File not found.");
        } catch (LapackException e) {
            e.printStackTrace();
        } catch (Exception e){
            System.out.println("Does the csv start with a header?");
        }

        assert !collumns.isEmpty();
        assert floatMatrix != null;



        int[] revOptionsIndices = new int[possible_options.length];
        int[] unitsIndices = new int[collumns.size()];


        Pattern pattern = Pattern.compile("(\\w+)( [ (\\[](\\w+)[)\\] ])?");

        for(int i = 0; i < collumns.size(); i ++){
            String header = collumns.get(i);
            Matcher matcher = pattern.matcher(header);
            if(matcher.find()){
                revOptionsIndices[getTheClosestMatch(possible_options, matcher.group(1))] = i;
                unitsIndices[i] = getTheClosestMatch(unit_prefixes, matcher.group(3));
            }
        }
        //frame, x, y
        final FloatMatrix data = floatMatrix.getColumns(new int[]{revOptionsIndices[1], revOptionsIndices[2], revOptionsIndices[3]});

        
        final int frames = (int) data.getColumn(0).max(); //1 indexed

        csvTime = System.nanoTime() - csvTime;
        System.out.println("Loading CSV took " + String.format("%.3f", csvTime/1000000000) + " s");
        double processingTime = System.nanoTime();

        int id = 0;

        // Currently doesnt give good results when a situations might arise where:
        // 1 - 2
        // 1 - 3
        // 3 - 4
        // This would stop the chain after 1-2, instead of finding 1-3-4

        int totalCollumns = 2 + (orders * 5);
        FloatMatrix finalPossibilities = new FloatMatrix(0, totalCollumns);
        for(int frame = 1; frame <= frames; frame++){
            IJ.showProgress(frame, frames);
            IJ.showStatus(frame + "/" + frames);

            int[] frameIndicices = data.getColumn(0).eq(frame).findIndices();

            FloatMatrix frameData = data.getRows(frameIndicices);

            FloatMatrix subtractedX = makeSubstractedMatrix(frameData.getColumn(1));
            FloatMatrix subtractedY = makeSubstractedMatrix(frameData.getColumn(2));

            FloatMatrix distances = Distance(subtractedX, subtractedY);
            FloatMatrix angles = atan2(subtractedX, subtractedY, distances);


            int[] correctAngleAndDistance = distances.gt(distRange[0]).and(distances.lt(distRange[1]))
                    .and(angles.gt(angRange[0])).and(angles.lt(angRange[1])).findIndices();


            if(correctAngleAndDistance.length > 1){

                FloatMatrix possibilities = new FloatMatrix(correctAngleAndDistance.length, totalCollumns);
                FloatMatrix intermediateFinalPossibilities = new FloatMatrix(0, totalCollumns);

                for(int i = 0; i < correctAngleAndDistance.length; i++) {
                    int index = correctAngleAndDistance[i];
                    possibilities.putRow(i, extend(new FloatMatrix(1, 10,
                            id++,                                                           //0
                            frame,                                                          //1
                            Math.floorDiv(index, distances.rows),                           //2
                            frameData.get(index / distances.rows, 1),    //3
                            frameData.get(index / distances.rows, 2),    //4
                            index % distances.rows,                                         //5
                            frameData.get(index % distances.rows, 1),    //6
                            frameData.get(index % distances.rows, 2),    //7
                            distances.get(index),                                           //8
                            angles.get(index)                                               //9
                    ),1, totalCollumns));
                }
                FloatMatrix connectedIds = possibilities.getColumn(5);
                List<Integer> toSkip = new ArrayList<>();
                for(int i = 0; i < possibilities.rows; i++){
                    if(!toSkip.contains(i)){
                        FloatMatrix identicalIds = connectedIds.eq(connectedIds.get(i));
                        float identicalIdSum = identicalIds.sum();
                        if(identicalIdSum > 1.0f){

                            float cumX = possibilities.get(i, 3);
                            float cumY = possibilities.get(i, 4);

                            for(int index2 : identicalIds.findIndices()) {
                                if(index2 != i) {
                                    cumX += possibilities.get(index2, 3);
                                    cumY += possibilities.get(index2, 4);
                                    toSkip.add(index2);
                                }
                            }

                            FloatMatrix temp_arr = possibilities.getRow(i).dup();
                            temp_arr.put(3, cumX/identicalIdSum);
                            temp_arr.put(4, cumY/identicalIdSum);

                            intermediateFinalPossibilities = FloatMatrix.concatVertically(intermediateFinalPossibilities, temp_arr);

                        } else {
                            intermediateFinalPossibilities = FloatMatrix.concatVertically(intermediateFinalPossibilities, possibilities.getRow(i));
                        }
                    }
                }
                List<Integer> toKeep = new ArrayList<>();
                for(int row = 0; row < intermediateFinalPossibilities.rows; row ++){
                    float index = intermediateFinalPossibilities.get(row, 2);
                    int curr_row = row;
                    if(intermediateFinalPossibilities.getColumn(5).eq(index).sum() == 0.0f) { //Start of chain
                        toKeep.add(row);
                        for(int order = 2; order < orders; order++){

                            int connected_to = (int) intermediateFinalPossibilities.get(curr_row, 5);
                            int[] connected_indices = intermediateFinalPossibilities.getColumn(2).eq(connected_to).findIndices();
                            if(connected_indices.length > 0) {
                                //Do only one for now
                                curr_row = connected_indices[0];

                                //i.e. 10: id, 11: x, 12: y, 13: distance, 14: angle
                                int[] target_range = new int[]{2 + (order * 4), 3 + (order * 4), 4 + (order * 4), 5 + (order * 4), 6 + (order * 4)};
                                FloatMatrix target_row = intermediateFinalPossibilities.getRow(curr_row);
                                intermediateFinalPossibilities.put(row, target_range, target_row.getColumns(new int[]{5, 6, 7, 8, 9}));

                            } else {
                                break;
                            }
                        }
                        int connected_to = (int) intermediateFinalPossibilities.get(curr_row, 5);
                        int[] connected_indices = intermediateFinalPossibilities.getColumn(2).eq(connected_to).findIndices();
                        if(connected_indices.length > 0) {
                            System.out.println("There seem to be more orders than " + orders);
                        }
                    }
                }

                int[] indices = toKeep.stream().mapToInt(i->i).toArray();
                finalPossibilities = FloatMatrix.concatVertically(finalPossibilities, intermediateFinalPossibilities.getRows(indices));

            } else if(correctAngleAndDistance.length > 0){
                int index = correctAngleAndDistance[0];
                finalPossibilities = FloatMatrix.concatVertically(finalPossibilities, extend(new FloatMatrix(1, 10,
                        id++,
                        frame,
                        Math.floorDiv(index, distances.rows),
                        frameData.get( index/distances.rows, 1),
                        frameData.get( index/distances.rows, 2),
                        index%distances.rows,
                        frameData.get( index%distances.rows, 1),
                        frameData.get( index%distances.rows, 2),
                        distances.get(index),
                        angles.get(index)
                ), 1, totalCollumns));
            }
        }

        //Fix ids
        id = 0;
        List<String> Headers = new ArrayList<>();
        Headers.add("id");
        Headers.add("frame");
        String distanceUnit =  (unit_prefixes[unitsIndices[revOptionsIndices[2]]].equals(unit_prefixes[unitsIndices[revOptionsIndices[3]]]) ? unit_prefixes[unitsIndices[revOptionsIndices[2]]]:("(" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "*" +unit_prefixes[unitsIndices[revOptionsIndices[2]]] + ")^½"));
        for(int i = 0; i < orders; i++){
            Headers.add(i + " index");
            Headers.add(i + " x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]");
            Headers.add(i + " y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]");
            if(i>0){
                Headers.add((i-1) + "-" + i + " distance" + distanceUnit);
                Headers.add((i-1) + "-" + i + "angle");
            }
        }
        //System.out.println(Arrays.toString(Header));
        for(int i = 0; i < finalPossibilities.rows; i++){
            finalPossibilities.put(i, 0, id++);
            //System.out.println(finalPossibilities.getRow(i));
        }

        processingTime = System.nanoTime() - processingTime;
        System.out.println("\nProcessing data took " + String.format("%.3f", processingTime/1000000000) + " s");

        FloatProcessor histData = new FloatProcessor(new float[][]{finalPossibilities.getColumn(8).toArray()});
        ImagePlus dummy = new ImagePlus("WHY", histData);
        HistogramWindow hist = new HistogramWindow("0th-1st", dummy,(int) ((finalPossibilities.getColumn(8).max()-finalPossibilities.getColumn(8).min())/10), distRange[0], distRange[1]);


        FloatMatrix nonZero = finalPossibilities.getColumn(13).get(finalPossibilities.getColumn(13).ne(0.0f).findIndices());
        FloatProcessor histData2 = new FloatProcessor(nonZero.toArray2());
        ImagePlus dummy2 = new ImagePlus("WHY", histData2);
        HistogramWindow hist2 = new HistogramWindow("1st-2nd", dummy2,(int) ((nonZero.max()-nonZero.min())/10), distRange[0], distRange[1]);


        Plot plot = new Plot("Points", "x", "y");

        String[] colors = {"blue", "red", "green", "black"};
        String[] shapes = {"circle", "cross", "box", "diamond"};//  "line", "connected circle", "filled", "bar", "separated bar", "circle", "box", "triangle", "diamond", "cross", "x", "dot", "error bars" or "xerror bars"


        plot.setColor(colors[0]); //0th
        plot.add(shapes[0],toDouble(finalPossibilities.getColumn(3)), toDouble(finalPossibilities.getColumn(4)));
        for(int i = 1; i < orders ; i++){
            plot.setColor(colors[i]);
            plot.add(shapes[i],toDouble(finalPossibilities.getColumn(1 + (i * 5))), toDouble(finalPossibilities.getColumn(2 + (i * 5))));
        } //6, 7//11, 12


        plot.setLegend("0th\t1st\t2nd\t3rd", Plot.AUTO_POSITION);
        plot.show();

        if(saveSCV) SaveCSV(finalPossibilities, Headers, CSV_FILE_NAME);


    }





    public static void main(String[] args) {
        new ImageJ();
        IJ.runPlugIn(sSMLMA.class.getName(), "");

    }
}

