package com.wurgobes.sSMLMAnalyzer;



import ij.*;
import ij.gui.HistogramWindow;
import ij.plugin.FFT;
import ij.process.FloatProcessor;
import ij.process.ShortProcessor;
import net.imagej.ImageJ;
import net.imglib2.*;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import org.jblas.FloatMatrix;
import org.jblas.exceptions.LapackException;
import org.scijava.command.Command;
import org.scijava.plugin.Plugin;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static com.wurgobes.sSMLMAnalyzer.Util.*;
import static com.wurgobes.sSMLMAnalyzer.levenshtein.getTheClosestMatch;

@Plugin(type = Command.class, menuPath = "Plugins>Spectral Analyzer>Analyze Angles and Distance")
public class AngleAnalyzer implements Command {

    private final String[] possible_options = {"id", "frame", "x", "y", "z", "intensity", "offset", "bkgstd", "sigma1", "sigma2", "uncertainty", "detections", "chi"};
    private final String[] unit_prefixes = {"null", "photons", "m", "cm", "mm", "um", "nm", "ang", "pm", "fm"};

    private String filePath = "F:\\ThesisData\\output\\4_grating_drift.csv";

    @Override
    public void run() {

        double csvTime = System.nanoTime();

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


        Pattern pattern = Pattern.compile("(\\w+)( [ (\\[](\\w+)[)\\] ])?");

        for(int i = 0; i < collumns.size(); i ++){
            String header = collumns.get(i);
            Matcher matcher = pattern.matcher(header);
            if(matcher.find()){
                revOptionsIndices[getTheClosestMatch(possible_options, matcher.group(1))] = i;

            }
        }
        //frame, x, y, intensity
        final FloatMatrix data = floatMatrix.getColumns(new int[]{revOptionsIndices[1], revOptionsIndices[2], revOptionsIndices[3], revOptionsIndices[5]});

        final int frames = (int) data.getColumn(0).max(); //1 indexed

        csvTime = System.nanoTime() - csvTime;
        System.out.println("Loading CSV took " + String.format("%.3f", csvTime / 1000000000) + " s");
        double processingTime = System.nanoTime();

        // Implement FFT to guess angle and distance

        System.out.println("Total Frames: " + frames);
        System.out.println("Total Points: " + (int) floatMatrix.getColumn(revOptionsIndices[0]).max());


        int[] dimensions = new int[] {512, 512};
        final int[] reduce = new int[] {(int) Math.ceil((data.getColumn(1).max() + 1)/512), (int) Math.ceil((data.getColumn(2).max() + 1)/512)};

        float[] sum = new float[512*512];

        for (int frame = 1; frame <= 100; frame++) {
            int[] frameIndicices = data.getColumn(0).eq(frame).findIndices();

            FloatMatrix frameData = data.getRows(frameIndicices);

            FloatProcessor temp = getImageFromPoints(frameData.getColumns(new int[]{1, 2, 3}), reduce, dimensions[0], dimensions[1]);

            ImagePlus FFTImage = FFT.forward(new ImagePlus("", temp));

            byte[] vals = (byte[]) FFTImage.getProcessor().getPixelsCopy();
            for(int i = 0; i < vals.length; i++){
                sum[i] +=  vals[i];
            }
        }


        new ImagePlus("Sum", new FloatProcessor(512, 512, sum, null)).show();

    }

    public static void main(String[] args) {
        ImageJ ij = new ImageJ();
        ij.ui().showUI();

        ij.command().run(AngleAnalyzer.class, true);
    }
}
